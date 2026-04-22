// Intent: scan reads, match k-mers, score panels, and write rank outputs for `find`.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "FastqReader.hpp"
#include "PanelFindInternal.hpp"
#include "utils.hpp"

constexpr double kFscoreBeta = 2.0;

// Threshold ratio to apply a panel vs panel specific second pass.
// this happens when two gene panels have very similar sequence content
constexpr double kSecondPassRatio = 0.95;

// Store the pair-specific second-pass scores for the top two panels.
struct PanelPairSecondPassResult {
    double firstScore = 0.0;
    double secondScore = 0.0;
    std::uint64_t firstCovered = 0;
    std::uint64_t secondCovered = 0;
    std::uint64_t firstTotal = 0;
    std::uint64_t secondTotal = 0;
};

// Keep the mutable scan state for one sample while reads are processed.
struct SampleScanState {
    std::vector<std::uint64_t> panelMatchedKmers;
    std::unordered_set<std::uint64_t> matchedLookupKmers;
    std::uint64_t totalReads = 0;
};

// Store per-panel support collected from the sample's matched lookup k-mers.
struct PanelSupportStats {
    std::vector<std::uint64_t> panelCoveredKmers;
    std::vector<double> panelCoveredSpecificityMass;
    double totalMatchedLookupKmers = 0.0;
};

// Hold the score and summary metrics for one candidate panel.
struct PanelCandidateMetrics {
    std::size_t panelId = 0;
    double candidateScore = 0.0;
    double indexUniqueToPanelPercent = 0.0;
    double indexMeanPanelsPerKmer = 0.0;
    double panelKmerExtentPercent = 0.0;
    std::uint64_t matchedKmers = 0;
    std::uint64_t coveredPanelKmers = 0;
};

// Sort candidate panels by score, then coverage, then panel name.
struct PanelCandidateComparator {
    const std::vector<std::string>& panelNames;

    bool operator()(const PanelCandidateMetrics& left, const PanelCandidateMetrics& right) const {
        if (left.candidateScore != right.candidateScore) {
            return left.candidateScore > right.candidateScore;
        }
        if (left.coveredPanelKmers != right.coveredPanelKmers) {
            return left.coveredPanelKmers > right.coveredPanelKmers;
        }
        return panelNames[left.panelId] < panelNames[right.panelId];
    }
};

// Compute the Shannon entropy of one encoded k-mer from its nucleotide composition.
// Higher entropy means the k-mer has a more diverse base composition. Lower entropy means lower complexity.
// The cached table avoids recomputing log2 terms for every k-mer.
static double calculateKmerEntropy(
    std::uint64_t encodedKmer, std::size_t k,
    const std::array<double, 33>& entropyContributionByCount) {
    std::array<std::uint8_t, 4> baseCounts = {0, 0, 0, 0};
    for (std::size_t i = 0; i < k; i = i + 1) {
        std::uint8_t baseCode = static_cast<std::uint8_t>(encodedKmer & 0x3ULL);
        baseCounts[baseCode] = static_cast<std::uint8_t>(baseCounts[baseCode] + 1);
        encodedKmer = encodedKmer >> 2;
    }

    return entropyContributionByCount[baseCounts[0]] + entropyContributionByCount[baseCounts[1]] +
           entropyContributionByCount[baseCounts[2]] + entropyContributionByCount[baseCounts[3]];
}

// Return true when a k-mer is considered low complexity.
static bool isLowComplexityKmer(std::uint64_t encodedKmer,
                                const KmerScanSettings& scanSettings) {
    if (!scanSettings.useLowComplexityFilter) {
        return false;
    }

    double entropy =
        calculateKmerEntropy(encodedKmer, scanSettings.k,
                             scanSettings.entropyContributionByCount);
    return entropy < scanSettings.minKmerEntropy;
}

// Update per-panel counters for one k-mer hit.
static void processEncodedKmerHit(std::uint64_t encodedKmer,
    const PanelIndexStruct& panelIndex, const KmerScanSettings& scanSettings, SampleScanState& state) {

    // Ignore low-complexity k-mers
    if (isLowComplexityKmer(encodedKmer, scanSettings)) {
        return;
    }

    // Skip k-mers that do not appear in any loaded panel index.
    auto hit = panelIndex.globalLookup.find(encodedKmer);
    if (hit == panelIndex.globalLookup.end()) {
        return;
    }

    // Record that this sample matched this lookup k-mer at least once.
    state.matchedLookupKmers.insert(encodedKmer);

    // Count one hit for every panel that contains this k-mer.
    const std::vector<std::size_t>& panelIds = hit->second;
    for (std::size_t j = 0; j < panelIds.size(); j = j + 1) {
        std::size_t panelId = panelIds[j];
        state.panelMatchedKmers[panelId] = state.panelMatchedKmers[panelId] + 1;
    }
}

// Process one contiguous k-mer segment (all k-mers or minimizers only).
static void processKmerSegment(const std::vector<std::uint64_t>& segmentKmers,
    const PanelIndexStruct& panelIndex, const KmerScanSettings& scanSettings, SampleScanState& state) {
    // Nothing to do when the current valid segment has no k-mers.
    if (segmentKmers.empty()) {
        return;
    }

    // Without minimizers, pass every k-mer in the segment to panel matching.
    if (!scanSettings.useMinimizers || segmentKmers.size() == 1) {
        for (std::size_t j = 0; j < segmentKmers.size(); j = j + 1) {
            processEncodedKmerHit(segmentKmers[j], panelIndex, scanSettings, state);
        }
        return;
    }

    // If the segment is shorter than one minimizer window, keep only its smallest k-mer.
    if (segmentKmers.size() < scanSettings.minimizerWindow) {
        std::uint64_t minValue = segmentKmers[0];
        for (std::size_t j = 1; j < segmentKmers.size(); j = j + 1) {
            if (segmentKmers[j] < minValue) {
                minValue = segmentKmers[j];
            }
        }
        processEncodedKmerHit(minValue, panelIndex, scanSettings, state);
        return;
    }

    const std::size_t noPickedIndex = static_cast<std::size_t>(-1);
    std::size_t lastPickedIndex = noPickedIndex;
    // Slide the minimizer window across the segment and pick the smallest k-mer per window.
    for (std::size_t start = 0; start + scanSettings.minimizerWindow <= segmentKmers.size();
         start = start + 1) {
        std::size_t minIndex = start;
        std::uint64_t minValue = segmentKmers[start];
        std::size_t end = start + scanSettings.minimizerWindow;
        for (std::size_t j = start + 1; j < end; j = j + 1) {
            if (segmentKmers[j] < minValue) {
                minValue = segmentKmers[j];
                minIndex = j;
            }
        }

        // Skip duplicate picks when adjacent windows share the same minimizer position.
        if (minIndex == lastPickedIndex) {
            continue;
        }

        processEncodedKmerHit(minValue, panelIndex, scanSettings, state);
        lastPickedIndex = minIndex;
    }
}

// Scan one read sequence and feed valid k-mer segments to panel matching.
static void processSequence(const std::string& sequence, const PanelIndexStruct& panelIndex,
                            const KmerScanSettings& scanSettings, SampleScanState& state) {
    if (sequence.size() < scanSettings.k) {
        return;
    }

    std::uint64_t rollingKmer = 0;
    std::size_t validRun = 0;
    std::vector<std::uint64_t> segmentKmers;
    segmentKmers.reserve(sequence.size() - scanSettings.k + 1);

    for (std::size_t i = 0; i < sequence.size(); i = i + 1) {
        // Convert the current base character to its 2-bit value.
        std::uint8_t baseCode = scanSettings.lookup[static_cast<unsigned char>(sequence[i])];
        // Values above 3 mark non-ACGT bases, which break the current valid k-mer run.
        if (baseCode > 3) {
            // Finish the current valid segment before resetting at this invalid base.
            processKmerSegment(segmentKmers, panelIndex, scanSettings, state);
            segmentKmers.clear();
            rollingKmer = 0;
            validRun = 0;
            continue;
        }

        // Shift in the new base and keep only the last k bases.
        rollingKmer = ((rollingKmer << 2) | static_cast<std::uint64_t>(baseCode)) &
                      scanSettings.rollingMask;
        validRun = validRun + 1;
        // Wait until we have seen enough consecutive valid bases to form a full k-mer.
        if (validRun < scanSettings.k) {
            continue;
        }

        // Store this encoded k-mer as part of the current valid segment.
        segmentKmers.push_back(rollingKmer);
    }

    // Process the final valid segment left at the end of the read.
    processKmerSegment(segmentKmers, panelIndex, scanSettings, state);
}

// Check if a panel id is present in a k-mer panel-id list.
static bool panelSetContains(const std::vector<std::size_t>& panelIds, std::size_t panelId) {
    for (std::size_t i = 0; i < panelIds.size(); i = i + 1) {
        if (panelIds[i] == panelId) {
            return true;
        }
    }
    return false;
}

// Re-score top-2 panels using k-mers present in exactly one of the two.
static PanelPairSecondPassResult runPairSpecificSecondPass(
    std::size_t firstPanelId, std::size_t secondPanelId,
    const std::unordered_map<std::uint64_t, std::vector<std::size_t>>& globalLookup,
    const std::unordered_set<std::uint64_t>& matchedLookupKmers) {
    PanelPairSecondPassResult result;

    for (const auto& kmerPanelsEntry : globalLookup) {
        const std::vector<std::size_t>& panelIds = kmerPanelsEntry.second;
        bool inFirst = panelSetContains(panelIds, firstPanelId);
        bool inSecond = panelSetContains(panelIds, secondPanelId);
        if (inFirst == inSecond) {
            continue;
        }
        if (inFirst) {
            result.firstTotal = result.firstTotal + 1;
        }
        else {
            result.secondTotal = result.secondTotal + 1;
        }
    }

    for (const auto& encodedKmer : matchedLookupKmers) {
        auto hit = globalLookup.find(encodedKmer);
        if (hit == globalLookup.end()) {
            continue;
        }
        const std::vector<std::size_t>& panelIds = hit->second;
        bool inFirst = panelSetContains(panelIds, firstPanelId);
        bool inSecond = panelSetContains(panelIds, secondPanelId);
        if (inFirst == inSecond) {
            continue;
        }
        if (inFirst) {
            result.firstCovered = result.firstCovered + 1;
        }
        else {
            result.secondCovered = result.secondCovered + 1;
        }
    }

    if (result.firstTotal > 0) {
        result.firstScore =
            static_cast<double>(result.firstCovered) / static_cast<double>(result.firstTotal);
    }
    if (result.secondTotal > 0) {
        result.secondScore =
            static_cast<double>(result.secondCovered) / static_cast<double>(result.secondTotal);
    }

    return result;
}

// Stream one FASTQ file and update panel matching state.
static bool processFastq(const std::string& path, bool reverseComplementReads,
                         const PanelIndexStruct& panelIndex,
                         const KmerScanSettings& scanSettings,
                         const ScanLimits& scanLimits,
                         SampleScanState& state,
                         bool ignoreMaxReadsLimit) {
    fastq::FastqReader reader(path);
    if (!reader.isOpen()) {
        std::cerr << "Unable to open FASTQ file: " << path << "\n";
        return false;
    }

    fastq::FastqRecord record;
    while (true) {
        if (!ignoreMaxReadsLimit && state.totalReads >= scanLimits.maxReads) {
            break;
        }

        fastq::FastqReadStatus status = reader.readNext(record);
        if (status == fastq::FastqReadStatus::EndOfFile) {
            break;
        }
        if (status == fastq::FastqReadStatus::OpenError) {
            std::cerr << "Unable to open FASTQ file: " << path << "\n";
            return false;
        }
        if (status == fastq::FastqReadStatus::TruncatedRecord) {
            std::cerr << "Truncated FASTQ in file: " << path << "\n";
            return false;
        }

        if (reverseComplementReads) {
            record.sequence = utils::reverseComplement(record.sequence);
        }
        state.totalReads = state.totalReads + 1;
        processSequence(record.sequence, panelIndex, scanSettings, state);
    }

    return true;
}

static bool scanSampleReads(const SampleFqEntry& sampleFq,
                            const PanelIndexStruct& panelIndex,
                            const KmerScanSettings& scanSettings,
                            const ScanLimits& scanLimits,
                            SampleScanState& state) {
    if (!sampleFq.fq1.empty()) {
        std::uint64_t readsBefore = state.totalReads;
        bool ok = processFastq(sampleFq.fq1, false, panelIndex, scanSettings, scanLimits, state,
                               false);
        if (!ok) {
            return false;
        }
        if (sampleFq.hasFq2) {
            std::cout << " INFO: fq1_reads=" << (state.totalReads - readsBefore) << std::endl;
        }
    }

    if (sampleFq.hasFq2) {
        if (state.totalReads >= scanLimits.maxReads && !scanLimits.forcePaired) {
            std::cout << " INFO: fq2_skipped=max_reads_reached" << std::endl;
        }
        else {
            bool ignoreMaxReadsForFq2 =
                scanLimits.forcePaired && (state.totalReads >= scanLimits.maxReads);
            std::uint64_t readsBefore = state.totalReads;
            bool ok = processFastq(sampleFq.fq2, true, panelIndex, scanSettings, scanLimits,
                                   state, ignoreMaxReadsForFq2);
            if (!ok) {
                return false;
            }
            std::cout << " INFO: fq2_reads=" << (state.totalReads - readsBefore) << std::endl;
        }
    }

    return true;
}

static PanelSupportStats buildPanelSupportStats(const PanelIndexStruct& panelIndex,
                                                const SampleScanState& state) {
    PanelSupportStats supportStats;
    supportStats.panelCoveredKmers.assign(panelIndex.panelNames.size(), 0);
    supportStats.panelCoveredSpecificityMass.assign(panelIndex.panelNames.size(), 0.0);
    supportStats.totalMatchedLookupKmers =
        static_cast<double>(state.matchedLookupKmers.size());

    for (const auto& encodedKmer : state.matchedLookupKmers) {
        auto hit = panelIndex.globalLookup.find(encodedKmer);
        if (hit == panelIndex.globalLookup.end()) {
            continue;
        }

        const std::vector<std::size_t>& panelIds = hit->second;
        double specificityWeight = 1.0 / static_cast<double>(panelIds.size());
        for (std::size_t j = 0; j < panelIds.size(); j = j + 1) {
            std::size_t panelId = panelIds[j];
            supportStats.panelCoveredKmers[panelId] =
                supportStats.panelCoveredKmers[panelId] + 1;
            supportStats.panelCoveredSpecificityMass[panelId] =
                supportStats.panelCoveredSpecificityMass[panelId] + specificityWeight;
        }
    }

    return supportStats;
}

static std::vector<PanelCandidateMetrics> buildPanelCandidates(
    const PanelIndexStruct& panelIndex, const SampleScanState& state,
    const PanelSupportStats& supportStats) {
    std::vector<PanelCandidateMetrics> panelCandidates;
    panelCandidates.reserve(panelIndex.panelNames.size());

    for (std::size_t i = 0; i < panelIndex.panelNames.size(); i = i + 1) {
        PanelCandidateMetrics metrics;
        metrics.panelId = i;
        metrics.matchedKmers = state.panelMatchedKmers[i];
        metrics.coveredPanelKmers = supportStats.panelCoveredKmers[i];

        // Fraction of this panel's indexed k-mers that were covered by the sample.
        double panelCoverage = 0.0;
        if (panelIndex.panelIndexUniqueKmers[i] > 0) {
            panelCoverage = static_cast<double>(metrics.coveredPanelKmers) /
                            static_cast<double>(panelIndex.panelIndexUniqueKmers[i]);
        }

        double panelSpecificSupport = 0.0;
        if (supportStats.totalMatchedLookupKmers > 0.0) {
            panelSpecificSupport =
                supportStats.panelCoveredSpecificityMass[i] /
                supportStats.totalMatchedLookupKmers;
        }

        if (panelIndex.panelIndexUniqueKmers[i] > 0) {
            metrics.panelKmerExtentPercent =
                (100.0 * static_cast<double>(metrics.coveredPanelKmers)) /
                static_cast<double>(panelIndex.panelIndexUniqueKmers[i]);
        }

        const double betaSquared = kFscoreBeta * kFscoreBeta;
        const double denominator = (betaSquared * panelSpecificSupport) + panelCoverage;
        if (denominator > 0.0) {
            metrics.candidateScore =
                100.0 * (1.0 + betaSquared) * panelSpecificSupport * panelCoverage / denominator;
        }
        else {
            metrics.candidateScore = 0.0;
        }

        if (panelIndex.panelIndexUniqueKmers[i] > 0) {
            metrics.indexUniqueToPanelPercent =
                (100.0 * static_cast<double>(panelIndex.panelIndexUniqueToPanel[i])) /
                static_cast<double>(panelIndex.panelIndexUniqueKmers[i]);
            metrics.indexMeanPanelsPerKmer =
                static_cast<double>(panelIndex.panelIndexPanelsPerKmerSum[i]) /
                static_cast<double>(panelIndex.panelIndexUniqueKmers[i]);
        }

        panelCandidates.push_back(metrics);
    }

    std::sort(panelCandidates.begin(), panelCandidates.end(),
              PanelCandidateComparator{panelIndex.panelNames});
    return panelCandidates;
}

static bool scoreTopPanels(const std::string& sampleName,
                           const PanelIndexStruct& panelIndex,
                           const SampleScanState& state,
                           std::vector<PanelCandidateMetrics>& panelCandidates) {
    if (panelCandidates.size() <= 1 || panelCandidates[0].candidateScore <= 0.0 ||
        panelCandidates[1].candidateScore <= 0.0 ||
        panelCandidates[1].candidateScore <
            (panelCandidates[0].candidateScore * kSecondPassRatio)) {
        return false;
    }

    std::size_t firstPanelId = panelCandidates[0].panelId;
    std::size_t secondPanelId = panelCandidates[1].panelId;
    PanelPairSecondPassResult panelPairSecondPassResult =
        runPairSpecificSecondPass(firstPanelId, secondPanelId,
                                  panelIndex.globalLookup, state.matchedLookupKmers);

    if (panelPairSecondPassResult.firstTotal == 0) {
        std::cout << " WARNING: \t" << sampleName
                  << " second_pass=pair_specific"
                  << " panel_a=" << panelIndex.panelNames[firstPanelId]
                  << " no_pair_specific_kmers_found=true"
                  << std::endl;
    }
    if (panelPairSecondPassResult.secondTotal == 0) {
        std::cout << " WARNING: \t" << sampleName
                  << " second_pass=pair_specific"
                  << " panel_b=" << panelIndex.panelNames[secondPanelId]
                  << " no_pair_specific_kmers_found=true"
                  << std::endl;
    }

    if (panelPairSecondPassResult.firstTotal == 0 && panelPairSecondPassResult.secondTotal == 0) {
        return false;
    }

    panelCandidates[0].candidateScore = panelPairSecondPassResult.firstScore;
    panelCandidates[1].candidateScore = panelPairSecondPassResult.secondScore;
    if (panelCandidates[1].candidateScore > panelCandidates[0].candidateScore) {
        std::swap(panelCandidates[0], panelCandidates[1]);
    }

    std::cout << " INFO: \t" << sampleName
              << " second_pass=pair_specific"
              << " panel_a=" << panelIndex.panelNames[panelCandidates[0].panelId]
              << " panel_a_specific_score=" << panelCandidates[0].candidateScore
              << " panel_b=" << panelIndex.panelNames[panelCandidates[1].panelId]
              << " panel_b_specific_score=" << panelCandidates[1].candidateScore
              << std::endl;
    return true;
}

static void storeRankResults(const std::string& sampleName, bool minReadsReached,
                             bool maxReadsReached, const SampleScanState& state,
                             const PanelIndexStruct& panelIndex,
                             const std::vector<PanelCandidateMetrics>& panelCandidates,
                             double minPanelScore, bool usedTopPanelRescore,
                             PanelRankRecord& rankRecord) {
    rankRecord.sample = sampleName;
    rankRecord.totalReads = state.totalReads;
    rankRecord.minReadsReached = minReadsReached;
    rankRecord.maxReadsReached = maxReadsReached;

    if (panelCandidates.empty()) {
        std::cout << " INFO: \t" << sampleName
                  << " best_panel=none"
                  << " scanned_reads=" << rankRecord.totalReads << std::endl;
        return;
    }

    const PanelCandidateMetrics& best = panelCandidates[0];
    std::string nextPanel = "none";
    double nextScore = 0.0;
    if (panelCandidates.size() > 1) {
        nextPanel = panelIndex.panelNames[panelCandidates[1].panelId];
        nextScore = panelCandidates[1].candidateScore;
    }

    rankRecord.bestScore = best.candidateScore;
    rankRecord.scoreMarginVsNext = best.candidateScore - nextScore;
    rankRecord.bestMatchedKmers = best.matchedKmers;
    rankRecord.bestCoveredPanelKmers = best.coveredPanelKmers;
    rankRecord.bestPanelKmerExtentPercent = best.panelKmerExtentPercent;
    rankRecord.bestIndexUniqueToPanelPercent = best.indexUniqueToPanelPercent;
    rankRecord.bestIndexMeanPanelsPerKmer = best.indexMeanPanelsPerKmer;

    // Keep low-confidence evidence in the TSV, but do not call it a valid panel.
    bool scoreIsTooLow = best.candidateScore <= 0.0;
    if (minPanelScore > 0.0 && best.candidateScore < minPanelScore) {
        scoreIsTooLow = true;
    }

    if (scoreIsTooLow) {
        std::cout << " INFO: \t" << sampleName
                  << " best_panel=none"
                  << " best_candidate=" << panelIndex.panelNames[best.panelId]
                  << " best_score=" << rankRecord.bestScore
                  << " min_score=" << minPanelScore
                  << " scanned_reads=" << rankRecord.totalReads
                  << std::endl;
        return;
    }

    rankRecord.bestPanel = panelIndex.panelNames[best.panelId];

    std::cout << " INFO: \t" << sampleName
              << " best_panel=" << rankRecord.bestPanel
              << " best_score=" << rankRecord.bestScore
              << " second_panel=" << nextPanel
              << " second_score=" << nextScore
              << " scoring_pass=" << (usedTopPanelRescore ? "pair_specific_second_pass" : "primary")
              << std::endl;
}

// Analyze one sample entry and compute panel ranking metrics.
static bool scanSample(const SampleFqEntry& sampleFq, std::size_t sampleNum,
                         std::size_t totalSamples, const PanelIndexStruct& panelIndex,
                         const KmerScanSettings& scanSettings,
                         const ScanLimits& scanLimits,
                         PanelRankRecord& rankRecord) {
    SampleScanState state;
    state.panelMatchedKmers.assign(panelIndex.panelNames.size(), 0);

    std::string sampleName = sampleFq.sampleKey.empty() ? sampleFq.label : sampleFq.sampleKey;
    std::cout << " INFO: (" << sampleNum << "/" << totalSamples << ") " << sampleName
              << std::endl;
    if (sampleFq.hasFq2) {
        std::cout << " INFO:  Scanning " << sampleFq.fq1 << std::endl;
        std::cout << " INFO:  Scanning " << sampleFq.fq2 << std::endl;
    }
    else if (!sampleFq.fq1.empty()) {
        std::cout << " INFO:  Scanning " << sampleFq.fq1 << std::endl;
    }
    else if (!sampleFq.fq2.empty()) {
        std::cout << " INFO:  Scanning " << sampleFq.fq2 << std::endl;
    }
    bool ok = scanSampleReads(sampleFq, panelIndex, scanSettings, scanLimits, state);
    if (!ok) {
        return false;
    }

    bool minReadsReached = false;
    if (state.totalReads >= scanLimits.minReads) {
        minReadsReached = true;
    }
    else {
        minReadsReached = false;
    }

    bool maxReadsReached = false;
    if (state.totalReads >= scanLimits.maxReads) {
        maxReadsReached = true;
    }
    else {
        maxReadsReached = false;
    }
    if (!minReadsReached) {
        std::cerr << " INFO: Warning: scanned reads (" << state.totalReads
                  << ") below min_reads (" << scanLimits.minReads << ")\n";
    }

    // Summarize how much support each panel gets from the matched lookup k-mers
    PanelSupportStats supportStats = buildPanelSupportStats(panelIndex, state);

    // Turn that support into sortable panel scores and coverage metrics
    std::vector<PanelCandidateMetrics> panelCandidates =
        buildPanelCandidates(panelIndex, state, supportStats);

    // Score the top panels when the first-pass result is too close (score of 0.95)
    bool rescoreTopPanels = scoreTopPanels(sampleName, panelIndex, state, panelCandidates);

    // Save the best-panel result and print the final summary for a given sample
    storeRankResults(sampleName, minReadsReached, maxReadsReached, state, panelIndex,
                     panelCandidates, scanLimits.minPanelScore, rescoreTopPanels, rankRecord);
    return true;
}

static void buildFastqPairOrder(
    const std::vector<SampleFqEntry>& sampleFqs,
    std::unordered_map<std::string, std::size_t>& fastqPairOrderByKey,
    std::size_t& totalSamples) {
    fastqPairOrderByKey.clear();
    fastqPairOrderByKey.reserve(sampleFqs.size());
    totalSamples = 0;

    for (std::size_t i = 0; i < sampleFqs.size(); i = i + 1) {
        const std::string& key = sampleFqs[i].sampleKey;
        if (fastqPairOrderByKey.find(key) == fastqPairOrderByKey.end()) {
            totalSamples = totalSamples + 1;
            fastqPairOrderByKey[key] = totalSamples;
        }
    }
}

std::array<double, 33> buildEntropyContributionByCount(std::size_t k, bool useLowComplexityFilter) {
    std::array<double, 33> entropyContributionByCount = {};
    if (!useLowComplexityFilter) {
        return entropyContributionByCount;
    }

    double kAsDouble = static_cast<double>(k);
    for (std::size_t count = 1; count <= k; count = count + 1) {
        double probability = static_cast<double>(count) / kAsDouble;
        entropyContributionByCount[count] = -probability * std::log2(probability);
    }
    return entropyContributionByCount;
}

bool scanAllSamples(const std::vector<SampleFqEntry>& sampleFqs,
                    const PanelIndexStruct& panelIndex,
                    const KmerScanSettings& scanSettings,
                    const ScanLimits& scanLimits,
                    std::vector<PanelRankRecord>& panelRankRecords) {
    std::unordered_map<std::string, std::size_t> fastqPairOrderByKey;
    std::size_t totalSamples = 0;
    buildFastqPairOrder(sampleFqs, fastqPairOrderByKey, totalSamples);

    panelRankRecords.clear();
    panelRankRecords.reserve(sampleFqs.size());
    std::unordered_set<std::string> samplesAtMaxReads;

    // Scan each sample, score panels, and keep the best-panel summary.
    for (std::size_t sampleIndex = 0; sampleIndex < sampleFqs.size(); sampleIndex = sampleIndex + 1) {
        const SampleFqEntry& sampleFq = sampleFqs[sampleIndex];

        // Skip later FASTQ entries for this sample once the read cap was already reached.
        if (!scanLimits.forcePaired &&
            samplesAtMaxReads.find(sampleFq.sampleKey) != samplesAtMaxReads.end()) {
            if (sampleFq.hasFq2) {
                std::cout << " INFO:  Skipping FASTQ pair: " << sampleFq.fq1 << " "
                          << sampleFq.fq2 << " (reached max_reads)" << std::endl;
            }
            else if (!sampleFq.fq1.empty()) {
                std::cout << " INFO:  Skipping FASTQ: " << sampleFq.fq1
                          << " (reached max_reads)" << std::endl;
            }
            else {
                std::cout << " INFO:  Skipping FASTQ: " << sampleFq.fq2
                          << " (reached max_reads)" << std::endl;
            }
            continue;
        }

        PanelRankRecord rankRecord;
        std::size_t sampleNum = fastqPairOrderByKey[sampleFq.sampleKey];
        bool ok = scanSample(sampleFq, sampleNum, totalSamples, panelIndex, scanSettings, scanLimits, rankRecord);
        if (!ok) {
            return false;
        }

        // Keep one output record per scanned sample entry.
        panelRankRecords.push_back(rankRecord);

        // Mark this sample so later entries can be skipped once max_reads was reached.
        if (!scanLimits.forcePaired && rankRecord.maxReadsReached) {
            samplesAtMaxReads.insert(sampleFq.sampleKey);
        }
    }

    return true;
}

bool writePanelRanksFile(const std::string& panelRanksPath, const std::vector<PanelRankRecord>& panelRankRecords) {
    std::ofstream panelRanksOut(panelRanksPath);
    if (!panelRanksOut.is_open()) {
        std::cerr << "Unable to write panel ranks file: " << panelRanksPath << "\n";
        return false;
    }

    panelRanksOut
        << "sample\tscanned_reads\tbest_panel\tbest_score\tscore_margin_vs_next\tbest_panel_covered_kmers\tbest_panel_covered_kmers_pct\n";
    for (std::size_t i = 0; i < panelRankRecords.size(); i = i + 1) {
        const PanelRankRecord& record = panelRankRecords[i];
        panelRanksOut << record.sample << "\t" << record.totalReads << "\t" << record.bestPanel
                      << "\t" << record.bestScore << "\t" << record.scoreMarginVsNext << "\t"
                      << record.bestCoveredPanelKmers << "\t"
                      << record.bestPanelKmerExtentPercent << "\n";
    }

    if (!panelRanksOut.good()) {
        std::cerr << "Error writing panel ranks file: " << panelRanksPath << "\n";
        return false;
    }

    return true;
}
