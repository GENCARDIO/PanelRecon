#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "FastqReader.hpp"
#include "parser.hpp"
#include "utils.hpp"



constexpr double kFscoreBeta = 2.0;

// Threshold to apply a pair-specific second pass
constexpr double kSecondPassScoreRatioTrigger = 0.95;

struct PairSpecificPassScores {
    double firstScore = 0.0;
    double secondScore = 0.0;
    std::uint64_t firstCovered = 0;
    std::uint64_t secondCovered = 0;
    std::uint64_t firstTotal = 0;
    std::uint64_t secondTotal = 0;
};

char complementBase(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
        case 'U':
        case 'u':
            return 'A';
        default:
            return 'N';
    }
}

void reverseComplementInPlace(std::string& sequence) {
    if (sequence.empty()) {
        return;
    }

    std::size_t left = 0;
    std::size_t right = sequence.size() - 1;
    while (left < right) {
        char leftComplement = complementBase(sequence[left]);
        char rightComplement = complementBase(sequence[right]);
        sequence[left] = rightComplement;
        sequence[right] = leftComplement;
        left = left + 1;
        right = right - 1;
    }

    if (left == right) {
        sequence[left] = complementBase(sequence[left]);
    }
}

struct SampleFqEntry {
    std::string label;
    std::string sampleKey;
    std::string fq1;
    std::string fq2;
    bool hasFq2 = false;
};

bool endsWith(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string deriveSampleKeyFromFastqPath(const std::string& fastqPath) {
    std::string name = std::filesystem::path(fastqPath).filename().string();
    if (endsWith(name, ".gz")) {
        name = name.substr(0, name.size() - 3);
    }
    if (endsWith(name, ".fastq")) {
        name = name.substr(0, name.size() - 6);
    } 
    else if (endsWith(name, ".fq")) {
        name = name.substr(0, name.size() - 3);
    }

    bool strippedIlluminaReadToken = false;
    std::size_t readMarkerPos = name.rfind("_R1_");
    if (readMarkerPos == std::string::npos) {
        readMarkerPos = name.rfind("_R2_");
    }
    if (readMarkerPos != std::string::npos) {
        name = name.substr(0, readMarkerPos);
        strippedIlluminaReadToken = true;
    } 
    else if (endsWith(name, "_R1") || endsWith(name, "_R2")) {
        name = name.substr(0, name.size() - 3);
        strippedIlluminaReadToken = true;
    }

    bool strippedIlluminaLaneToken = false;
    if (name.size() >= 5) {
        std::size_t lanePos = name.size() - 5;
        if (name[lanePos] == '_' && name[lanePos + 1] == 'L' &&
            std::isdigit(static_cast<unsigned char>(name[lanePos + 2])) != 0 &&
            std::isdigit(static_cast<unsigned char>(name[lanePos + 3])) != 0 &&
            std::isdigit(static_cast<unsigned char>(name[lanePos + 4])) != 0) {
            name = name.substr(0, lanePos);
            strippedIlluminaLaneToken = true;
        }
    }

    if ((strippedIlluminaReadToken || strippedIlluminaLaneToken) && name.size() > 2) {
        std::size_t sampleNumberPos = name.rfind("_S");
        if (sampleNumberPos != std::string::npos && sampleNumberPos + 2 < name.size()) {
            bool onlyDigits = true;
            for (std::size_t i = sampleNumberPos + 2; i < name.size(); i = i + 1) {
                if (std::isdigit(static_cast<unsigned char>(name[i])) == 0) {
                    onlyDigits = false;
                    break;
                }
            }
            if (onlyDigits) {
                name = name.substr(0, sampleNumberPos);
            }
        }
    }

    return name;
}

bool loadFastqListFile(const std::string& listPath, std::vector<SampleFqEntry>& entries) {
    namespace fs = std::filesystem;
    entries.clear();

    std::ifstream in(listPath);
    if (!in.is_open()) {
        std::cerr << "Unable to open FASTQ list file: " << listPath << "\n";
        return false;
    }

    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(in, line)) {
        lineNumber = lineNumber + 1;
        std::string trimmed = utils::trimWhitespace(line);
        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }

        std::istringstream lineStream(trimmed);
        std::vector<std::string> tokens;
        std::string token;
        while (lineStream >> token) {
            tokens.push_back(token);
        }

        if (tokens.empty()) {
            continue;
        }
        if (tokens.size() > 2) {
            std::cerr << "Invalid FASTQ list line " << lineNumber
                      << ": expected one or two paths, got " << tokens.size() << "\n";
            return false;
        }

        SampleFqEntry entry;
        entry.fq1 = utils::resolveListPathEntry(listPath, tokens[0]);
        if (tokens.size() == 2) {
            entry.fq2 = utils::resolveListPathEntry(listPath, tokens[1]);
            entry.hasFq2 = true;
            entry.label = fs::path(entry.fq1).stem().string() + "+" + fs::path(entry.fq2).stem().string();
        } else {
            entry.label = fs::path(entry.fq1).stem().string();
        }
        entry.sampleKey = deriveSampleKeyFromFastqPath(entry.fq1);
        if (entry.sampleKey.empty()) {
            entry.sampleKey = entry.label;
        }
        if (entry.label.empty()) {
            entry.label = "line_" + std::to_string(lineNumber);
        }
        entries.push_back(entry);
    }

    if (entries.empty()) {
        std::cerr << "No FASTQ paths found in list file: " << listPath << "\n";
        return false;
    }

    std::sort(entries.begin(), entries.end(), [](const SampleFqEntry& left,
                                                 const SampleFqEntry& right) {
        if (left.sampleKey != right.sampleKey) {
            return left.sampleKey < right.sampleKey;
        }
        if (left.fq1 != right.fq1) {
            return left.fq1 < right.fq1;
        }
        if (left.fq2 != right.fq2) {
            return left.fq2 < right.fq2;
        }
        if (left.hasFq2 != right.hasFq2) {
            return left.hasFq2;
        }
        return left.label < right.label;
    });

    return true;
}

std::size_t countInputFastqFiles(const std::vector<SampleFqEntry>& sampleFqs) {
    std::size_t totalFastqFiles = 0;
    for (std::size_t i = 0; i < sampleFqs.size(); i = i + 1) {
        if (!sampleFqs[i].fq1.empty()) {
            totalFastqFiles = totalFastqFiles + 1;
        }
        if (!sampleFqs[i].fq2.empty()) {
            totalFastqFiles = totalFastqFiles + 1;
        }
    }
    return totalFastqFiles;
}

struct PanelRankRecord {
    std::string sample;
    std::string bestPanel = "none";
    double bestScore = 0.0;
    double scoreMarginVsNext = 0.0;
    std::uint64_t totalReads = 0;
    bool minReadsReached = false;
    bool maxReadsReached = false;
    std::uint64_t bestMatchedKmers = 0;
    std::uint64_t bestReadsWithMatch = 0;
    double bestReadSupportPercent = 0.0;
    std::uint64_t bestCoveredIndexKmers = 0;
    double bestPanelKmerExtentPercent = 0.0;
    double bestIndexUniqueToPanelPercent = 0.0;
    double bestIndexMeanPanelsPerKmer = 0.0;
};

struct PanelCandidateMetrics {
    std::size_t panelId = 0;
    double candidateScore = 0.0;
    double readSupportPercent = 0.0;
    double indexUniqueToPanelPercent = 0.0;
    double indexMeanPanelsPerKmer = 0.0;
    double panelKmerExtentPercent = 0.0;
    std::uint64_t matchedKmers = 0;
    std::uint64_t readsWithMatch = 0;
    std::uint64_t coveredIndexKmers = 0;
};

struct PanelCandidateComparator {
    const std::vector<std::string>& panelNames;

    bool operator()(const PanelCandidateMetrics& left, const PanelCandidateMetrics& right) const {
        if (left.candidateScore != right.candidateScore) {
            return left.candidateScore > right.candidateScore;
        }
        if (left.coveredIndexKmers != right.coveredIndexKmers) {
            return left.coveredIndexKmers > right.coveredIndexKmers;
        }
        if (left.readsWithMatch != right.readsWithMatch) {
            return left.readsWithMatch > right.readsWithMatch;
        }
        return panelNames[left.panelId] < panelNames[right.panelId];
    }
};

struct FindSampleContext {
    const std::vector<std::string>& panelNames;
    const std::unordered_map<std::uint64_t, std::vector<std::size_t>>& globalLookup;
    const std::vector<std::uint64_t>& panelIndexUniqueKmers;
    const std::vector<std::uint64_t>& panelIndexUniqueToPanel;
    const std::vector<std::uint64_t>& panelIndexPanelMultiplicitySum;
    const std::array<std::uint8_t, 256>& lookup;
    const std::string& indexDir;
    std::size_t k = 0;
    std::uint64_t rollingMask = 0;
    bool useMinimizers = false;
    std::size_t minimizerWindow = 1;
    std::uint64_t minReads = 0;
    std::uint64_t maxReads = 0;
    bool useLowComplexityFilter = false;
    double minKmerEntropy = 0.0;
    std::array<double, 33> entropyContributionByCount = {};
    bool forcePaired = false;
};

struct FindSampleState {
    std::vector<std::uint64_t> panelMatchedKmers;
    std::vector<std::uint64_t> panelReadsWithMatch;
    std::vector<std::uint32_t> panelSeenEpoch;
    std::unordered_set<std::uint64_t> matchedLookupKmers;
    std::uint64_t totalReads = 0;
    std::uint32_t readEpoch = 0;
};

void advanceReadEpoch(FindSampleState& state) {
    if (state.readEpoch == UINT32_MAX) {
        std::fill(state.panelSeenEpoch.begin(), state.panelSeenEpoch.end(), 0);
        state.readEpoch = 1;
    } else {
        state.readEpoch = state.readEpoch + 1;
    }
}

double calculateKmerEntropy(std::uint64_t encodedKmer, std::size_t k, const std::array<double, 33>& entropyContributionByCount) {
    std::array<std::uint8_t, 4> baseCounts = {0, 0, 0, 0};
    for (std::size_t i = 0; i < k; i = i + 1) {
        std::uint8_t code = static_cast<std::uint8_t>(encodedKmer & 0x3ULL);
        baseCounts[code] = static_cast<std::uint8_t>(baseCounts[code] + 1);
        encodedKmer = encodedKmer >> 2;
    }

    return entropyContributionByCount[baseCounts[0]] + entropyContributionByCount[baseCounts[1]] +
           entropyContributionByCount[baseCounts[2]] + entropyContributionByCount[baseCounts[3]];
}

bool shouldSkipKmerAsLowComplexity(std::uint64_t encodedKmer, const FindSampleContext& context) {
    if (!context.useLowComplexityFilter) {
        return false;
    }
    double entropy =
        calculateKmerEntropy(encodedKmer, context.k, context.entropyContributionByCount);
    return entropy < context.minKmerEntropy;
}

void processEncodedKmerHit(std::uint64_t encodedKmer, const FindSampleContext& context,
                           FindSampleState& state) {
    if (shouldSkipKmerAsLowComplexity(encodedKmer, context)) {
        return;
    }

    auto hit = context.globalLookup.find(encodedKmer);
    if (hit == context.globalLookup.end()) {
        return;
    }

    state.matchedLookupKmers.insert(encodedKmer);

    const std::vector<std::size_t>& panelIds = hit->second;
    for (std::size_t j = 0; j < panelIds.size(); j = j + 1) {
        std::size_t panelId = panelIds[j];
        state.panelMatchedKmers[panelId] = state.panelMatchedKmers[panelId] + 1;
        if (state.panelSeenEpoch[panelId] != state.readEpoch) {
            state.panelSeenEpoch[panelId] = state.readEpoch;
            state.panelReadsWithMatch[panelId] = state.panelReadsWithMatch[panelId] + 1;
        }
    }
}

void processKmerSegment(const std::vector<std::uint64_t>& segmentKmers, const FindSampleContext& context,
                        FindSampleState& state) {
    if (segmentKmers.empty()) {
        return;
    }

    if (!context.useMinimizers || segmentKmers.size() == 1) {
        for (std::size_t j = 0; j < segmentKmers.size(); j = j + 1) {
            processEncodedKmerHit(segmentKmers[j], context, state);
        }
        return;
    }

    if (segmentKmers.size() < context.minimizerWindow) {
        std::uint64_t minValue = segmentKmers[0];
        for (std::size_t j = 1; j < segmentKmers.size(); j = j + 1) {
            if (segmentKmers[j] < minValue) {
                minValue = segmentKmers[j];
            }
        }
        processEncodedKmerHit(minValue, context, state);
        return;
    }

    const std::size_t noPickedIndex = static_cast<std::size_t>(-1);
    std::size_t lastPickedIndex = noPickedIndex;
    for (std::size_t start = 0; start + context.minimizerWindow <= segmentKmers.size();
         start = start + 1) {
        std::size_t minIndex = start;
        std::uint64_t minValue = segmentKmers[start];
        std::size_t end = start + context.minimizerWindow;
        for (std::size_t j = start + 1; j < end; j = j + 1) {
            if (segmentKmers[j] < minValue) {
                minValue = segmentKmers[j];
                minIndex = j;
            }
        }

        if (minIndex == lastPickedIndex) {
            continue;
        }

        processEncodedKmerHit(minValue, context, state);
        lastPickedIndex = minIndex;
    }
}

void processSequenceForPanels(const std::string& sequence, const FindSampleContext& context,
                              FindSampleState& state) {
    if (sequence.size() < context.k) {
        return;
    }

    advanceReadEpoch(state);

    std::uint64_t rollingKmer = 0;
    std::size_t validRun = 0;
    std::vector<std::uint64_t> segmentKmers;
    segmentKmers.reserve(sequence.size() - context.k + 1);

    for (std::size_t i = 0; i < sequence.size(); i = i + 1) {
        std::uint8_t code = context.lookup[static_cast<unsigned char>(sequence[i])];
        if (code > 3) {
            processKmerSegment(segmentKmers, context, state);
            segmentKmers.clear();
            rollingKmer = 0;
            validRun = 0;
            continue;
        }

        rollingKmer = ((rollingKmer << 2) | static_cast<std::uint64_t>(code)) & context.rollingMask;
        validRun = validRun + 1;
        if (validRun < context.k) {
            continue;
        }

        segmentKmers.push_back(rollingKmer);
    }

    processKmerSegment(segmentKmers, context, state);
}

bool panelSetContains(const std::vector<std::size_t>& panelIds, std::size_t panelId) {
    for (std::size_t i = 0; i < panelIds.size(); i = i + 1) {
        if (panelIds[i] == panelId) {
            return true;
        }
    }
    return false;
}

PairSpecificPassScores runPairSpecificSecondPass(
    std::size_t firstPanelId, std::size_t secondPanelId,
    const std::unordered_map<std::uint64_t, std::vector<std::size_t>>& globalLookup,
    const std::unordered_set<std::uint64_t>& matchedLookupKmers) {
    PairSpecificPassScores result;

    for (const auto& entry : globalLookup) {
        const std::vector<std::size_t>& panelIds = entry.second;
        bool inFirst = panelSetContains(panelIds, firstPanelId);
        bool inSecond = panelSetContains(panelIds, secondPanelId);
        if (inFirst == inSecond) {
            continue;
        }
        if (inFirst) {
            result.firstTotal = result.firstTotal + 1;
        } else {
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
        } else {
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

bool processFastqForPanels(const std::string& path, bool reverseComplementReads,
                           const FindSampleContext& context, FindSampleState& state,
                           bool ignoreMaxReadsLimit) {
    fastq::FastqReader reader(path);
    if (!reader.isOpen()) {
        std::cerr << "Unable to open FASTQ file: " << path << "\n";
        return false;
    }

    fastq::FastqRecord record;
    // Simple loop for beginners: one FASTQ record at a time.
    while (true) {
        if (!ignoreMaxReadsLimit && state.totalReads >= context.maxReads) {
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
            reverseComplementInPlace(record.sequence);
        }
        state.totalReads = state.totalReads + 1;
        processSequenceForPanels(record.sequence, context, state);
    }

    return true;
}

bool scanOneSampleFq(const SampleFqEntry& sampleFq, std::size_t sampleOrdinal,
                     std::size_t totalSamples,
                     const FindSampleContext& context, PanelRankRecord& rankRecord) {
    FindSampleState state;
    state.panelMatchedKmers.assign(context.panelNames.size(), 0);
    state.panelReadsWithMatch.assign(context.panelNames.size(), 0);
    state.panelSeenEpoch.assign(context.panelNames.size(), 0);

    std::string sampleName = sampleFq.sampleKey.empty() ? sampleFq.label : sampleFq.sampleKey;
    std::cout << " INFO: (" << sampleOrdinal << "/" << totalSamples << ") " << sampleName
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

    if (!sampleFq.fq1.empty()) {
        std::uint64_t readsBefore = state.totalReads;
        bool ok = processFastqForPanels(sampleFq.fq1, false, context, state, false);
        if (!ok) {
            return false;
        }
        if (sampleFq.hasFq2) {
            std::cout << " INFO: fq1_reads=" << (state.totalReads - readsBefore) << std::endl;
        }
    }

    if (sampleFq.hasFq2) {
        if (state.totalReads >= context.maxReads && !context.forcePaired) {
            std::cout << " INFO: fq2_skipped=max_reads_reached" << std::endl;
        } else {
            bool ignoreMaxReadsForFq2 = context.forcePaired && (state.totalReads >= context.maxReads);
            std::uint64_t readsBefore = state.totalReads;
            bool ok = processFastqForPanels(sampleFq.fq2, true, context, state, ignoreMaxReadsForFq2);
            if (!ok) {
                return false;
            }
            std::cout << " INFO: fq2_reads=" << (state.totalReads - readsBefore) << std::endl;
        }
    }

    bool minReadsReached = state.totalReads >= context.minReads;
    bool maxReadsReached = state.totalReads >= context.maxReads;
    if (!minReadsReached) {
        std::cerr << " INFO: Warning: scanned reads (" << state.totalReads << ") below min_reads ("
                  << context.minReads << ")\n";
    }

    std::vector<std::uint64_t> panelCoveredIndexKmers(context.panelNames.size(), 0);
    std::vector<double> panelCoveredSpecificityMass(context.panelNames.size(), 0.0);
    double totalMatchedLookupKmers = static_cast<double>(state.matchedLookupKmers.size());

    for (const auto& encodedKmer : state.matchedLookupKmers) {
        auto hit = context.globalLookup.find(encodedKmer);
        if (hit == context.globalLookup.end()) {
            continue;
        }

        const std::vector<std::size_t>& panelIds = hit->second;
        double specificityWeight = 1.0 / static_cast<double>(panelIds.size());
        for (std::size_t j = 0; j < panelIds.size(); j = j + 1) {
            std::size_t panelId = panelIds[j];
            panelCoveredIndexKmers[panelId] = panelCoveredIndexKmers[panelId] + 1;
            panelCoveredSpecificityMass[panelId] =
                panelCoveredSpecificityMass[panelId] + specificityWeight;
        }
    }

    std::vector<PanelCandidateMetrics> panelCandidates;
    panelCandidates.reserve(context.panelNames.size());

    for (std::size_t i = 0; i < context.panelNames.size(); i = i + 1) {
        PanelCandidateMetrics metrics;
        metrics.panelId = i;
        metrics.matchedKmers = state.panelMatchedKmers[i];
        metrics.readsWithMatch = state.panelReadsWithMatch[i];
        metrics.coveredIndexKmers = panelCoveredIndexKmers[i];

        double panelCoverage = 0.0;
        if (context.panelIndexUniqueKmers[i] > 0) {
            panelCoverage = static_cast<double>(metrics.coveredIndexKmers) /
                            static_cast<double>(context.panelIndexUniqueKmers[i]);
        }

        double specificityPrecision = 0.0;
        if (totalMatchedLookupKmers > 0.0) {
            specificityPrecision = panelCoveredSpecificityMass[i] / totalMatchedLookupKmers;
        }

        // shared k-mers count less!
        if (context.panelIndexUniqueKmers[i] > 0) {
            metrics.panelKmerExtentPercent =
                (100.0 * static_cast<double>(metrics.coveredIndexKmers)) /
                static_cast<double>(context.panelIndexUniqueKmers[i]);
        }
        const double betaSquared = kFscoreBeta * kFscoreBeta;
        const double denominator = (betaSquared * specificityPrecision) + panelCoverage;
        if (denominator > 0.0) {
            metrics.candidateScore =
                100.0 * (1.0 + betaSquared) * specificityPrecision * panelCoverage / denominator;
        } else {
            metrics.candidateScore = 0.0;
        }

        if (state.totalReads > 0) {
            metrics.readSupportPercent =
                (100.0 * static_cast<double>(metrics.readsWithMatch)) /
                static_cast<double>(state.totalReads);
        }
        if (context.panelIndexUniqueKmers[i] > 0) {
            metrics.indexUniqueToPanelPercent =
                (100.0 * static_cast<double>(context.panelIndexUniqueToPanel[i])) /
                static_cast<double>(context.panelIndexUniqueKmers[i]);
            metrics.indexMeanPanelsPerKmer =
                static_cast<double>(context.panelIndexPanelMultiplicitySum[i]) /
                static_cast<double>(context.panelIndexUniqueKmers[i]);
        }

        panelCandidates.push_back(metrics);
    }

    std::sort(panelCandidates.begin(), panelCandidates.end(),
              PanelCandidateComparator{context.panelNames});

    bool ranSecondPass = false;
    if (panelCandidates.size() > 1 && panelCandidates[0].candidateScore > 0.0 &&
        panelCandidates[1].candidateScore > 0.0 &&
        panelCandidates[1].candidateScore >=
            (panelCandidates[0].candidateScore * kSecondPassScoreRatioTrigger)) {
        PairSpecificPassScores secondPassScores =
            runPairSpecificSecondPass(panelCandidates[0].panelId, panelCandidates[1].panelId,
                                     context.globalLookup, state.matchedLookupKmers);
        if (secondPassScores.firstTotal > 0 || secondPassScores.secondTotal > 0) {
            panelCandidates[0].candidateScore = secondPassScores.firstScore;
            panelCandidates[1].candidateScore = secondPassScores.secondScore;
            if (panelCandidates[1].candidateScore > panelCandidates[0].candidateScore) {
                std::swap(panelCandidates[0], panelCandidates[1]);
            }
            ranSecondPass = true;
            std::cout << " INFO: \t" << sampleName
                      << " second_pass=pair_specific"
                      << " panel_a=" << context.panelNames[panelCandidates[0].panelId]
                      << " panel_a_specific_score=" << panelCandidates[0].candidateScore
                      << " panel_b=" << context.panelNames[panelCandidates[1].panelId]
                      << " panel_b_specific_score=" << panelCandidates[1].candidateScore
                      << std::endl;
        }
    }

    rankRecord.sample = sampleName;
    rankRecord.totalReads = state.totalReads;
    rankRecord.minReadsReached = minReadsReached;
    rankRecord.maxReadsReached = maxReadsReached;

    if (!panelCandidates.empty() && panelCandidates[0].candidateScore > 0.0) {
        const PanelCandidateMetrics& best = panelCandidates[0];
        std::string nextPanel = "none";
        double nextScore = 0.0;
        if (panelCandidates.size() > 1) {
            nextPanel = context.panelNames[panelCandidates[1].panelId];
            nextScore = panelCandidates[1].candidateScore;
        }
        double scoreMargin = best.candidateScore - nextScore;

        rankRecord.bestPanel = context.panelNames[best.panelId];
        rankRecord.bestScore = best.candidateScore;
        rankRecord.scoreMarginVsNext = scoreMargin;
        rankRecord.bestMatchedKmers = best.matchedKmers;
        rankRecord.bestReadsWithMatch = best.readsWithMatch;
        rankRecord.bestReadSupportPercent = best.readSupportPercent;
        rankRecord.bestCoveredIndexKmers = best.coveredIndexKmers;
        rankRecord.bestPanelKmerExtentPercent = best.panelKmerExtentPercent;
        rankRecord.bestIndexUniqueToPanelPercent = best.indexUniqueToPanelPercent;
        rankRecord.bestIndexMeanPanelsPerKmer = best.indexMeanPanelsPerKmer;

        std::cout << " INFO: \t" << sampleName
                  << " best_panel=" << rankRecord.bestPanel
                  << " best_score=" << rankRecord.bestScore
                  << " second_panel=" << nextPanel
                  << " second_score=" << nextScore
                  << " scoring_pass=" << (ranSecondPass ? "pair_specific_second_pass" : "primary")
                  << std::endl;
    } 
    else {
        std::cout << " INFO: \t" << sampleName 
                  << " best_panel=none"
                  << " scanned_reads=" << rankRecord.totalReads << std::endl;
    }

    return true;
}

int runFindCommand(int argc, char** argv) {
    namespace fs = std::filesystem;

    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_dir", "Directory containing .bit/.2bit panel index files", "--index_dir", true);
    parser.add("fq1", "FASTQ read1 file (plain or .gz)", "--fq1", false);
    parser.add("fq2", "FASTQ read2 file (plain or .gz)", "--fq2", false);
    parser.add("fastq_list", "Text file with one FASTQ path (or fq1 fq2 pair) per line",
               "--fastq_list", false);
    parser.add("min_reads", "Minimum reads target to scan (default: 50000)", "--min_reads", false);
    parser.add("max_reads", "Maximum reads to scan before stopping (default: 1000000)", "--max_reads",
               false);
    parser.add("minimizer_window",
               "Minimizer window size in k-mers (default: 10, set 1 to disable minimizer subsampling)",
               "--minimizer_window", false);
    parser.add("min_kmer_entropy",
               "Skip k-mers with entropy below this value (range: 0.0 to 2.0, default: 0.0 disabled)",
               "--min_kmer_entropy", false);
    parser.add("output_path", "Output TSV path (default: panel_ranks.tsv)", "--output", false);
    parser.add("force_paired",
               "If set, do not skip paired/sample-mate FASTQ files when max_reads is reached",
               "--force_paired", false, true);

    if (!parser.parse()) {
        return 1;
    }

    bool hasFq1 = parser.parsed("fq1");
    bool hasFq2 = parser.parsed("fq2");
    bool hasFastqList = parser.parsed("fastq_list");
    if (hasFastqList && (hasFq1 || hasFq2)) {
        std::cerr << "find cannot combine --fastq_list with --fq1/--fq2\n";
        return 1;
    }
    if (!hasFastqList && !hasFq1 && !hasFq2) {
        std::cerr << "find requires --fastq_list or at least one input FASTQ (--fq1 and/or --fq2)\n";
        return 1;
    }

    std::vector<SampleFqEntry> sampleFqs;
    if (hasFastqList) {
        std::string fastqListPath = parser.get<std::string>("fastq_list");
        bool ok = loadFastqListFile(fastqListPath, sampleFqs);
        if (!ok) {
            return 1;
        }
    } else {
        SampleFqEntry singleSampleFq;
        singleSampleFq.label = "single_run";
        if (hasFq1) {
            singleSampleFq.fq1 = parser.get<std::string>("fq1");
        }
        if (hasFq2) {
            singleSampleFq.fq2 = parser.get<std::string>("fq2");
            singleSampleFq.hasFq2 = true;
        }
        if (singleSampleFq.fq1.empty() && singleSampleFq.hasFq2) {
            singleSampleFq.label = std::filesystem::path(singleSampleFq.fq2).stem().string();
        } 
        else if (!singleSampleFq.fq1.empty() && singleSampleFq.hasFq2) {
            singleSampleFq.label = std::filesystem::path(singleSampleFq.fq1).stem().string() + "+" +
                              std::filesystem::path(singleSampleFq.fq2).stem().string();
        } 
        else if (!singleSampleFq.fq1.empty()) {
            singleSampleFq.label = std::filesystem::path(singleSampleFq.fq1).stem().string();
        }
        if (singleSampleFq.label.empty()) {
            singleSampleFq.label = "single_run";
        }
        if (!singleSampleFq.fq1.empty()) {
            singleSampleFq.sampleKey = deriveSampleKeyFromFastqPath(singleSampleFq.fq1);
        } 
        else if (singleSampleFq.hasFq2) {
            singleSampleFq.sampleKey = deriveSampleKeyFromFastqPath(singleSampleFq.fq2);
        }
        if (singleSampleFq.sampleKey.empty()) {
            singleSampleFq.sampleKey = singleSampleFq.label;
        }
        sampleFqs.push_back(singleSampleFq);
    }

    std::uint64_t minReads = 50000;
    std::uint64_t maxReads = 1000000;
    std::size_t minimizerWindow = 10;
    double minKmerEntropy = 0.0;
    std::string panelRanksPath = "panel_ranks.tsv";

    if (parser.parsed("min_reads")) {
        long long value = parser.get<long long>("min_reads");
        if (value < 1) {
            std::cerr << "min_reads must be >= 1\n";
            return 1;
        }
        minReads = static_cast<std::uint64_t>(value);
    }

    if (parser.parsed("max_reads")) {
        long long value = parser.get<long long>("max_reads");
        if (value < 1) {
            std::cerr << "max_reads must be >= 1\n";
            return 1;
        }
        maxReads = static_cast<std::uint64_t>(value);
    }

    if (minReads > maxReads) {
        std::cerr << "min_reads cannot be greater than max_reads\n";
        return 1;
    }

    if (parser.parsed("minimizer_window")) {
        long long value = parser.get<long long>("minimizer_window");
        if (value < 1) {
            std::cerr << "minimizer_window must be >= 1\n";
            return 1;
        }
        minimizerWindow = static_cast<std::size_t>(value);
    }
    if (parser.parsed("min_kmer_entropy")) {
        double value = parser.get<double>("min_kmer_entropy");
        if (value < 0.0 || value > 2.0) {
            std::cerr << "min_kmer_entropy must be between 0.0 and 2.0\n";
            return 1;
        }
        minKmerEntropy = value;
    }
    if (parser.parsed("output_path")) {
        panelRanksPath = parser.get<std::string>("output_path");
    }
    bool useMinimizers = minimizerWindow > 1;
    bool useLowComplexityFilter = minKmerEntropy > 0.0;
    bool forcePaired = parser.get<bool>("force_paired");

    std::string indexDir = parser.get<std::string>("index_dir");

    std::error_code error;
    if (!fs::exists(indexDir, error) || !fs::is_directory(indexDir, error)) {
        std::cerr << "Invalid index_dir: " << indexDir << "\n";
        return 1;
    }

    std::cout << " INFO: Params --min_reads=" << minReads << " --max_reads=" << maxReads
              << " --minimizer_window=" << minimizerWindow
              << " --min_kmer_entropy=" << minKmerEntropy
              << " --force_paired=" << (forcePaired ? "true" : "false") << std::endl;

    std::size_t totalInputFastqFiles = countInputFastqFiles(sampleFqs);
    std::cout << " INFO: Found FASTQ files: " << totalInputFastqFiles << std::endl;
    std::cout << " INFO: Reading panel index from " << indexDir << std::endl;


    std::vector<std::string> panelNames;
    int kmerSize = -1;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> globalLookup;

    std::vector<std::string> indexFiles;
    for (const auto& entry : fs::directory_iterator(indexDir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        std::string ext = entry.path().extension().string();
        if (ext == ".bit" || ext == ".2bit") {
            indexFiles.push_back(entry.path().string());
        }
    }

    std::sort(indexFiles.begin(), indexFiles.end());

    if (indexFiles.empty()) {
        std::cerr << "No .bit/.2bit files found in " << indexDir << "\n";
        return 1;
    }

    // here allocating the new vector for unique index counts
    std::vector<std::uint64_t> indexUniqueCounts(indexFiles.size(), 0);
    std::uint64_t estimatedTotalUniqueKmers = 0;

    for (std::size_t i = 0; i < indexFiles.size(); i = i + 1) {

        std::ifstream in(indexFiles[i], std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Unable to open index file: " << indexFiles[i] << "\n";
            return 1;
        }

        std::int32_t kmerSize32 = 0;
        std::uint64_t uniqueCount = 0;
        in.read(reinterpret_cast<char*>(&kmerSize32), sizeof(kmerSize32));
        in.read(reinterpret_cast<char*>(&uniqueCount), sizeof(uniqueCount));

        if (!in.good()) {
            std::cerr << "Invalid index header: " << indexFiles[i] << "\n";
            return 1;
        }

        if (kmerSize32 < 1 || kmerSize32 > 32) {
            std::cerr << "Invalid kmer size in index file: " << indexFiles[i] << "\n";
            return 1;
        }
        if (kmerSize == -1) {
            kmerSize = kmerSize32;
        }
        else if (kmerSize32 != kmerSize) {
            std::cerr << "All index files must have same kmer size\n";
            return 1;
        }

        indexUniqueCounts[i] = uniqueCount;
        estimatedTotalUniqueKmers = estimatedTotalUniqueKmers + uniqueCount;
    }

    std::size_t reserveTarget = 0;
    const std::size_t maxSizeT = static_cast<std::size_t>(-1);
    if (estimatedTotalUniqueKmers > static_cast<std::uint64_t>(maxSizeT)) {
        reserveTarget = maxSizeT;
    } 
    else {
        reserveTarget = static_cast<std::size_t>(estimatedTotalUniqueKmers);
    }
    globalLookup.reserve(reserveTarget);
    panelNames.reserve(indexFiles.size());

    for (std::size_t i = 0; i < indexFiles.size(); i = i + 1) {
        std::ifstream in(indexFiles[i], std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Unable to open index file: " << indexFiles[i] << "\n";
            return 1;
        }

        std::cout << " INFO: Loading "
                  << fs::path(indexFiles[i]).filename().string() << std::endl;


        std::int32_t kmerSize32 = 0;
        std::uint64_t uniqueCount = 0;
        in.read(reinterpret_cast<char*>(&kmerSize32), sizeof(kmerSize32));
        in.read(reinterpret_cast<char*>(&uniqueCount), sizeof(uniqueCount));

        if (!in.good()) {
            std::cerr << "Invalid index header: " << indexFiles[i] << "\n";
            return 1;
        }
        if (kmerSize32 != kmerSize || uniqueCount != indexUniqueCounts[i]) {
            std::cerr << "Index metadata changed during load: " << indexFiles[i] << "\n";
            return 1;
        }

        std::size_t panelId = panelNames.size();
        panelNames.push_back(fs::path(indexFiles[i]).filename().string());

        for (std::uint64_t j = 0; j < uniqueCount; j = j + 1) {
            std::uint64_t encodedKmer = 0;
            in.read(reinterpret_cast<char*>(&encodedKmer), sizeof(encodedKmer));
            if (!in.good()) {
                std::cerr << "Unexpected EOF in index file: " << indexFiles[i] << "\n";
                return 1;
            }
            globalLookup[encodedKmer].push_back(panelId);
        }

    }

    std::vector<std::uint64_t> panelIndexUniqueKmers = indexUniqueCounts;
    std::vector<std::uint64_t> panelIndexUniqueToPanel(panelNames.size(), 0);
    std::vector<std::uint64_t> panelIndexPanelMultiplicitySum(panelNames.size(), 0);
    for (const auto& entry : globalLookup) {
        std::size_t panelsForKmer = entry.second.size();

        for (std::size_t j = 0; j < entry.second.size(); j = j + 1) {
            std::size_t panelId = entry.second[j];
            panelIndexPanelMultiplicitySum[panelId] =
                panelIndexPanelMultiplicitySum[panelId] +
                static_cast<std::uint64_t>(panelsForKmer);
            if (panelsForKmer == 1) {
                panelIndexUniqueToPanel[panelId] = panelIndexUniqueToPanel[panelId] + 1;
            }
        }
    }

    const std::size_t k = static_cast<std::size_t>(kmerSize);
    const std::array<std::uint8_t, 256>& lookup = utils::baseTo2bitLookup();
    std::uint64_t rollingMask = 0;
    if (k == 32) {
        rollingMask = ~static_cast<std::uint64_t>(0);
    } 
    else {
        rollingMask = (static_cast<std::uint64_t>(1) << (2 * k)) - 1;
    }
    std::array<double, 33> entropyContributionByCount = {};
    if (useLowComplexityFilter) {
        double kAsDouble = static_cast<double>(k);
        for (std::size_t count = 1; count <= k; count = count + 1) {
            double probability = static_cast<double>(count) / kAsDouble;
            entropyContributionByCount[count] = -probability * std::log2(probability);
        }
    }

    std::vector<PanelRankRecord> panelRankRecords;
    panelRankRecords.reserve(sampleFqs.size());
    std::unordered_set<std::string> samplesAtMaxReads;

    FindSampleContext sampleContext = {panelNames,
                                 globalLookup,
                                 panelIndexUniqueKmers,
                                 panelIndexUniqueToPanel,
                                 panelIndexPanelMultiplicitySum,
                                 lookup,
                                 indexDir,
                                 k,
                                 rollingMask,
                                 useMinimizers,
                                 minimizerWindow,
                                 minReads,
                                 maxReads,
                                 useLowComplexityFilter,
                                 minKmerEntropy,
                                 entropyContributionByCount,
                                 forcePaired};

    std::unordered_map<std::string, std::size_t> sampleOrdinalByKey;
    sampleOrdinalByKey.reserve(sampleFqs.size());
    std::size_t totalSamples = 0;
    for (std::size_t i = 0; i < sampleFqs.size(); i = i + 1) {
        const std::string& key = sampleFqs[i].sampleKey;
        if (sampleOrdinalByKey.find(key) == sampleOrdinalByKey.end()) {
            totalSamples = totalSamples + 1;
            sampleOrdinalByKey[key] = totalSamples;
        }
    }

    for (std::size_t sampleIndex = 0; sampleIndex < sampleFqs.size(); sampleIndex = sampleIndex + 1) {
        const SampleFqEntry& sampleFq = sampleFqs[sampleIndex];
        if (!forcePaired && samplesAtMaxReads.find(sampleFq.sampleKey) != samplesAtMaxReads.end()) {
            if (sampleFq.hasFq2) {
                std::cout << " INFO:  Skipping FASTQ pair: " << sampleFq.fq1 << " " << sampleFq.fq2
                          << " (reached max_reads)" << std::endl;
            }
            else if (!sampleFq.fq1.empty()) {
                std::cout << " INFO:  Skipping FASTQ: " << sampleFq.fq1 << " (reached max_reads)" << std::endl;
            } 
            else {
                std::cout << " INFO:  Skipping FASTQ: " << sampleFq.fq2 << " (reached max_reads)" << std::endl;
            }
            continue;
        }

        PanelRankRecord rankRecord;
        std::size_t sampleOrdinal = sampleOrdinalByKey[sampleFq.sampleKey];
        bool ok = scanOneSampleFq(sampleFq, sampleOrdinal, totalSamples, sampleContext, rankRecord);
        if (!ok) {
            return 1;
        }
        panelRankRecords.push_back(rankRecord);
        if (!forcePaired && rankRecord.maxReadsReached) {
            samplesAtMaxReads.insert(sampleFq.sampleKey);
        }
    }

    std::ofstream panelRanksOut(panelRanksPath);
    if (!panelRanksOut.is_open()) {
        std::cerr << "Unable to write panel ranks file: " << panelRanksPath << "\n";
        return 1;
    }
    panelRanksOut
        << "sample\tscanned_reads\tbest_panel\tbest_score\tscore_margin_vs_next\tbest_panel_covered_kmers\tbest_panel_covered_kmers_pct\n";
    for (std::size_t i = 0; i < panelRankRecords.size(); i = i + 1) {
        const PanelRankRecord& record = panelRankRecords[i];
        panelRanksOut << record.sample << "\t" << record.totalReads << "\t" << record.bestPanel
                      << "\t" << record.bestScore << "\t" << record.scoreMarginVsNext << "\t"
                      << record.bestCoveredIndexKmers << "\t"
                      << record.bestPanelKmerExtentPercent << "\n";
    }
    if (!panelRanksOut.good()) {
        std::cerr << "Error writing panel ranks file: " << panelRanksPath << "\n";
        return 1;
    }
    std::cout << " INFO: Finished! " << std::endl;

    return 0;
}
