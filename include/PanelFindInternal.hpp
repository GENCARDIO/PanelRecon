#ifndef PANEL_FIND_INTERNAL_HPP
#define PANEL_FIND_INTERNAL_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

struct SampleFqEntry {
    std::string label;
    std::string sampleKey;
    std::string fq1;
    std::string fq2;
    bool hasFq2 = false;
};

struct PanelRankRecord {
    std::string sample;
    std::string bestPanel = "none";
    double bestScore = 0.0;
    double scoreMarginVsNext = 0.0;
    std::uint64_t totalReads = 0;
    bool minReadsReached = false;
    bool maxReadsReached = false;
    std::uint64_t bestMatchedKmers = 0;
    std::uint64_t bestCoveredPanelKmers = 0;
    double bestPanelKmerExtentPercent = 0.0;
    double bestIndexUniqueToPanelPercent = 0.0;
    double bestIndexMeanPanelsPerKmer = 0.0;
};

struct PanelIndexStruct {
    const std::vector<std::string>& panelNames;
    // Map each encoded k-mer to the panel IDs whose indexes contain it.
    const std::unordered_map<std::uint64_t, std::vector<std::size_t>>& globalLookup;
    const std::vector<std::uint64_t>& panelIndexUniqueKmers;
    const std::vector<std::uint64_t>& panelIndexUniqueToPanel;
    const std::vector<std::uint64_t>& panelIndexPanelsPerKmerSum;
};

struct KmerScanSettings {
    const std::array<std::uint8_t, 256>& lookup;
    std::size_t k = 0;
    std::uint64_t rollingMask = 0;
    bool useMinimizers = false;
    std::size_t minimizerWindow = 1;
    bool useLowComplexityFilter = false;
    double minKmerEntropy = 0.0;
    std::array<double, 33> entropyContributionByCount = {};
};

struct ScanLimits {
    std::uint64_t minReads = 0;
    std::uint64_t maxReads = 0;
    double minPanelScore = 0.0;
    bool forcePaired = false;
};

struct LoadedPanelIndexData {
    std::vector<std::string> panelNames;
    int kmerSize = -1;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> globalLookup;
    std::vector<std::uint64_t> panelIndexUniqueKmers;
    std::vector<std::uint64_t> panelIndexUniqueToPanel;
    std::vector<std::uint64_t> panelIndexPanelsPerKmerSum;
};

std::string sampleKeyFromFastqPath(const std::string& fastqPath);
bool loadFastqListFile(const std::string& listPath, std::vector<SampleFqEntry>& entries);
bool loadFastqTsvFile(const std::string& tsvPath, std::vector<SampleFqEntry>& entries);
std::size_t countInputFastqFiles(const std::vector<SampleFqEntry>& sampleFqs);

bool loadPanelIndexData(const std::string& indexDir, LoadedPanelIndexData& panelIndexData);
std::array<double, 33> buildEntropyContributionByCount(std::size_t k,
                                                       bool useLowComplexityFilter);
bool scanAllSamples(const std::vector<SampleFqEntry>& sampleFqs,
                    const PanelIndexStruct& panelIndex,
                    const KmerScanSettings& scanSettings,
                    const ScanLimits& scanLimits,
                    std::vector<PanelRankRecord>& panelRankRecords);
bool writePanelRanksFile(const std::string& panelRanksPath,
                         const std::vector<PanelRankRecord>& panelRankRecords);

#endif
