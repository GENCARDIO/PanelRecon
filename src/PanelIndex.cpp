#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "RefFasta.h"
#include "parser.hpp"
#include "utils.hpp"


// namespace {

bool baseTo2bit(char base, std::uint64_t& code) {
    const std::array<std::uint8_t, 256>& lookup = utils::baseTo2bitLookup();
    unsigned char baseIndex = base;
    std::uint8_t value = lookup[baseIndex];
    if (value > 3) {
        return false;
    }
    code = value;
    return true;
}

bool encodeKmerAt(const std::string& sequence, std::size_t startPos, std::size_t kmerSize,
                  std::uint64_t& encodedKmer) {
    encodedKmer = 0;

    for (std::size_t i = 0; i < kmerSize; i = i + 1) {
        std::uint64_t code = 0;
        bool ok = baseTo2bit(sequence[startPos + i], code);
        if (!ok) {
            return false;
        }

        encodedKmer = (encodedKmer << 2) | code;
    }

    return true;
}

struct KmerBuildResult {
    std::unordered_set<std::uint64_t> uniqueKmers;
    std::size_t regionsWithKmers = 0;
    std::size_t totalKmers = 0;
    std::size_t singletonKmers = 0;
    std::size_t sharedKmers = 0;
    std::size_t maxKmerMultiplicity = 0;
    std::size_t kmersInOneRegion = 0;
    std::size_t kmersInMultipleRegions = 0;
    std::size_t maxRegionsPerKmer = 0;
};

void addSelectedKmer(std::uint64_t encodedKmer,
                     std::unordered_map<std::uint64_t, std::uint32_t>& kmerOccurrenceCount,
                     std::unordered_set<std::uint64_t>& regionSeenKmers,
                     std::size_t& totalKmers) {
    kmerOccurrenceCount[encodedKmer] = kmerOccurrenceCount[encodedKmer] + 1;
    totalKmers = totalKmers + 1;
    regionSeenKmers.insert(encodedKmer);
}

void addSelectedKmersFromSegment(
    const std::vector<std::uint64_t>& segmentKmers, std::size_t minimizerWindow,
    std::unordered_map<std::uint64_t, std::uint32_t>& kmerOccurrenceCount,
    std::unordered_set<std::uint64_t>& regionSeenKmers, std::size_t& totalKmers) {
    if (segmentKmers.empty()) {
        return;
    }

    if (minimizerWindow <= 1 || segmentKmers.size() == 1) {
        for (std::size_t j = 0; j < segmentKmers.size(); j = j + 1) {
            addSelectedKmer(segmentKmers[j], kmerOccurrenceCount, regionSeenKmers, totalKmers);
        }
        return;
    }

    if (segmentKmers.size() < minimizerWindow) {
        std::uint64_t minValue = segmentKmers[0];
        for (std::size_t j = 1; j < segmentKmers.size(); j = j + 1) {
            if (segmentKmers[j] < minValue) {
                minValue = segmentKmers[j];
            }
        }
        addSelectedKmer(minValue, kmerOccurrenceCount, regionSeenKmers, totalKmers);
        return;
    }

    const std::size_t noPickedIndex = static_cast<std::size_t>(-1);
    std::size_t lastPickedIndex = noPickedIndex;
    for (std::size_t start = 0; start + minimizerWindow <= segmentKmers.size();
         start = start + 1) {
        std::size_t minIndex = start;
        std::uint64_t minValue = segmentKmers[start];
        std::size_t end = start + minimizerWindow;
        for (std::size_t j = start + 1; j < end; j = j + 1) {
            if (segmentKmers[j] < minValue) {
                minValue = segmentKmers[j];
                minIndex = j;
            }
        }

        if (minIndex == lastPickedIndex) {
            continue;
        }

        addSelectedKmer(minValue, kmerOccurrenceCount, regionSeenKmers, totalKmers);
        lastPickedIndex = minIndex;
    }
}

KmerBuildResult buildUniqueKmerLookup(const std::string& bedFile, RefFasta& reference,
                                      int kmerSize, std::size_t minimizerWindow) {
    KmerBuildResult result;
    std::unordered_map<std::uint64_t, std::uint32_t> kmerOccurrenceCount;
    std::unordered_map<std::uint64_t, std::uint32_t> kmerRegionCount;

    std::ifstream bedStream(bedFile);
    if (!bedStream.is_open()) {
        std::cerr << "Unable to open BED file: " << bedFile << "\n";
        return result;
    }

    std::size_t k = static_cast<std::size_t>(kmerSize);
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(bedStream, line)) {
        lineNumber = lineNumber + 1;

        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            continue;
        }

        std::istringstream lineStream(line);
        std::string chrom;
        std::int64_t start = 0;
        std::int64_t end = 0;

        if (!(lineStream >> chrom >> start >> end)) {
            std::cerr << "Skipping invalid BED line " << lineNumber << ": " << line << "\n";
            continue;
        }

        if (start < 0 || end <= start) {
            std::cerr << "Skipping invalid interval at line " << lineNumber << ": " << line << "\n";
            continue;
        }

        if (start >= static_cast<std::int64_t>(INT_MAX) ||
            end >= static_cast<std::int64_t>(INT_MAX)) {
            std::cerr << "Skipping out-of-range interval at line " << lineNumber << ": " << line
                      << "\n";
            continue;
        }

        int oneBasedStart = static_cast<int>(start) + 1;
        int oneBasedEnd = static_cast<int>(end);
        std::string sequence = reference.fetchSequence(chrom, oneBasedStart, oneBasedEnd);

        if (sequence.empty()) {
            std::cerr << "No sequence fetched at line " << lineNumber << ": " << line << "\n";
            continue;
        }

        if (sequence.size() < k) {
            continue;
        }

        result.regionsWithKmers = result.regionsWithKmers + 1;
        std::unordered_set<std::uint64_t> regionSeenKmers;
        regionSeenKmers.reserve(sequence.size() - k + 1);
        std::vector<std::uint64_t> segmentKmers;
        segmentKmers.reserve(sequence.size() - k + 1);

        for (std::size_t startPos = 0; startPos + k <= sequence.size(); startPos = startPos + 1) {
            std::uint64_t encodedKmer = 0;
            bool ok = encodeKmerAt(sequence, startPos, k, encodedKmer);
            if (!ok) {
                addSelectedKmersFromSegment(segmentKmers, minimizerWindow, kmerOccurrenceCount,
                                            regionSeenKmers, result.totalKmers);
                segmentKmers.clear();
                continue;
            }

            segmentKmers.push_back(encodedKmer);
        }
        addSelectedKmersFromSegment(segmentKmers, minimizerWindow, kmerOccurrenceCount,
                                    regionSeenKmers, result.totalKmers);

        for (const auto& encodedKmer : regionSeenKmers) {
            kmerRegionCount[encodedKmer] = kmerRegionCount[encodedKmer] + 1;
        }
    }

    result.uniqueKmers.reserve(kmerOccurrenceCount.size());
    for (const auto& entry : kmerOccurrenceCount) {
        result.uniqueKmers.insert(entry.first);
        std::size_t occurrences = static_cast<std::size_t>(entry.second);
        if (occurrences == 1) {
            result.singletonKmers = result.singletonKmers + 1;
        } else {
            result.sharedKmers = result.sharedKmers + 1;
        }
        if (occurrences > result.maxKmerMultiplicity) {
            result.maxKmerMultiplicity = occurrences;
        }
    }

    for (const auto& entry : kmerRegionCount) {
        std::size_t regions = static_cast<std::size_t>(entry.second);
        if (regions == 1) {
            result.kmersInOneRegion = result.kmersInOneRegion + 1;
        } else {
            result.kmersInMultipleRegions = result.kmersInMultipleRegions + 1;
        }
        if (regions > result.maxRegionsPerKmer) {
            result.maxRegionsPerKmer = regions;
        }
    }

    return result;
}

bool writeBitIndexFile(const std::unordered_set<std::uint64_t>& kmerLookup, int kmerSize,
                       const std::string& bedFile, const std::string& outputDir,
                       std::string& outputPath) {
    namespace fs = std::filesystem;
    std::error_code error;

    bool exists = fs::exists(outputDir, error);
    if (error) {
        std::cerr << "Unable to check output directory: " << outputDir << "\n";
        return false;
    }

    if (!exists) {
        bool created = fs::create_directories(outputDir, error);
        if (!created || error) {
            std::cerr << "Unable to create output directory: " << outputDir << "\n";
            return false;
        }
    }

    std::string panelName = fs::path(bedFile).stem().string();
    if (panelName.empty()) {
        panelName = "panel";
    }

    outputPath = (fs::path(outputDir) / (panelName + ".2bit")).string();

    std::ofstream out(outputPath, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Unable to open output file: " << outputPath << "\n";
        return false;
    }

    std::int32_t kmerSize32 = static_cast<std::int32_t>(kmerSize);
    std::uint64_t uniqueCount = static_cast<std::uint64_t>(kmerLookup.size());
    out.write(reinterpret_cast<const char*>(&kmerSize32), sizeof(kmerSize32));
    out.write(reinterpret_cast<const char*>(&uniqueCount), sizeof(uniqueCount));

    std::vector<std::uint64_t> sortedKmers;
    sortedKmers.reserve(kmerLookup.size());
    for (const auto& encodedKmer : kmerLookup) {
        sortedKmers.push_back(encodedKmer);
    }
    std::sort(sortedKmers.begin(), sortedKmers.end());

    for (std::size_t i = 0; i < sortedKmers.size(); i = i + 1) {
        std::uint64_t encodedKmer = sortedKmers[i];
        out.write(reinterpret_cast<const char*>(&encodedKmer), sizeof(encodedKmer));
    }

    if (!out.good()) {
        std::cerr << "Error writing output file: " << outputPath << "\n";
        return false;
    }

    return true;
}

bool loadBedListFile(const std::string& listPath, std::vector<std::string>& bedPaths) {
    bedPaths.clear();

    std::ifstream in(listPath);
    if (!in.is_open()) {
        std::cerr << "Unable to open BED list file: " << listPath << "\n";
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        std::string trimmed = utils::trimWhitespace(line);
        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }
        bedPaths.push_back(utils::resolveListPathEntry(listPath, trimmed));
    }

    if (bedPaths.empty()) {
        std::cerr << "No BED paths found in list file: " << listPath << "\n";
        return false;
    }

    return true;
}

std::string panelIndexPathForBed(const std::string& bedPath, const std::string& outputDir) {
    namespace fs = std::filesystem;
    std::string panelName = fs::path(bedPath).stem().string();
    if (panelName.empty()) {
        panelName = "panel";
    }
    return (fs::path(outputDir) / (panelName + ".2bit")).string();
}

// }

int runIndexCommand(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("bed", "Input BED file path", "--bed", false);
    parser.add("bed_list", "Text file with one BED path per line", "--bed_list", false);
    parser.add("fasta", "Input FASTA file path", "--fasta", true);
    parser.add("output_dir", "Directory for outputs", "--output_dir", true);
    parser.add("kmer_size", "K-mer size (positive integer, default: 31)", "--kmer_size", false);
    parser.add("minimizer_window",
               "Minimizer window size in k-mers for indexing (default: 1, disabled)",
               "--minimizer_window", false);

    if (!parser.parse()) {
        return 1;
    }

    bool hasBed = parser.parsed("bed");
    bool hasBedList = parser.parsed("bed_list");
    if (hasBed == hasBedList) {
        std::cerr << "index requires exactly one of --bed or --bed_list\n";
        return 1;
    }

    std::string fasta = parser.get<std::string>("fasta");
    std::string outputDir = parser.get<std::string>("output_dir");
    int kmerSize = 31;
    std::size_t minimizerWindow = 1;
    if (parser.parsed("kmer_size")) {
        kmerSize = parser.get<int>("kmer_size");
    }
    if (parser.parsed("minimizer_window")) {
        long long value = parser.get<long long>("minimizer_window");
        if (value < 1) {
            std::cerr << "minimizer_window must be >= 1\n";
            return 1;
        }
        minimizerWindow = static_cast<std::size_t>(value);
    }

    if (kmerSize < 1 || kmerSize > 32) {
        std::cerr << "kmer_size must be between 1 and 32\n";
        return 1;
    }

    std::vector<std::string> bedPaths;
    if (hasBed) {
        bedPaths.push_back(parser.get<std::string>("bed"));
    } else {
        std::string bedListPath = parser.get<std::string>("bed_list");
        bool ok = loadBedListFile(bedListPath, bedPaths);
        if (!ok) {
            return 1;
        }
    }

    std::cout << "[index] start fasta=" << fasta << " k=" << kmerSize
              << " minimizer_window=" << minimizerWindow
              << " panels=" << bedPaths.size() << " output_dir=" << outputDir << "\n";

    RefFasta reference(fasta);

    std::size_t totalPanelsIndexed = 0;
    std::size_t totalPanelsSkipped = 0;

    for (std::size_t bedId = 0; bedId < bedPaths.size(); bedId = bedId + 1) {
        const std::string& bedPath = bedPaths[bedId];
        std::string expectedIndexPath = panelIndexPathForBed(bedPath, outputDir);

        if (hasBedList) {
            std::error_code existsError;
            if (std::filesystem::exists(expectedIndexPath, existsError)) {
                std::cout << "[index] skip panel=" << bedPath << " (already indexed)\n";
                totalPanelsSkipped = totalPanelsSkipped + 1;
                continue;
            }
        }

        std::cout << "[index] panel " << (bedId + 1) << "/" << bedPaths.size()
                  << " " << bedPath << "\n";

        KmerBuildResult buildResult =
            buildUniqueKmerLookup(bedPath, reference, kmerSize, minimizerWindow);

        std::string bitFilePath;
        bool wroteFile =
            writeBitIndexFile(buildResult.uniqueKmers, kmerSize, bedPath, outputDir, bitFilePath);
        if (!wroteFile) {
            return 1;
        }

        std::cout << "[index] done panel=" << bedPath
                  << " unique_kmers=" << buildResult.uniqueKmers.size()
                  << " index_file=" << bitFilePath << "\n";

        totalPanelsIndexed = totalPanelsIndexed + 1;
    }

    std::cout << "[index] done indexed=" << totalPanelsIndexed
              << " skipped_existing=" << totalPanelsSkipped
              << " total=" << bedPaths.size() << "\n";
    return 0;
}
