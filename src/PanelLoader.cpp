// Intent: load panel index files into the shared lookup structures for `find`.
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "PanelFindInternal.hpp"

// Load all panel index files and build the shared lookup tables used during scanning.
bool loadPanelIndexData(const std::string& indexDir, LoadedPanelIndexData& panelIndexData) {
    namespace fs = std::filesystem;

    panelIndexData = LoadedPanelIndexData{};

    std::vector<std::string> indexFiles;
    for (const auto& dirEntry : fs::directory_iterator(indexDir)) {
        if (!dirEntry.is_regular_file()) {
            continue;
        }
        const fs::path& indexFile = dirEntry.path();
        std::string ext = indexFile.extension().string();
        if (ext == ".bit" || ext == ".2bit") {
            indexFiles.push_back(indexFile.string());
        }
    }

    std::sort(indexFiles.begin(), indexFiles.end());
    if (indexFiles.empty()) {
        std::cerr << "No .bit/.2bit files found in " << indexDir << "\n";
        return false;
    }

    std::vector<std::uint64_t> indexUniqueCounts(indexFiles.size(), 0);
    std::uint64_t estimatedTotalUniqueKmers = 0;
    for (std::size_t i = 0; i < indexFiles.size(); i = i + 1) {
        std::ifstream in(indexFiles[i], std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Unable to open index file: " << indexFiles[i] << "\n";
            return false;
        }

        std::int32_t kmerSize32 = 0;
        std::uint64_t uniqueCount = 0;
        in.read(reinterpret_cast<char*>(&kmerSize32), sizeof(kmerSize32));
        in.read(reinterpret_cast<char*>(&uniqueCount), sizeof(uniqueCount));
        if (!in.good()) {
            std::cerr << "Invalid index header: " << indexFiles[i] << "\n";
            return false;
        }

        if (kmerSize32 < 1 || kmerSize32 > 32) {
            std::cerr << "Invalid kmer size in index file: " << indexFiles[i] << "\n";
            return false;
        }
        if (panelIndexData.kmerSize == -1) {
            panelIndexData.kmerSize = kmerSize32;
        }
        else if (kmerSize32 != panelIndexData.kmerSize) {
            std::cerr << "All index files must have same kmer size\n";
            return false;
        }

        indexUniqueCounts[i] = uniqueCount;
        estimatedTotalUniqueKmers = estimatedTotalUniqueKmers + uniqueCount;
    }

    // Clamp the reserve size to what size_t can represent on this build.
    std::size_t reserveTarget = 0;
    const std::size_t maxSizeT = static_cast<std::size_t>(-1);
    if (estimatedTotalUniqueKmers > static_cast<std::uint64_t>(maxSizeT)) {
        reserveTarget = maxSizeT;
    }
    else {
        reserveTarget = static_cast<std::size_t>(estimatedTotalUniqueKmers);
    }

    // Reserve preallocates storage so these containers grow less during index loading.
    panelIndexData.globalLookup.reserve(reserveTarget);
    panelIndexData.panelNames.reserve(indexFiles.size());

    for (std::size_t i = 0; i < indexFiles.size(); i = i + 1) {
        std::ifstream in(indexFiles[i], std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Unable to open index file: " << indexFiles[i] << "\n";
            return false;
        }

        std::cout << " INFO: Loading " << fs::path(indexFiles[i]).filename().string() << std::endl;

        // Read the fixed header written by the index command.
        // kmerSize32 stores the index k-mer size.
        // uniqueCount stores how many unique k-mers are in this file.
        std::int32_t kmerSize32 = 0;
        std::uint64_t uniqueCount = 0;
        in.read(reinterpret_cast<char*>(&kmerSize32), sizeof(kmerSize32));
        in.read(reinterpret_cast<char*>(&uniqueCount), sizeof(uniqueCount));
        if (!in.good()) {
            std::cerr << "Invalid index header: " << indexFiles[i] << "\n";
            return false;
        }
        if (kmerSize32 != panelIndexData.kmerSize || uniqueCount != indexUniqueCounts[i]) {
            std::cerr << "Index metadata changed during load: " << indexFiles[i] << "\n";
            return false;
        }

        std::size_t panelId = panelIndexData.panelNames.size();
        panelIndexData.panelNames.push_back(fs::path(indexFiles[i]).filename().string());

        for (std::uint64_t j = 0; j < uniqueCount; j = j + 1) {
            // Read one encoded k-mer from the index body.
            std::uint64_t encodedKmer = 0;
            in.read(reinterpret_cast<char*>(&encodedKmer), sizeof(encodedKmer));
            if (!in.good()) {
                std::cerr << "Unexpected EOF in index file: " << indexFiles[i] << "\n";
                return false;
            }
            panelIndexData.globalLookup[encodedKmer].push_back(panelId);
        }
    }

    // Total unique k-mers stored for each panel index.
    panelIndexData.panelIndexUniqueKmers = indexUniqueCounts;
    // Count of k-mers that appear in only one panel.
    panelIndexData.panelIndexUniqueToPanel.assign(panelIndexData.panelNames.size(), 0);
    // Sum of how many panels each indexed k-mer belongs to.
    panelIndexData.panelIndexPanelsPerKmerSum.assign(panelIndexData.panelNames.size(), 0);

    for (const auto& kmerPanelsEntry : panelIndexData.globalLookup) {
        std::size_t panelsForKmer = kmerPanelsEntry.second.size();
        for (std::size_t j = 0; j < kmerPanelsEntry.second.size(); j = j + 1) {
            std::size_t panelId = kmerPanelsEntry.second[j];
            panelIndexData.panelIndexPanelsPerKmerSum[panelId] =
                panelIndexData.panelIndexPanelsPerKmerSum[panelId] +
                static_cast<std::uint64_t>(panelsForKmer);
            if (panelsForKmer == 1) {
                panelIndexData.panelIndexUniqueToPanel[panelId] =
                    panelIndexData.panelIndexUniqueToPanel[panelId] + 1;
            }
        }
    }

    return true;
}
