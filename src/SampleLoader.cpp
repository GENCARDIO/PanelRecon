// Intent: normalize FASTQ inputs and load sample manifests for `find`.
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "PanelFindInternal.hpp"
#include "utils.hpp"

namespace {

bool hasIlluminaLaneSuffix(const std::string& name, std::size_t lanePos) {
    if (name[lanePos] != '_' || name[lanePos + 1] != 'L') {
        return false;
    }

    for (std::size_t i = lanePos + 2; i < lanePos + 5; i = i + 1) {
        if (std::isdigit(static_cast<unsigned char>(name[i])) == 0) {
            return false;
        }
    }

    return true;
}

void sortSampleFqEntries(std::vector<SampleFqEntry>& entries) {
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
}

}  // namespace

// Return a normalized sample key from a FASTQ filename.
std::string sampleKeyFromFastqPath(const std::string& fastqPath) {
    std::string name = std::filesystem::path(fastqPath).filename().string();
    if (utils::endsWith(name, ".gz")) {
        name = name.substr(0, name.size() - 3);
    }
    if (utils::endsWith(name, ".fastq")) {
        name = name.substr(0, name.size() - 6);
    }
    else if (utils::endsWith(name, ".fq")) {
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
    else if (utils::endsWith(name, "_R1") || utils::endsWith(name, "_R2")) {
        name = name.substr(0, name.size() - 3);
        strippedIlluminaReadToken = true;
    }

    bool strippedIlluminaLaneToken = false;
    if (name.size() >= 5) {
        std::size_t lanePos = name.size() - 5;
        if (hasIlluminaLaneSuffix(name, lanePos)) {
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

// Read the FASTQ list file (exactly one FASTQ path per non-comment line).
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
        if (tokens.size() != 1) {
            std::cerr << "Invalid FASTQ list line " << lineNumber
                      << ": expected exactly one path, got " << tokens.size() << "\n";
            return false;
        }

        SampleFqEntry entry;
        entry.fq1 = utils::resolveListPathEntry(listPath, tokens[0]);
        entry.label = fs::path(entry.fq1).stem().string();
        entry.sampleKey = sampleKeyFromFastqPath(entry.fq1);
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

    sortSampleFqEntries(entries);
    return true;
}

// Read a TSV manifest with columns: sample, fq1, optional fq2.
bool loadFastqTsvFile(const std::string& tsvPath, std::vector<SampleFqEntry>& entries) {
    entries.clear();

    std::ifstream in(tsvPath);
    if (!in.is_open()) {
        std::cerr << "Unable to open FASTQ TSV file: " << tsvPath << "\n";
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

        std::vector<std::string> fields;
        std::istringstream lineStream(line);
        std::string field;
        while (std::getline(lineStream, field, '\t')) {
            fields.push_back(utils::trimWhitespace(field));
        }

        if (fields.size() < 2 || fields.size() > 3) {
            std::cerr << "Invalid FASTQ TSV line " << lineNumber
                      << ": expected 2 or 3 tab-separated columns, got " << fields.size()
                      << "\n";
            return false;
        }

        SampleFqEntry entry;
        entry.label = fields[0];
        entry.sampleKey = fields[0];

        if (entry.sampleKey.empty()) {
            std::cerr << "Invalid FASTQ TSV line " << lineNumber
                      << ": sample name cannot be empty\n";
            return false;
        }

        entry.fq1 = utils::resolveListPathEntry(tsvPath, fields[1]);
        if (entry.fq1.empty()) {
            std::cerr << "Invalid FASTQ TSV line " << lineNumber
                      << ": fq1 cannot be empty\n";
            return false;
        }

        if (fields.size() == 3 && !fields[2].empty()) {
            entry.fq2 = utils::resolveListPathEntry(tsvPath, fields[2]);
            entry.hasFq2 = true;
        }

        entries.push_back(entry);
    }

    if (entries.empty()) {
        std::cerr << "No FASTQ entries found in TSV file: " << tsvPath << "\n";
        return false;
    }

    sortSampleFqEntries(entries);
    return true;
}

// Count how many FASTQ files are present across all input entries.
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
