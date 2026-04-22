// Intent: hold shared small helpers used across indexing and finding.
#include "utils.hpp"

#include <cctype>
#include <filesystem>

namespace utils {

namespace {

std::array<std::uint8_t, 256> createBaseTo2bitLookup() {
    std::array<std::uint8_t, 256> lookup{};
    lookup.fill(255);
    lookup['A'] = 0;
    lookup['a'] = 0;
    lookup['C'] = 1;
    lookup['c'] = 1;
    lookup['G'] = 2;
    lookup['g'] = 2;
    lookup['T'] = 3;
    lookup['t'] = 3;
    lookup['U'] = 3;
    lookup['u'] = 3;
    return lookup;
}

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

}  // namespace

const std::array<std::uint8_t, 256>& baseTo2bitLookup() {
    static const std::array<std::uint8_t, 256> lookup = createBaseTo2bitLookup();
    return lookup;
}

bool endsWith(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string reverseComplement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());

    for (std::size_t i = sequence.size(); i > 0; i = i - 1) {
        result.push_back(complementBase(sequence[i - 1]));
    }

    return result;
}

std::string trimWhitespace(const std::string& value) {
    std::size_t begin = 0;
    while (begin < value.size() &&
           std::isspace(static_cast<unsigned char>(value[begin])) != 0) {
        begin = begin + 1;
    }

    if (begin == value.size()) {
        return "";
    }

    std::size_t end = value.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(value[end - 1])) != 0) {
        end = end - 1;
    }

    return value.substr(begin, end - begin);
}

std::string resolveListPathEntry(const std::string& listFilePath, const std::string& entryPath) {
    namespace fs = std::filesystem;
    fs::path path(entryPath);
    if (path.is_absolute()) {
        return path.lexically_normal().string();
    }
    fs::path resolved = fs::path(listFilePath).parent_path() / path;
    return resolved.lexically_normal().string();
}

}  // namespace utils
