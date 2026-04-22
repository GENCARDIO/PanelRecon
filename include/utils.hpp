#ifndef UTILS_HPP
#define UTILS_HPP

#include <array>
#include <cstdint>
#include <string>

namespace utils {

    const std::array<std::uint8_t, 256>& baseTo2bitLookup();

    bool endsWith(const std::string& value, const std::string& suffix);

    std::string reverseComplement(const std::string& sequence);

    std::string trimWhitespace(const std::string& value);

    std::string resolveListPathEntry(const std::string& listFilePath,
                                    const std::string& entryPath);

}

#endif
