#ifndef UTILS_HPP
#define UTILS_HPP

#include <array>
#include <cstdint>
#include <string>

namespace utils {

    const std::array<std::uint8_t, 256>& baseTo2bitLookup();

    std::string trimWhitespace(const std::string& value);

    std::string resolveListPathEntry(const std::string& listFilePath,
                                    const std::string& entryPath);

}

#endif
