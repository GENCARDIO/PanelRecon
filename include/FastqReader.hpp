#ifndef FASTQ_READER_HPP
#define FASTQ_READER_HPP

#include <fstream>
#include <string>

#include <zlib.h>

namespace fastq {

enum class FastqReadStatus {
    RecordRead,
    EndOfFile,
    TruncatedRecord,
    OpenError,
};

struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

class FastqReader {
public:
    explicit FastqReader(const std::string& path);
    ~FastqReader();

    bool isOpen() const;
    FastqReadStatus readNext(FastqRecord& record);

private:
    bool readLine(std::string& line);
    static void trimLineEnding(std::string& line);

    bool isGzip_ = false;
    std::ifstream plainFile_;
    gzFile gzFile_ = nullptr;
};

}  // namespace fastq

#endif
