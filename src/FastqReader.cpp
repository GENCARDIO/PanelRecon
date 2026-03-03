#include "FastqReader.hpp"

#include <filesystem>

namespace fastq {

FastqReader::FastqReader(const std::string& path) {
    isGzip_ = std::filesystem::path(path).extension() == ".gz";
    if (isGzip_) {
        gzFile_ = gzopen(path.c_str(), "rb");
    } else {
        plainFile_.open(path);
    }
}

FastqReader::~FastqReader() {
    if (gzFile_ != nullptr) {
        gzclose(gzFile_);
    }
}

bool FastqReader::isOpen() const {
    if (isGzip_) {
        return gzFile_ != nullptr;
    }
    return plainFile_.is_open();
}

FastqReadStatus FastqReader::readNext(FastqRecord& record) {
    if (!isOpen()) {
        return FastqReadStatus::OpenError;
    }

    if (!readLine(record.header)) {
        return FastqReadStatus::EndOfFile;
    }
    if (!readLine(record.sequence)) {
        return FastqReadStatus::TruncatedRecord;
    }
    if (!readLine(record.plus)) {
        return FastqReadStatus::TruncatedRecord;
    }
    if (!readLine(record.quality)) {
        return FastqReadStatus::TruncatedRecord;
    }

    return FastqReadStatus::RecordRead;
}

bool FastqReader::readLine(std::string& line) {
    // Reset the output string before reading a new line
    line.clear();

    if (!isGzip_) {
        // If reading fails (EOF or error), signal no more data
        if (!std::getline(plainFile_, line)) {
            return false;
        }

        // Remove trailing '\r' if present
        trimLineEnding(line);
        return true;
    }

    // Gzipped path: read chunk by chunk because gzgets reads into a C buffer
    char buffer[4096];
    while (true) {
        // Try to read up to sizeof(buffer)-1 chars from the gzip stream
        char* status = gzgets(gzFile_, buffer, sizeof(buffer));

        //this handles files that end without a final newline
        if (status == nullptr) {
            return !line.empty();
        }

        // Append the new chunk to the growing line
        line += buffer;

        // Stop once we have read the newline character
        if (!line.empty() && line.back() == '\n') {
            break;
        }
    }

    // Remove trailing newline / carriage return
    trimLineEnding(line);
    return true;
}

void FastqReader::trimLineEnding(std::string& line) {
    if (!line.empty() && line.back() == '\n') {
        line.pop_back();
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
}

} 
