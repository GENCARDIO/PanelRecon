#include <iostream>
#include <string>

#include "commands.hpp"

void printMainUsage(const std::string& programName) {
    std::cerr << "Usage:\n";
    std::cerr << "  " << programName
              << " index [--bed <panel.bed> | --bed_list <beds.txt>] --fasta <ref.fa> --output_dir <out_dir> [--kmer_size <k>]\n";
    std::cerr << "  " << programName
              << " find --index_dir <dir_with_2bit_or_bit_files> ([--fq1 <reads_1.fq(.gz)>] [--fq2 <reads_2.fq(.gz)>] | [--fastq_list <fastqs.txt>]) [--min_reads <n>] [--max_reads <n>] [--minimizer_window <w>] [--min_kmer_entropy <e>] [--output <path>] [--force_paired]\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printMainUsage(argv[0]);
        return 1;
    }

    std::string command = argv[1];
    if (command == "index") {
        return runIndexCommand(argc - 1, argv + 1);
    }

    if (command == "find") {
        return runFindCommand(argc - 1, argv + 1);
    }

    std::cerr << "Unknown command: " << command << "\n";
    printMainUsage(argv[0]);
    return 1;
}
