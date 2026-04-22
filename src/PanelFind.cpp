// Intent: orchestrate the `find` command from CLI parsing to final output.
#include <array>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "PanelFindInternal.hpp"
#include "parser.hpp"
#include "utils.hpp"

// Entry point for the `find` command.
int runFindCommand(int argc, char** argv) {
    namespace fs = std::filesystem;

    // Define and parse command-line options for the find workflow.
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_dir", "Directory containing .bit/.2bit panel index files", "--index_dir", true);
    parser.add("fq1", "FASTQ read1 file (plain or .gz)", "--fq1", false);
    parser.add("fq2", "FASTQ read2 file (plain or .gz)", "--fq2", false);
    parser.add("fastq_list", "Text file with exactly one FASTQ path per line",
               "--fastq_list", false);
    parser.add("fastq_tsv", "TSV file with columns: sample, fq1, optional fq2",
               "--fastq_tsv", false);
    parser.add("min_reads", "Minimum reads target to scan (default: 50000)", "--min_reads", false);
    parser.add("max_reads", "Maximum reads to scan before stopping (default: 1000000)", "--max_reads",
               false);
    parser.add("minimizer_window",
               "Minimizer window size in k-mers (default: 10, set 1 to disable minimizer subsampling)",
               "--minimizer_window", false);
    parser.add("min_kmer_entropy",
               "Skip k-mers with entropy below this value (range: 0.0 to 2.0, default: 0.0 disabled)",
               "--min_kmer_entropy", false);
    parser.add("min_score",
               "Minimum final score required to call a panel (default: 0.001, set 0 to disable)",
               "--min_score", false);
    parser.add("output_path", "Output TSV path (default: panel_ranks.tsv)", "--output", false);
    parser.add("force_paired",
               "If set, do not skip paired/sample-mate FASTQ files when max_reads is reached",
               "--force_paired", false, true);

    if (!parser.parse()) {
        return 1;
    }

    // Validate which input mode is being used: explicit FASTQs or a FASTQ list.
    bool hasFq1 = parser.parsed("fq1");
    bool hasFq2 = parser.parsed("fq2");
    bool hasFastqList = parser.parsed("fastq_list");
    bool hasFastqTsv = parser.parsed("fastq_tsv");
    int inputModeCount = 0;
    if (hasFastqList) {
        inputModeCount = inputModeCount + 1;
    }
    if (hasFastqTsv) {
        inputModeCount = inputModeCount + 1;
    }
    if (hasFq1 || hasFq2) {
        inputModeCount = inputModeCount + 1;
    }
    if (inputModeCount > 1) {
        std::cerr << "find cannot combine --fastq_list, --fastq_tsv, and --fq1/--fq2\n";
        return 1;
    }
    if (inputModeCount == 0) {
        std::cerr << "find requires --fastq_list, --fastq_tsv, or at least one input FASTQ (--fq1 and/or --fq2)\n";
        return 1;
    }

    // Build the sample list that will be scanned.
    std::vector<SampleFqEntry> sampleFqs;
    if (hasFastqList) {
        std::string fastqListPath = parser.get<std::string>("fastq_list");
        bool ok = loadFastqListFile(fastqListPath, sampleFqs);
        if (!ok) {
            return 1;
        }
    }
    else if (hasFastqTsv) {
        std::string fastqTsvPath = parser.get<std::string>("fastq_tsv");
        bool ok = loadFastqTsvFile(fastqTsvPath, sampleFqs);
        if (!ok) {
            return 1;
        }
    }
    else {
        // Build one synthetic sample entry from the directly provided FASTQ arguments.
        SampleFqEntry singleSampleFq;
        singleSampleFq.label = "single_run";
        if (hasFq1) {
            singleSampleFq.fq1 = parser.get<std::string>("fq1");
        }
        if (hasFq2) {
            singleSampleFq.fq2 = parser.get<std::string>("fq2");
            singleSampleFq.hasFq2 = true;
        }
        if (singleSampleFq.fq1.empty() && singleSampleFq.hasFq2) {
            singleSampleFq.label = fs::path(singleSampleFq.fq2).stem().string();
        }
        else if (!singleSampleFq.fq1.empty() && singleSampleFq.hasFq2) {
            singleSampleFq.label = fs::path(singleSampleFq.fq1).stem().string() + "+" +
                                   fs::path(singleSampleFq.fq2).stem().string();
        }
        else if (!singleSampleFq.fq1.empty()) {
            singleSampleFq.label = fs::path(singleSampleFq.fq1).stem().string();
        }
        if (singleSampleFq.label.empty()) {
            singleSampleFq.label = "single_run";
        }
        if (!singleSampleFq.fq1.empty()) {
            singleSampleFq.sampleKey = sampleKeyFromFastqPath(singleSampleFq.fq1);
        }
        else if (singleSampleFq.hasFq2) {
            singleSampleFq.sampleKey = sampleKeyFromFastqPath(singleSampleFq.fq2);
        }
        if (singleSampleFq.sampleKey.empty()) {
            singleSampleFq.sampleKey = singleSampleFq.label;
        }
        sampleFqs.push_back(singleSampleFq);
    }

    // Read optional settings and check their values.
    std::uint64_t minReads = 50000;
    std::uint64_t maxReads = 1000000;
    std::size_t minimizerWindow = 10;
    double minKmerEntropy = 0.0;
    double minPanelScore = 0.001;
    std::string panelRanksPath = "panel_ranks.tsv";

    if (parser.parsed("min_reads")) {
        long long value = parser.get<long long>("min_reads");
        if (value < 1) {
            std::cerr << "min_reads must be >= 1\n";
            return 1;
        }
        minReads = static_cast<std::uint64_t>(value);
    }

    if (parser.parsed("max_reads")) {
        long long value = parser.get<long long>("max_reads");
        if (value < 1) {
            std::cerr << "max_reads must be >= 1\n";
            return 1;
        }
        maxReads = static_cast<std::uint64_t>(value);
    }

    if (minReads > maxReads) {
        std::cerr << "min_reads cannot be greater than max_reads\n";
        return 1;
    }

    if (parser.parsed("minimizer_window")) {
        long long value = parser.get<long long>("minimizer_window");
        if (value < 1) {
            std::cerr << "minimizer_window must be >= 1\n";
            return 1;
        }
        minimizerWindow = static_cast<std::size_t>(value);
    }
    if (parser.parsed("min_kmer_entropy")) {
        double value = parser.get<double>("min_kmer_entropy");
        if (value < 0.0 || value > 2.0) {
            std::cerr << "min_kmer_entropy must be between 0.0 and 2.0\n";
            return 1;
        }
        minKmerEntropy = value;
    }
    if (parser.parsed("min_score")) {
        double value = parser.get<double>("min_score");
        if (value < 0.0) {
            std::cerr << "min_score must be >= 0.0\n";
            return 1;
        }
        minPanelScore = value;
    }
    if (parser.parsed("output_path")) {
        panelRanksPath = parser.get<std::string>("output_path");
    }

    // Precompute simple feature flags used later in the scanning path.
    bool useMinimizers = minimizerWindow > 1;
    bool useLowComplexityFilter = minKmerEntropy > 0.0;
    bool forcePaired = parser.get<bool>("force_paired");
    std::string indexDir = parser.get<std::string>("index_dir");

    // Validate the index directory and report the effective parameters.
    std::error_code error;
    if (!fs::exists(indexDir, error) || !fs::is_directory(indexDir, error)) {
        std::cerr << "Invalid index_dir: " << indexDir << "\n";
        return 1;
    }

    std::cout << " INFO: Params --min_reads=" << minReads << " --max_reads=" << maxReads
              << " --minimizer_window=" << minimizerWindow
              << " --min_kmer_entropy=" << minKmerEntropy
              << " --min_score=" << minPanelScore
              << " --force_paired=" << (forcePaired ? "true" : "false") << std::endl;

    std::size_t totalInputFastqFiles = countInputFastqFiles(sampleFqs);
    std::cout << " INFO: Found FASTQ files: " << totalInputFastqFiles << std::endl;
    std::cout << " INFO: Reading panel index from " << indexDir << std::endl;

    LoadedPanelIndexData panelIndexData;
    if (!loadPanelIndexData(indexDir, panelIndexData)) {
        return 1;
    }

    const std::size_t k = static_cast<std::size_t>(panelIndexData.kmerSize);
    // Map each base character to its 2-bit DNA encoding.
    const std::array<std::uint8_t, 256>& lookup = utils::baseTo2bitLookup();
    // Bit mask used to keep only the current k-mer bits in the rolling encoding.
    std::uint64_t rollingMask = 0;
    if (k == 32) {
        rollingMask = ~static_cast<std::uint64_t>(0);
    }
    else {
        rollingMask = (static_cast<std::uint64_t>(1) << (2 * k)) - 1;
    }

    std::array<double, 33> entropyContributionByCount =
        buildEntropyContributionByCount(k, useLowComplexityFilter);

    // Build the shared panel index data, scan settings, and read limits once.
    PanelIndexStruct panelIndex = {panelIndexData.panelNames,
                                   panelIndexData.globalLookup,
                                   panelIndexData.panelIndexUniqueKmers,
                                   panelIndexData.panelIndexUniqueToPanel,
                                   panelIndexData.panelIndexPanelsPerKmerSum};
    KmerScanSettings scanSettings = {lookup,
                                     k,
                                     rollingMask,
                                     useMinimizers,
                                     minimizerWindow,
                                     useLowComplexityFilter,
                                     minKmerEntropy,
                                     entropyContributionByCount};
    ScanLimits scanLimits = {minReads, maxReads, minPanelScore, forcePaired};

    std::vector<PanelRankRecord> panelRankRecords;
    if (!scanAllSamples(sampleFqs, panelIndex, scanSettings, scanLimits, panelRankRecords)) {
        return 1;
    }

    if (!writePanelRanksFile(panelRanksPath, panelRankRecords)) {
        return 1;
    }

    std::cout << " INFO: Finished! " << std::endl;
    return 0;
}
