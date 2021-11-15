#include <gtest/gtest.h>

#include <bamit/sample_functions.hpp>

TEST(sample_functions_test, sample_read_depth_test)
{
    std::filesystem::path input{DATADIR"simulated_chr1_small_golden.bam"};
    seqan3::sam_file_input input_file{input};

    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(input_file);
    bamit::EstimationResult<uint64_t> result = bamit::sample_read_depth(input_file, input_file.header(), node_list, 10);
    result.print();
}
