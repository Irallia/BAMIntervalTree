#pragma once

#include <algorithm>
#include <execution>
#include <random>
#include <numeric>

#include <bamit/Record.hpp>
#include <bamit/IntervalNode.hpp>

#include <seqan3/io/sam_file/input.hpp>

namespace bamit
{

template<typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

template <arithmetic value_type>
struct EstimationResult
{
    value_type mean{0}, median{0}, mode{0}, sd{0}, iqr{0};

    void print()
    {
        seqan3::debug_stream << "Mean: " << mean
                             << "\nMedian: " << median
                             << "\nMode: " << mode
                             << "\nSD: " << sd
                             << "\nIQR: " << iqr << "\n";
    }
};

template <typename traits_type, typename fields_type, typename format_type>
inline EstimationResult<uint64_t> sample_read_depth(seqan3::sam_file_input<traits_type, fields_type, format_type> & input_file,
                                          seqan3::sam_file_header<std::deque<std::string>> & header,
                                          std::vector<std::unique_ptr<IntervalNode>> const & bamit_index,
                                          uint64_t const & sample_value)
    {
        // Intialize vector of read depths at # of positions defined by sample_value.
        std::vector<uint64_t> read_depths(sample_value);
        std::map<uint64_t, uint64_t> depth_counts{};
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr_chr(0, header.ref_ids().size() - 1); // define the range

        uint64_t rand_pos, rand_chr;
        std::tuple<uint64_t, uint64_t> pos_tuple;
        std::streamoff position{};

        for (auto & i : read_depths)
        {
            rand_chr = distr_chr(gen);
            std::uniform_int_distribution<> distr_pos(0, std::get<0>(header.ref_id_info[rand_chr]) - 1);
            rand_pos = distr_pos(gen);
            pos_tuple = std::make_tuple(rand_chr, rand_pos);

            // Get file position of the above location.
            position = get_overlap_records(input_file, bamit_index, pos_tuple, pos_tuple);
            get_correct_position(input_file, pos_tuple, position);

            auto it = input_file.begin();
            it.seek_to(position);
            for (; std::make_tuple((*it).reference_id().value(), (*it).reference_position().value()) <= pos_tuple; ++it)
            {
                ++i;
            }
        }
        std::sort(read_depths.begin(), read_depths.end());
        std::for_each(read_depths.begin(), read_depths.end(), [&depth_counts](int const & value) {++(depth_counts[value]);});

        EstimationResult<uint64_t> result;
        result.mean = std::reduce(std::execution::par, read_depths.begin(), read_depths.end()) / read_depths.size();
        result.median = sample_value % 2 == 0 ? ((read_depths[(sample_value / 2) - 1] + read_depths[(sample_value / 2)]) / 2) :
                                                (read_depths[std::floor(sample_value / 2)]);
        result.mode = (std::max_element(depth_counts.begin(), depth_counts.end(),
                                       [](std::pair<uint64_t, uint64_t> const & a, std::pair<uint64_t, uint64_t> const & b) {return a.second < b.second;}))->second;

        return result;
    }
}
