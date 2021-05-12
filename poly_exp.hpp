
///////////////////////////////////////////////////////////////////////////////
// poly_exp.hpp
//
// Definitions for two algorithms that solve the Maximum Subarray Problem,
// and one algorithm that solves the Subset Sum Problem.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <functional>
#include <optional>
#include <vector>

namespace subarray {

    // A summed_span represents a non-empty range of indices inside of a vector of
    // ints, stored in a begin iterator and end iterator. The class also stores
    // the sum of the ints in that range.
    //
    // Just like in the rest of the C++ standard library, the range includes all
    // elements in [begin, end), or in other words the range includes begin, and all
    // elements up to BUT NOT INCLUDING end itself.
    class summed_span {
    public:
        using iterator = std::vector<int>::const_iterator;

    private:
        iterator begin_, end_;
        int sum_;

    public:

        // Constructor, given the begin iterator, end iterator, and sum of elements
        // in the range. begin must come before end. sum must be the total of all
        // elements in the range. O(1) time.
        summed_span(iterator begin, iterator end, int sum)
            : begin_(begin), end_(end), sum_(sum) {
            assert(begin < end);
        }

        // Constructor, given only the begin and end iterators. The sum is computed
        // in O(n) time.
        summed_span(iterator begin, iterator end)
            : summed_span(begin, end, std::accumulate(begin, end, 0)) {}

        // Equality tests, two spans are equal when each of their iterators are equal.
        bool operator== (const summed_span& rhs) const {
            return (begin_ == rhs.begin_) && (end_ == rhs.end_);
        }

        // Accessors.
        const iterator& begin() const { return begin_; }
        const iterator& end() const { return end_; }
        int sum() const { return sum_; }

        // Compute the number of elements in the span.
        size_t size() const { return end_ - begin_; }

        // Stream insertion operator, so this class is printable.
        friend std::ostream& operator<<(std::ostream& stream, const summed_span& rhs) {
            stream << "summed_span, size=" << rhs.size() << ", sum=" << rhs.sum();
            return stream;
        }
    };

    // Compute the maximum subarray of input; i.e. the non-empty contiguous span of
    // elements with the maximum sum. input must be nonempty. This function uses an
    // exhaustive search algorithm that takes O(n^3) time.
    summed_span max_subarray_exh(const std::vector<int>& input) {

        assert(!input.empty());

        int b = 0;
        int e = 1;
        int n = input.size();

        for (size_t i = 0; i <= n - 1; i++)
        {
            for (int j = i + 1; j <= n; j++)
            {
                int sumItoJ = 0;
                int sumBtoE = 0;
                for (int k = i; k < j; k++)
                {
                    sumItoJ += input[k];
                }
                for (int l = b; l < e; l++)
                {
                    sumBtoE += input[l];
                }
                if (sumItoJ > sumBtoE)
                {
                    b = i;
                    e = j;
                }
            }
        }
        if (input.begin() + e == input.end())
        {
            return summed_span(input.begin() + b, input.end());
        }
        else
            return summed_span(input.begin() + b, input.begin() + e);
    }

    // Compute the maximum subarray using a decrease-by-half algorithm that takes
    // O(n log n) time.

    //
    summed_span crossing_max_subarray(const std::vector<int>& input, int low, int middle, int high)
    {
        int sumL = -99999999;
        int sumR = sumL;
        int trueSum = 0;
        int b = 0;
        int e = 0;

        for (int i = middle; i >= low; i--)
        {
            trueSum += input[i];
            if (trueSum > sumL)
            {
                sumL = trueSum;
                b = i;
            }
        }
        trueSum = 0;
        for (int j = middle + 1; j <= high; j++)
        {
            trueSum += input[j];
            if (trueSum > sumR)
            {
                sumR = trueSum;
                e = j;
            }
        }
        if (input.begin() + e + 1 == input.end())
        {
            return summed_span(input.begin() + b, input.end());
        }
        else
        {
            return summed_span(input.begin() + b, input.begin() + e + 1);
        }
    }

    //Using a recursive helper function to do decrease-by-half
    summed_span maximum_subarray_recurse(const std::vector<int>& input, int low, int high)
    {
        if (low == high)
        {
            return summed_span(input.begin() + low, input.begin() + low + 1);
        }

        size_t avg = (low + high) / 2;

        summed_span onLeft = maximum_subarray_recurse(input, low, avg);
        summed_span onRight = maximum_subarray_recurse(input, avg + 1, high);
        summed_span crossing = crossing_max_subarray(input, low, avg, high);

        if ((crossing.sum() > onLeft.sum()) && (crossing.sum() > onRight.sum()))
        {
            return crossing;
        }
        if ((onLeft.sum() > onRight.sum()) && (onLeft.sum() > crossing.sum()))
        {
            return onLeft;
        }
        else
        {
            return onRight;
        }

    }
    summed_span max_subarray_dbh(const std::vector<int>& input) {

        assert(!input.empty());
        return maximum_subarray_recurse(input, 0, input.size() - 1);
        //return summed_span(input.begin(), input.end());
    }

    // Solve the subset sum problem: return a non-empty subset of input that adds
    // up to exactly target. If no such subset exists, return an empty optional.
    // input must not be empty, and must contain fewer than 64 elements.
    // Note that the returned subset must never be empty, even if target == 0.
    // This uses an exhaustive search algorithm that takes exponential O(n * 2^n)
    // time.
    std::optional<std::vector<int>>
        subset_sum_exh(const std::vector<int>& input, int target) {

        assert(!input.empty());
        assert(input.size() < 64);

        uint64_t largestInt = 1 << input.size();

        for (uint64_t i = 0; i < largestInt; i++)
        {
            std::vector<int> candidate = {};
            for (int j = 0; j < input.size(); j++)
            {
                if ((i >> j & 1) == 1)
                {
                    candidate.emplace_back(input[j]);
                }
                int x = 0;
                for (int k = 0; k < candidate.size(); k++)
                {
                    x += candidate[k];
                }
                if (candidate.size() > 0 && x == target)
                {
                    return candidate;
                }
            }
        }

        return std::nullopt;
    }

}
