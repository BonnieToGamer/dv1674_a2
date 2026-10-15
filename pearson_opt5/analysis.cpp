/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

namespace Analysis {

// Normalize vectors before sending to pearson - reduces recalculating multiple times for the same vector.
std::vector<Vector> normalize_vectors(const std::vector<Vector>& datasets)
{
    std::vector<Vector> normalized_vectors{};
    normalized_vectors.reserve(datasets.size());

    for(const auto& nVector: datasets)
    {
        auto mm = nVector - nVector.mean();
        normalized_vectors.push_back(mm/ mm.magnitude());
    }
    return normalized_vectors;

}

std::vector<double> correlation_coefficients(const std::vector<Vector>& datasets)
{
    auto normalized_vectors{normalize_vectors(datasets)};

    std::vector<double> result {};
    result.reserve(normalized_vectors.size() * (normalized_vectors.size() - 1)/2);

    for (auto sample1 { 0 }; sample1 < (normalized_vectors.size() - 1); sample1++) {
        for (auto sample2 { sample1 + 1 }; sample2 < normalized_vectors.size(); sample2++) {
            auto corr { pearson(normalized_vectors[sample1], normalized_vectors[sample2]) };
            result.push_back(corr);
        }
    }

    return result;
}

double pearson(const Vector& vec1, const Vector& vec2)
{
    // auto x_mean { vec1.mean() };
    // auto y_mean { vec2.mean() };

    // auto x_mm { vec1 - x_mean };
    // auto y_mm { vec2 - y_mean };

    // auto x_mag { x_mm.magnitude() };
    // auto y_mag { y_mm.magnitude() };

    // auto x_mm_over_x_mag { x_mm / x_mag };
    // auto y_mm_over_y_mag { y_mm / y_mag };

    // auto r { x_mm_over_x_mag.dot(y_mm_over_y_mag) };
    auto r {vec1.dot(vec2)};
    
    return std::max(std::min(r, 1.0), -1.0);
}
};
