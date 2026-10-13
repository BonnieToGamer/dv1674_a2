/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include <vector>

#if !defined(FILTERS_HPP)
#define FILTERS_HPP

namespace Filter
{

    namespace Gauss
    {
        constexpr unsigned max_radius{1000};
        constexpr float max_x{1.33};
        constexpr float pi{3.14159};

        std::vector<double> get_weights(int n);
        std::vector<double> get_normals(int size, int radius, const std::vector<double>& weights);
    }

    Matrix blur(Matrix m, const int radius, const int num_threads);

}

#endif