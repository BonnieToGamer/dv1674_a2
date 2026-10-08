/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>

namespace Filter
{
    namespace Gauss
    {
        void get_weights(int n, float *weights_out)
        {
            for (int i = 0; i <= n; i++)
            {
                float x{static_cast<float>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{m.get_x_size(), m.get_y_size()};
        auto dst{m};
        
        std::vector<float> w(radius + 1, 0.0);
        Gauss::get_weights(radius, w.data());

        auto dst_r_data = dst.get_R();
        auto dst_g_data = dst.get_G();
        auto dst_b_data = dst.get_B();
        
        auto scratch_r_data = scratch.get_R();
        auto scratch_g_data = scratch.get_G();
        auto scratch_b_data = scratch.get_B();
        
        auto size_x = dst.get_x_size();
        auto size_y = dst.get_y_size();
        
        for (int y = 0; y < size_y; y++)
        {
            for (int x = 0; x < size_x; x++)
            {
                // unsigned char Matrix::r(unsigned x, unsigned y) const
                // {
                //     return R[y * x_size + x];
                // }

                auto r{w[0] * dst.r(x, y)}, g{w[0] * dst.g(x, y)}, b{w[0] * dst.b(x, y)}, n{w[0]};
                
                for (int wi = 1; wi <= radius; wi++)
                {
                    float wc{w[wi]};
                    int x2{x - wi};
                    
                    if (x2 >= 0)
                    {
                        const int j = y * size_x + x2;
                        r += wc * dst_r_data[j];
                        g += wc * dst_g_data[j];
                        b += wc * dst_b_data[j];
                        n += wc;
                    }
                    x2 = x + wi;
                    
                    if (x2 < size_x)
                    {
                        const int j = y * size_x + x2;
                        r += wc * dst_r_data[j];
                        g += wc * dst_g_data[j];
                        b += wc * dst_b_data[j];
                        n += wc;
                    }
                }

                const auto i = y * size_x + x;
                scratch_r_data[i] = r / n;
                scratch_g_data[i] = g / n;
                scratch_b_data[i] = b / n;
            }
        }

        // transpose for better cache utilization
        scratch.transpose();

        // re-get the pointers
        scratch_r_data = scratch.get_R();
        scratch_g_data = scratch.get_G();
        scratch_b_data = scratch.get_B();

        size_x = scratch.get_x_size();
        size_y = scratch.get_y_size();

        const auto dst_x_size = dst.get_x_size();

        for (int y = 0; y < size_y; y++)
        {
            for (int x = 0; x < size_x; x++)
            {
                auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)}, n{w[0]};

                for (int wi = 1; wi <= radius; wi++)
                {
                    float wc{w[wi]};
                    int x2{x - wi};
                    
                    if (x2 >= 0)
                    {
                        const int j = y * size_x + x2;
                        r += wc * scratch_r_data[j];
                        g += wc * scratch_g_data[j];
                        b += wc * scratch_b_data[j];
                        n += wc;
                    }
                    x2 = x + wi;
                    
                    if (x2 < size_x)
                    {
                        const int j = y * size_x + x2;
                        r += wc * scratch_r_data[j];
                        g += wc * scratch_g_data[j];
                        b += wc * scratch_b_data[j];
                        n += wc;
                    }
                }

                const auto i = x * dst_x_size + y;
                
                dst_r_data[i] = r / n;
                dst_g_data[i] = g / n;
                dst_b_data[i] = b / n;
            }
        }

        return dst;
    }
}
