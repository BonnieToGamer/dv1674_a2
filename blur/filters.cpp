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
        std::vector<double> get_weights(int n)
        {
            std::vector<double>w(n + 1, 0);
            for (int i = 0; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                w[i] = exp(-x * x * pi);
            }

            return w;
        }

        std::vector<double> get_normals(int size, int radius, const std::vector<double>& weights) {
            std::vector<double> normals(size);
            for (int x = 0; x < size; x++) {
                double sum = weights[0];

                for (int wi = 1; wi <= radius; wi++) {
                    if (x >= wi)
                        sum += weights[wi];

                    if (x + wi < size)
                        sum += weights[wi];
                }

                normals[x] = 1.0 / sum;
            }

            return normals;
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{m.get_y_size(), m.get_x_size()};
        auto dst{m};
        
        auto w = Gauss::get_weights(radius);
        auto normals_x = Gauss::get_normals(m.get_x_size(), radius, w);
        auto normals_y = Gauss::get_normals(m.get_y_size(), radius, w);

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

                auto r{w[0] * dst.r(x, y)}, g{w[0] * dst.g(x, y)}, b{w[0] * dst.b(x, y)};
                
                for (int wi = 1; wi <= radius; wi++)
                {
                    double wc{w[wi]};
                    int x2{x - wi};
                    
                    if (x2 >= 0)
                    {
                        const int j = y * size_x + x2;
                        r += wc * dst_r_data[j];
                        g += wc * dst_g_data[j];
                        b += wc * dst_b_data[j];
                    }
                    x2 = x + wi;
                    
                    if (x2 < size_x)
                    {
                        const int j = y * size_x + x2;
                        r += wc * dst_r_data[j];
                        g += wc * dst_g_data[j];
                        b += wc * dst_b_data[j];
                    }
                }

                const auto i = x * size_y + y;
                const auto normal = normals_x[x];
                scratch_r_data[i] = r * normal;
                scratch_g_data[i] = g * normal;
                scratch_b_data[i] = b * normal;
            }
        }

        size_x = scratch.get_x_size();
        size_y = scratch.get_y_size();

        const auto dst_x_size = dst.get_x_size();

        for (int y = 0; y < size_y; y++)
        {
            
            for (int x = 0; x < size_x; x++)
            {
                auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)};
                
                for (int wi = 1; wi <= radius; wi++)
                {
                    double wc{w[wi]};
                    int x2{x - wi};
                    
                    if (x2 >= 0)
                    {
                        const int j = y * size_x + x2;
                        r += wc * scratch_r_data[j];
                        g += wc * scratch_g_data[j];
                        b += wc * scratch_b_data[j];
                    }
                    x2 = x + wi;
                    
                    if (x2 < size_x)
                    {
                        const int j = y * size_x + x2;
                        r += wc * scratch_r_data[j];
                        g += wc * scratch_g_data[j];
                        b += wc * scratch_b_data[j];
                    }
                }
                
                const auto i = x * dst_x_size + y;
                const auto normal = normals_y[x];
                dst_r_data[i] = r * normal;
                dst_g_data[i] = g * normal;
                dst_b_data[i] = b * normal;
            }
        }
        
        return dst;
    }
}
