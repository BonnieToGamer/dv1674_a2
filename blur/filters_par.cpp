/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters_par.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>
#include <pthread.h>

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

    struct BlurSectionArgs 
    {
        unsigned int start;
        unsigned int end;

        unsigned int size_x;
        unsigned int size_y;

        unsigned char* r_src_data;
        unsigned char* g_src_data;
        unsigned char* b_src_data;

        unsigned char* r_dst_data;
        unsigned char* g_dst_data;
        unsigned char* b_dst_data;

        int radius;
        const std::vector<double>* w;
        const std::vector<double>* normals;
    };
    
    void* blur_section(void* arg) 
    {
        auto& args = *static_cast<BlurSectionArgs*>(arg);

        const auto& w = *args.w;
        const auto& normals = *args.normals;

        for (int i = args.start; i < args.end; i++)
        {
            const int x = i % args.size_x;
            const int y = i / args.size_x;
            
            auto r{w[0] * args.r_src_data[i]}, g{w[0] * args.g_src_data[i]}, b{w[0] * args.b_src_data[i]};
            
            for (int wi = 1; wi <= args.radius; wi++)
            {
                double wc{w[wi]};
                int x2{x - wi};
                
                if (x2 >= 0)
                {
                    const int j = y * args.size_x + x2;
                    r += wc * args.r_src_data[j];
                    g += wc * args.g_src_data[j];
                    b += wc * args.b_src_data[j];
                }
                x2 = x + wi;
                
                if (x2 < args.size_x)
                {
                    const int j = y * args.size_x + x2;
                    r += wc * args.r_src_data[j];
                    g += wc * args.g_src_data[j];
                    b += wc * args.b_src_data[j];
                }
            }

            const auto j = x * args.size_y + y;
            
            const auto normal = normals[x];
            args.r_dst_data[j] = r * normal;
            args.g_dst_data[j] = g * normal;
            args.b_dst_data[j] = b * normal;
        }

        return nullptr;
    }

    struct ThreadWorker 
    {
        unsigned int id;
        unsigned int start;
        unsigned int end;
        pthread_t thread;
    };

    std::vector<ThreadWorker> create_workers(const int num_workers, unsigned int total_size) 
    {
        const auto base_size = total_size / num_workers;
        const auto remainder = total_size % num_workers;

        std::vector<ThreadWorker> workers;
        int start = 0;
        
        for (int i = 0; i < num_workers; i++) 
        {
            auto size = base_size + ((i < remainder) ? 1 : 0);
            auto end = start + size;

            workers.emplace_back(ThreadWorker{
                .id = static_cast<unsigned int>(i),
                .start = static_cast<unsigned int>(start),
                .end = static_cast<unsigned int>(end),
            });

            start = end;
        }

        return workers;
    }

    Matrix blur(Matrix m, const int radius, const int num_threads)
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
        
        auto workers = create_workers(num_threads, size_x * size_y);
        std::vector<BlurSectionArgs> args_list(workers.size());
        
        for (auto& worker : workers) 
        {
            auto& args = args_list[worker.id];

            args = {
                .start = worker.start,
                .end = worker.end,

                .size_x = size_x,
                .size_y = size_y,

                .r_src_data = dst_r_data,
                .g_src_data = dst_g_data,
                .b_src_data = dst_b_data,

                .r_dst_data = scratch_r_data,
                .g_dst_data = scratch_g_data,
                .b_dst_data = scratch_b_data,

                .radius = radius,
                .w = &w,
                .normals = &normals_x,
            };
            
            pthread_create(&worker.thread, nullptr, blur_section, &args);
        }

        for (auto& worker : workers) {
            pthread_join(worker.thread, nullptr);
        }

        size_x = scratch.get_x_size();
        size_y = scratch.get_y_size();

        const auto dst_x_size = dst.get_x_size();

        for (auto& worker : workers) 
        {
            auto& args = args_list[worker.id];

            args = {
                .start = worker.start,
                .end = worker.end,

                .size_x = size_x,
                .size_y = size_y,

                .r_src_data = scratch_r_data,
                .g_src_data = scratch_g_data,
                .b_src_data = scratch_b_data,

                .r_dst_data = dst_r_data,
                .g_dst_data = dst_g_data,
                .b_dst_data = dst_b_data,

                .radius = radius,
                .w = &w,
                .normals = &normals_y,
            };
            
            pthread_create(&worker.thread, nullptr, blur_section, &args);
        }

        for (auto& worker : workers) {
            pthread_join(worker.thread, nullptr);
        }
        
        return dst;
    }
}
