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

        struct NormalsThreadedWorker
        {
            const std::vector<double>* weights;
            std::vector<double>* normals;
            int start;
            int end;
            int size;
            int radius;
        };
        
        void* get_normals_threaded(void* arg) {
            auto& worker = *static_cast<NormalsThreadedWorker*>(arg);
            auto& weights = *worker.weights;
            auto& normals = *worker.normals;
            
            for (int x = worker.start; x < worker.end; x++) {
                double sum = weights[0];

                for (int wi = 1; wi <= worker.radius; wi++) {
                    if (x >= wi)
                        sum += weights[wi];

                    if (x + wi < worker.size)
                        sum += weights[wi];
                }

                normals[x] = 1.0 / sum;
            }

            return nullptr;
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

    void run_threaded_normals(const std::vector<double>& weights, std::vector<double>& x_normals, std::vector<double>& y_normals, int num_threads, int radius, int x_size, int y_size)
    {
        std::vector<pthread_t> threads(num_threads);
        std::vector<Gauss::NormalsThreadedWorker> workers(num_threads);

        x_normals.resize(x_size);
        y_normals.resize(y_size);

        // x normals
        const auto x_threads = num_threads / 2;
        const auto x_base_size = x_size / x_threads;
        const auto x_remainder = x_size % x_threads;
        
        int start = 0;
        int worker_index = 0;
        
        for (int i = 0; i < x_threads; i++) 
        {
            auto size = x_base_size + ((i < x_remainder) ? 1 : 0);
            auto end = start + size;
            
            auto& args = workers[worker_index];

            args = {
                .weights = &weights,
                .normals = &x_normals,
                .start = start,
                .end = end,
                .size = x_size,
                .radius = radius,
            };

            start = end;
            
            auto& thread = threads[worker_index];  
            pthread_create(&thread, nullptr, Gauss::get_normals_threaded, &args);
            worker_index++;
        }

        start = 0;

        // y normals
        const auto y_threads = num_threads - x_threads;
        const auto y_base_size = y_size / y_threads;
        const auto y_remainder = y_size % y_threads;
        
        for (int i = 0; i < y_threads; i++) 
        {            
            auto size = y_base_size + ((i < y_remainder) ? 1 : 0);
            auto end = start + size;
            
            auto& args = workers[worker_index];
            
            args = {
                .weights = &weights,
                .normals = &y_normals,
                .start = start,
                .end = end,
                .size = y_size,
                .radius = radius,
            };

            start = end;
            
            auto& thread = threads[worker_index];  
            pthread_create(&thread, nullptr, Gauss::get_normals_threaded, &args);
            worker_index++;
        }

        for (const auto& thread : threads)
            pthread_join(thread, nullptr);
    }

    Matrix blur(Matrix m, const int radius, const int num_threads)
    {        
        Matrix scratch{m.get_y_size(), m.get_x_size()};
        auto dst{m};
        
        std::vector<double> w;
        std::vector<double> normals_x;
        std::vector<double> normals_y;

        w = Gauss::get_weights(radius);
        
        if (num_threads < 2) 
        {
            normals_x = Gauss::get_normals(m.get_x_size(), radius, w);
            normals_y = Gauss::get_normals(m.get_y_size(), radius, w);
        }

        else
            run_threaded_normals(w, normals_x, normals_y, num_threads, radius, m.get_x_size(), m.get_y_size());

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
