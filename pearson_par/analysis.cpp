/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <pthread.h>


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
struct threadArgs 
{
    int id;
    int stepSize; 
    int numOfVectors;
    const std::vector<Vector>* vectors;
    double* result;
};

void* threadWorker(void* argument)
{
    auto* args = static_cast<threadArgs*>(argument);

    
    for (auto sample1 = args->id; sample1 < (args->numOfVectors - 1); sample1+= args->stepSize) {

        int index = sample1 * (args->numOfVectors - 1) - sample1 * (sample1 - 1) / 2;
        for (auto sample2 { sample1 + 1 }; sample2 < args->numOfVectors; sample2++) 
        {
            auto corr { pearson((*args->vectors)[sample1], (*args->vectors)[sample2]) };
            args->result[index++] = corr;
        }
    }
    return nullptr;
    
}
std::vector<double> correlation_coefficients(const std::vector<Vector>& datasets, const int& num_threads)
{


    auto normalized_vectors{normalize_vectors(datasets)};
    pthread_t thread[num_threads];
  
    threadArgs args[num_threads];
  
    std::vector<double> result {};
    result.resize(normalized_vectors.size() * (normalized_vectors.size() - 1)/2);

    for (int i = 0; i < num_threads; i++) 
    {
        args[i] = {i, num_threads, (int)normalized_vectors.size(),&normalized_vectors, result.data()};
        pthread_create(&thread[i],nullptr, threadWorker,&args[i]);
    }   

    for (int i = 0; i < num_threads; i++)
    {
        pthread_join(thread[i], nullptr);
        
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
