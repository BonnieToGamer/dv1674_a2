/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>


int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [threads]" << std::endl;
        std::exit(1);
    }

    auto datasets { Dataset::read(argv[1]) };
    int num_threads = 1;
    if(std::atoi(argv[3]) > 1)
    {
        num_threads = std::atoi(argv[3]);
    }
    
    auto corrs { Analysis::correlation_coefficients(datasets, num_threads) };
    Dataset::write(corrs, argv[2]);

    return 0;
}
