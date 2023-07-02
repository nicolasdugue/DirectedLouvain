/*! \file main.cpp
 *  \brief Header for class Graph (under [CSR](linktogithub) format)
 *         Base of the Directed Louvain community detection algorithm
 * 
 * ### Based on the articles:
 * + ["Fast unfolding of community hierarchies in large networks"](https://arxiv.org/abs/0803.0476)
 * Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
 * + ["Directed Louvain: maximizing modularity in directed networks"](https://hal.archives-ouvertes.fr/hal-01231784) 
 * N. Dugué, A.Perez 
 * This program must not be distributed without agreement of the above mentionned authors.
 *
 * ### Authors : 
 * + E. Lefebvre, adapted by J.-L. Guillaume 
 * + Adapted by Anthony Perez and Nicolas Dugué for handling directed graphs and modularity
 */

#include "../include/community.hpp"
#include "../include/utils.hpp" 
#include <chrono>

int main(int argc, char ** argv) {
    auto start = chrono::high_resolution_clock::now();
    parse_args(argc, argv);

    // Creating Community object
    Community *c = new Community(filename, precision, gamma, reproducibility, renumbering, randomized);

    auto end = chrono::high_resolution_clock::now();
    // Displaying time (seconds) needed to renumber and load graph
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cerr << "loading graph in: " << fixed << time_taken << " seconds" << endl;
   
    // Computing communities and keeping the number of levels computed
    int levels = c->run(verbose, display_level, filename_part,egc);
    // The last level of the hierarchical structure can be printed using the following
    //c->print_level(levels-1);

    cerr << levels << " levels computed" << endl;
    cerr << "modularity: " << c->modularity() << endl;
    end = chrono::high_resolution_clock::now();
    time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cerr << "computing communities in: " << fixed << time_taken << setprecision(9) << " seconds" << endl;

    // Releasing memory
    delete c;

    return EXIT_SUCCESS;
}
