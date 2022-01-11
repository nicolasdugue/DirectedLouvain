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
    ofstream foutput;

    Community *c = new Community(filename, weighted, precision, reproducibility, renumbering);
    if (filename_part != "")
        c->init_partition(filename_part);
    int level = 0;
    double mod = c->modularity(), new_mod;

    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cerr << "loading graph in: " << fixed << time_taken << " seconds" << endl;
    
    start = chrono::high_resolution_clock::now();
    do {
        const Graph *community_graph = c->get_graph();
        if (verbose) {
            cerr << "level " << level << ":\n";
            cerr << "  network size: " <<
                community_graph->get_nodes() << " nodes, " <<
                community_graph->get_arcs() << " arcs, " <<
                community_graph->get_total_weight() << " weight." << endl;
        }

        improvement = c->one_level();
        new_mod = c->modularity();
        if (++level == display_level)
            community_graph->display();
        if (display_level == -1)
            c->display_partition();
        c->partition_to_graph();
        if (verbose)
            cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

        mod = new_mod;
        if (filename_part != "" && level == 1) // do at least one more computation if partition is provided
            improvement = true;
    } while (improvement);

    delete c;
    cerr << "modularity: " << new_mod << endl;

    end = chrono::high_resolution_clock::now();
    time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

    time_taken *= 1e-9;

    cerr << "computing communities in: " << fixed << time_taken << setprecision(9) << " seconds" << endl;

}
