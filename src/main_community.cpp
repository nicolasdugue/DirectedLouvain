// File: main_community.cpp
// -- community detection, sample main file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and then Anthony Perez and Nicolas Dugué for directed modularity
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "../include/community.hpp"
#include "../include/utils.hpp" 
#include <chrono>

int main(int argc, char ** argv) {
    auto start = chrono::high_resolution_clock::now();
  
    // unsync the I/O of C and C++.
    /* FIXME: this generates a valgrind error: is this mandatory? */
    //ios_base::sync_with_stdio(false);

    parse_args(argc, argv);
    ofstream foutput;
    foutput.open("modularity_values_directed_louvain.txt", fstream::app | fstream::binary);

    Community *c = new Community(filename, weighted, -1, precision, reproducibility, renumbering);
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
        ++nb_pass;
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
    foutput << new_mod << endl;
    foutput.close();

    end = chrono::high_resolution_clock::now();
    time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

    time_taken *= 1e-9;

    cerr << "computing communities in: " << fixed << time_taken << setprecision(9) << " seconds" << endl;

}
