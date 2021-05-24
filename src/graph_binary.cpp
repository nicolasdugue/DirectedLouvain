// File: graph_binary.cpp
// -- graph handling source
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

#include <sys/mman.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>
#include "../include/graph_binary.h"

Graph::Graph() {
    nb_nodes = 0;
    nb_links_out = 0;
    nb_links_in = 0;
    total_weight = 0;
}

Graph::Graph(char * filename, char * filename_w, int type, bool renumbered) {
    ifstream finput;
    finput.open(filename, fstream:: in | fstream::binary);

    cerr << "number of nodes" << endl;
    // Read number of nodes on 4 bytes
    finput.read((char * ) & nb_nodes, sizeof(int));
    assert(finput.rdstate() == ios::goodbit);
    cerr << "done: " << nb_nodes << endl;

    // Read cumulative out degree sequence: 8 bytes for each node
    cerr << "degrees out" << endl;
    degrees_out.resize(nb_nodes);
    finput.read((char * ) & degrees_out[0], nb_nodes * sizeof(long));
    cerr << "done : " << degrees_out[nb_nodes - 1] << endl;

    // Read links_out: 4 bytes for each link
    nb_links_out = degrees_out[nb_nodes - 1];
    links.resize(nb_links_out);
    finput.read((char * )( & links[0]), (long) nb_links_out * sizeof(unsigned int));

    // Read cumulative in degree sequence: 8 bytes for each node
    cerr << "degrees in" << endl;
    degrees_in.resize(nb_nodes);
    finput.read((char * ) & degrees_in[0], nb_nodes * sizeof(long));
    cerr << "done : " << degrees_in[nb_nodes - 1] << endl;

    // Read links_in: 4 bytes for each link
    nb_links_in = degrees_in[nb_nodes - 1];
    links_in.resize(nb_links_in);
    finput.read((char * )( & links_in[0]), (long) nb_links_in * sizeof(unsigned int));

    // Read correspondance of labels
    if (renumbered) {
        cerr << "correspondance" << endl;
        correspondance.resize(nb_nodes);
        finput.read((char * )( & correspondance[0]), nb_nodes * sizeof(unsigned long long int));
    }

    // if weighted, read weights: 4 bytes for each link (each link is counted twice)
    weights.resize(0);
    weights_in.resize(0);

    total_weight = 0;
    if (type == WEIGHTED) {
        cerr << "Weights reading" << endl;
        ifstream finput_w;
        finput_w.open(filename_w, fstream:: in | fstream::binary);
        weights.resize(nb_links_out);
        finput_w.read((char * ) & weights[0], nb_links_out * sizeof(double));
        weights_in.resize(nb_links_in);
        finput_w.read((char * ) & weights_in[0], nb_links_in * sizeof(double));
        cerr << "Done" << endl;
    }

    // Compute total weight
    for (unsigned int i = 0; i < nb_nodes; i++) {
        total_weight += out_weighted_degree(i);
    }
}

void
Graph::display() {
    for (unsigned int node = 0; node < nb_nodes; node++) {
        pair < vector < unsigned int > ::iterator, vector < double > ::iterator > p = neighbors(node);
        cout << node << ":";
        for (unsigned int i = 0; i < nb_neighbors_out(node); i++) {
            if (true) {
                if (weights.size() != 0)
                    cout << " (" << * (p.first + i) << " " << * (p.second + i) << ")";
                else
                    cout << " " << * (p.first + i);
            }
        }
        cout << endl;
    }
}

void
Graph::writeFile(string outNeighbors, string inNeighbors) {

    ofstream foutput;
    foutput.open(outNeighbors.c_str(), fstream::out | fstream::binary);

    // fetching out-neighbors
    for (unsigned int node = 0; node < nb_nodes; node++) {

        pair < vector < unsigned int > ::iterator, vector < double > ::iterator > p = neighbors(node);
        for (unsigned int i = 0; i < nb_neighbors_out(node); i++) 
            foutput << correspondance[node] << " " << correspondance[ * (p.first + i)] << endl;

    }

    foutput.close();

    ofstream foutputIn;
    foutputIn.open(inNeighbors.c_str(), fstream::out | fstream::binary);

    // fetching in-neighbors
    for (unsigned int node = 0; node < nb_nodes; node++) {

        pair < vector < unsigned int > ::iterator, vector < double > ::iterator > p1 = in_neighbors(node);
        for (unsigned int i = 0; i < nb_neighbors_in(node); i++) 
            foutputIn << correspondance[node] << " " << correspondance[ * (p1.first + i)] << endl;

    }

}

bool
Graph::check_symmetry() {
    int error = 0;
    for (unsigned int node = 0; node < nb_nodes; node++) {
        pair < vector < unsigned int > ::iterator, vector < double > ::iterator > p = neighbors(node);
        for (unsigned int i = 0; i < nb_neighbors_out(node); i++) {
            unsigned int neigh = * (p.first + i);
            double weight = * (p.second + i);

            pair < vector < unsigned int > ::iterator, vector < double > ::iterator > p_neigh = neighbors(neigh);
            for (unsigned int j = 0; j < nb_neighbors_out(neigh); j++) {
                unsigned int neigh_neigh = * (p_neigh.first + j);
                double neigh_weight = * (p_neigh.second + j);

                if (node == neigh_neigh && weight != neigh_weight) {
                    cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
                    if (error++ == 10)
                        exit(0);
                }
            }
        }
    }
    return (error == 0);
}

void
Graph::display_binary(char * outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);

    foutput.write((char * )( & nb_nodes), 4);
    foutput.write((char * )( & degrees_out[0]), 4 * nb_nodes);
    foutput.write((char * )( & links[0]), 8 * nb_links_out);
    foutput.write((char * )( & degrees_in[0]), 4 * nb_nodes);
    foutput.write((char * )( & links_in[0]), 8 * nb_links_in);
}
