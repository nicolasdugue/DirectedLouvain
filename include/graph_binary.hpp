// File: graph_binary.h
// -- graph handling header file
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

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1
#define EMPTY -1 

typedef unsigned long long int ULLI;

using namespace std;

class Graph {
    public:
        unsigned int nb_nodes;
        unsigned int nb_links_out;
        unsigned int nb_links_in;
        float total_weight;  

        vector<unsigned long> degrees_out;
        vector<unsigned long> degrees_in;
        vector<unsigned int> links;
        vector<unsigned int> links_in;
        vector<float> weights, weights_in;

        vector<ULLI> correspondance;

        Graph();
        Graph(string in_filename, int type); 
        Graph(string filename, string filename_w, int type); 
        Graph(int nb_nodes, int nb_links, float total_weight, int *degrees_out, int *links, float *weights);

        void display(void);
        void display_reverse(void);
        void display_binary(char *outfile);
        bool check_symmetry();

        void writeFile(string outNeighbors, string inNeighbors);

        // return the number of out neighbors (degree) of the node
        inline unsigned int nb_neighbors_out(unsigned int node);

        // return the number of out neighbors (degree) of the node
        inline unsigned int nb_neighbors_in(unsigned int node);

        // return the number of self loops of the node
        inline float nb_selfloops(unsigned int node);

        // return the weighted degree of the node
        inline float out_weighted_degree(unsigned int node);

        // return the weighted in-degree of the node
        inline float in_weighted_degree(unsigned int node);

        // return the total degree
        inline float weighted_degree(unsigned int node);

        // return positions of the first out-neighbor and first weight of the node
        inline pair<size_t, size_t > neighbors(unsigned int node);

        // return pointers to the first in-neighbor and first weight of the node
        inline pair<size_t, size_t> in_neighbors(unsigned int node);
};


inline unsigned int
Graph::nb_neighbors_out(unsigned int node) {
    assert(node<nb_nodes);

    if (node==0)
        return degrees_out[0];
    else
        return degrees_out[node]-degrees_out[node-1];
}

inline unsigned int
Graph::nb_neighbors_in(unsigned int node) {
    assert(node<nb_nodes);

    if (node==0)
        return degrees_in[0];
    else
        return degrees_in[node]-degrees_in[node-1];
}

/* Out-neighbors 
 * DONE: SI LE DEGRE NE CHANGE PAS ON RETOURNE 0 
 */
inline pair<size_t, size_t>
Graph::neighbors(unsigned int node) {
    assert(node<nb_nodes);

    if (node==0)
        return make_pair(0,0);
    else if (weights.size()!=0)
        return make_pair(degrees_out[node-1], degrees_out[node-1]);
    else
        return make_pair(degrees_out[node-1], 0);
}

/* In-neighbors 
 * TODO: ADAPTER POIDS OUT/IN
 */
inline pair<size_t, size_t>
Graph::in_neighbors(unsigned int node) {
    assert(node<nb_nodes);

    if (node==0)
        return make_pair(0,0);
    else if (weights_in.size()!=0)
        return make_pair(degrees_in[node-1], degrees_in[node-1]);
    else 
        return make_pair(degrees_in[node-1], 0);
}

inline float
Graph::nb_selfloops(unsigned int node) {
    assert(node<nb_nodes);

    pair<size_t, size_t > p = neighbors(node);
    for (float i=0 ; i<nb_neighbors_out(node) ; ++i) {
        if (links[p.first+i]==node) {
            if (weights.size()!=0)
                return (float)weights[p.second+i];
            else 
                return 1.;
        }
    }

    return 0.;
}


inline float
Graph::out_weighted_degree(unsigned int node) {
    assert(node<nb_nodes);

    if (weights.size()==0)
        return (float)nb_neighbors_out(node);
    else {
        pair<size_t, size_t > p = neighbors(node);
        float res = 0;
        for (unsigned int i=0 ; i<nb_neighbors_out(node) ; i++) {
            res += (float)weights[p.second+i];
        }
        return res;
    }
}

inline float
Graph::in_weighted_degree(unsigned int node) {
    assert(node<nb_nodes);

    if (weights.size()==0)
        return (float)nb_neighbors_in(node);
    else {
        pair<size_t, size_t> p = in_neighbors(node);
        float res = 0;
        for (unsigned int i=0 ; i<nb_neighbors_in(node) ; i++) {
            res += (float)weights_in[p.second+i];
        }
        return res;
    }
}

inline float
Graph::weighted_degree(unsigned int node) {
    assert(node<nb_nodes);

    return out_weighted_degree(node) + in_weighted_degree(node);

}


#endif // GRAPH_HPP
