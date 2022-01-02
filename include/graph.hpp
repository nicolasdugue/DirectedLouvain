// File: graph.hpp
// -- graph handling header file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Authors   : E. Lefebvre, adapted by J.-L. Guillaume and then Anthony Perez and Nicolas Dugué for directed modularity
//-----------------------------------------------------------------------------

#pragma GCC optimize("O2")

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

typedef unsigned int ULI;

class Graph {
    private:
        friend class Community;
        bool weighted; 

        unsigned int nodes;
        unsigned int arcs;
        double total_weight;  

        vector<unsigned long> outdegrees;
        vector<unsigned long> indegrees;
        vector<unsigned int> outcoming_arcs;
        vector<unsigned int> incoming_arcs;
        vector<double> outcoming_weights; 
        vector<double> incoming_weights;

        vector<ULI> correspondance;

    public:
        Graph();
        Graph(string in_filename, bool weighted, bool reproducibility, bool renumbering); 
        Graph (const Graph& );
        friend void init_attributes(Graph &g, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN);

        unsigned int get_nodes() const { return this->nodes; }
        unsigned int get_arcs() const { return this->arcs; }
        double get_total_weight() const { return this->total_weight; }

        const vector<ULI> &get_correspondance() { return this->correspondance; }

        void display() const;
        void display_reverse();
        void load(string filename); 
        void write(string outfile);

        // return the number of self loops of the node
        double count_selfloops(unsigned int node);
        // return the weighted degree of the node
        double weighted_out_degree(unsigned int node);
        // return the weighted in-degree of the node
        double weighted_in_degree(unsigned int node);

        // return positions of the first out-neighbor and first weight of the node
        pair<size_t, size_t > out_neighbors(unsigned int node) const; 
        // return the number of out neighbors (degree) of the node
        unsigned int out_degree(unsigned int node) const;
        // return pointers to the first in-neighbor and first weight of the node
        pair<size_t, size_t> in_neighbors(unsigned int node);
        // return the number of out neighbors (degree) of the node
        unsigned int in_degree(unsigned int node);
        // return the total degree
        double weighted_degree(unsigned int node);

};

// return positions of the first out-neighbor and first weight of the node
inline pair<size_t, size_t > Graph::out_neighbors(unsigned int node) const {
    assert(node<this->nodes);
    if (node==0)
        return make_pair(0,0);
    else if (this->weighted)
        return make_pair(this->outdegrees[node-1], this->outdegrees[node-1]);
    /* FIXME: maybe useless? second member is probably used only if weighted==WEIGHTED*/
    else
        return make_pair(this->outdegrees[node-1], 0);
}

// return the number of out neighbors (degree) of the node
inline unsigned int Graph::out_degree(unsigned int node) const {
    assert(node<this->nodes);
    return (node==0 ? this->outdegrees[0] : this->outdegrees[node]-this->outdegrees[node-1]);
}

// return pointers to the first in-neighbor and first weight of the node
inline pair<size_t, size_t> Graph::in_neighbors(unsigned int node) {
    assert(node<this->nodes);
    if (node==0)
        return make_pair(0,0);
    else if (this->weighted)
        return make_pair(this->indegrees[node-1], this->indegrees[node-1]);
    /* FIXME: maybe useless? */
    else 
        return make_pair(this->indegrees[node-1], 0);
}

// return the number of out neighbors (degree) of the node
inline unsigned int Graph::in_degree(unsigned int node) {
    assert(node<this->nodes);
    return (node==0 ? this->indegrees[0] : this->indegrees[node]-this->indegrees[node-1]);
}

// return the total degree
inline double Graph::weighted_degree(unsigned int node) {
    assert(node<this->nodes);
    return this->weighted_out_degree(node) + this->weighted_in_degree(node);
}

#endif // GRAPH_HPP
