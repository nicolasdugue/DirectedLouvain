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

// To see the working of controlled
// optimization "Ofast"
#pragma GCC optimize("Ofast")

typedef unsigned int ULI;

using namespace std;

class Graph {
    private:
        friend class Community;
        short type; 

        unsigned int nodes;
        unsigned int arcs;
        double total_weight;  

        vector<unsigned long> outdegrees;
        vector<unsigned long> indegrees;
        vector<unsigned int> outcoming_arcs;
        vector<unsigned int> incoming_arcs;
        vector<double> outcoming_weights, incoming_weights;

        vector<ULI> correspondance;

    public:
        Graph();
        Graph(string in_filename, short type, bool reproducibility, bool renumbering); 
        Graph (const Graph& );

        unsigned int get_nodes() const { return this->nodes; }
        unsigned int get_arcs() const { return this->arcs; }

        double get_total_weight() const { return this->total_weight; }

        vector<ULI> get_correspondance() { return this->correspondance; }

        friend void init_attributes(Graph &g, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN);

        void display() const;
        void display_reverse();
        void load(string filename); 
        void write(string outfile);
        bool check_symmetry();

        void writeFile(string outNeighbors, string inNeighbors);

        // return the number of out neighbors (degree) of the node
        inline unsigned int out_degree(unsigned int node) const;
        // return the number of out neighbors (degree) of the node
        inline unsigned int in_degree(unsigned int node);
        // return the number of self loops of the node
        inline double nb_selfloops(unsigned int node);
        // return the weighted degree of the node
        inline double weighted_out_degree(unsigned int node);
        // return the weighted in-degree of the node
        inline double weighted_in_degree(unsigned int node);
        // return the total degree
        inline double weighted_degree(unsigned int node);
        // return positions of the first out-neighbor and first weight of the node
        inline pair<size_t, size_t > neighbors(unsigned int node) const;
        // return pointers to the first in-neighbor and first weight of the node
        inline pair<size_t, size_t> in_neighbors(unsigned int node);
};


inline unsigned int
Graph::out_degree(unsigned int node) const {
    assert(node<nodes);

    return (node==0 ? outdegrees[0] : outdegrees[node]-outdegrees[node-1]);
}

inline unsigned int
Graph::in_degree(unsigned int node) {
    assert(node<nodes);

    return (node==0 ? indegrees[0] : indegrees[node]-indegrees[node-1]);
}

/* Out-neighbors 
 * TODO: use iterator again ---for optimization sake? 
 */
inline pair<size_t, size_t>
Graph::neighbors(unsigned int node) const {
    assert(node<nodes);

    if (node==0)
        return make_pair(0,0);
    else if (outcoming_weights.size()!=0)
        return make_pair(outdegrees[node-1], outdegrees[node-1]);
    /* FIXME: maybe useless? */
    else
        return make_pair(outdegrees[node-1], 0);
}

/* In-neighbors 
 * TODO: ADAPTER POIDS OUT/IN
 */
inline pair<size_t, size_t>
Graph::in_neighbors(unsigned int node) {
    assert(node<nodes);

    if (node==0)
        return make_pair(0,0);
    else if (incoming_weights.size()!=0)
        return make_pair(indegrees[node-1], indegrees[node-1]);
    /* FIXME: maybe useless? */
    else 
        return make_pair(indegrees[node-1], 0);
}

inline double
Graph::nb_selfloops(unsigned int node) {
    assert(node<nodes);

    pair<size_t, size_t > p = neighbors(node);
    for (double i=0 ; i<out_degree(node) ; ++i) {
        if (outcoming_arcs[p.first+i]==node) {
            if (outcoming_weights.size()!=0)
                return (double)outcoming_weights[p.second+i];
            else 
                return 1.;
        }
    }

    return 0.;
}

inline double
Graph::weighted_out_degree(unsigned int node) {
    assert(node<nodes);
    if (outcoming_weights.size()==0)
        return (double)out_degree(node);
    else {
        pair<size_t, size_t > p = neighbors(node);
        double res = 0;
        for (unsigned int i=0 ; i<out_degree(node) ; ++i) 
            res += (double)outcoming_weights[p.second+i];
        return res;
    }
}

inline double
Graph::weighted_in_degree(unsigned int node) {
    assert(node<nodes);
    if (outcoming_weights.size()==0)
        return (double)in_degree(node);
    else {
        pair<size_t, size_t> p = in_neighbors(node);
        double res = 0;
        for (unsigned int i=0 ; i<in_degree(node) ; ++i) 
            res += (double)incoming_weights[p.second+i];
        return res;
    }
}

inline double
Graph::weighted_degree(unsigned int node) {
    assert(node<nodes);
    return weighted_out_degree(node) + weighted_in_degree(node);

}

#endif // GRAPH_HPP
