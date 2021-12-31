// File: community.hpp
// -- community detection header file
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

#ifndef COMMUNITY_HPP
#define COMMUNITY_HPP

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

#include "graph.hpp"

class Community {
    private:
        vector<double> neigh_weight;
        vector<unsigned int> neigh_pos;
        unsigned int neigh_last;

        Graph* g; // network to compute communities for
        unsigned int size; // nummber of nodes in the network and size of all vectors
        vector<int> n2c; // community to which each node belongs
        vector<double> in, tot_in, tot_out, tot; // used to compute the modularity participation of each community

        // number of pass for one level computation
        // if -1, compute as many pass as needed to increase modularity
        int nb_pass;

        // a new pass is computed if the last one has generated an increase 
        // greater than min_modularity
        // if 0. even a minor increase is enough to go for one more pass
        double min_modularity;

    public:
        // constructors:
        // reads graph from file using graph constructor
        // type defined the weighted/unweighted status of the graph file
        Community (string in_filename, int type, int nb_pass, double min_modularity, bool reproducibility, bool renumbering);
        Community (Graph* g, int nb_pass, double min_modularity);
        ~Community() { delete g; }

        const Graph *get_graph() { return this->g; }
        const vector<int>& get_n2c() const { return this->n2c; }
        unsigned int get_size() const { return this->size; }

        // initiliazes the partition with something else than all nodes alone
        void init_partition(string filename_part);

        // display the community of each node
        void display();

        // remove the node from its current community with which it has dnodecomm links
        inline void remove(unsigned int node, unsigned int comm, double dnodecomm);

        // insert the node in comm with which it shares dnodecomm links
        inline void insert(unsigned int node, unsigned int comm, double dnodecomm);

        // compute the gain of modularity if node where inserted in comm
        // given that node has dnodecomm links to comm.  The formula is:
        // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
        // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
        // where In(comm)    = number of half-links strictly inside comm
        //       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
        //       d(node,com) = number of links from node to comm
        //       deg(node)   = node degree
        //       m           = number of links
        inline double modularity_gain(unsigned int node, unsigned int comm, double dnodecomm, double w_degree_out, double w_degree_in);

        // compute the set of neighboring communities of node
        // for each community, gives the number of links from node to comm
        void neigh_comm(unsigned int node);

        // compute the modularity of the current partition
        double modularity();

        // displays the graph of communities as computed by one_level
        void partition2graph();
        // displays the current partition (with communities renumbered from 0 to k-1)
        void display_partition();

        // generates the binary graph of communities as computed by one_level
        Graph* partition2graph_binary();

        // compute communities of the graph for one level
        // return true if some nodes have been moved
        bool one_level();
};

inline void
Community::remove(unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<size);

    tot_out[comm]   -= (*g).out_weighted_degree(node);
    tot_in[comm]    -= (*g).in_weighted_degree(node);
    tot[comm]       -= (*g).weighted_degree(node);
    in[comm]        -= dnodecomm + (*g).nb_selfloops(node);
    n2c[node]       = -1;
}

inline void
Community::insert(unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<size);

    tot_out[comm]   += (*g).out_weighted_degree(node);
    tot_in[comm]    += (*g).in_weighted_degree(node);
    tot[comm]       += (*g).weighted_degree(node);
    in[comm]        += dnodecomm + (*g).nb_selfloops(node);
    n2c[node]       = comm;
}

inline double
Community::modularity_gain(unsigned int node, unsigned int comm, double dnodecomm, double w_degree_out, double w_degree_in) {
    assert(node<size);

    double totc_out     = tot_out[comm];
    double totc_in      = tot_in[comm];
    double degc_out     = w_degree_out;
    double degc_in      = w_degree_in;
    double m2           = g->get_total_weight();
    double dnc          = dnodecomm;

    return (dnc/m2 - ((degc_out*totc_in + degc_in*totc_out)/(m2*m2)));
}

#endif // COMMUNITY_HPP
