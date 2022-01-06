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

#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include "graph.hpp"

typedef struct count count;

class Community {
    private:
        Graph* g; // network to compute communities for
        unsigned int size; // nummber of nodes in the network and size of all vectors

        /* FIXME: are those really attributes of the Community class? 
         * They seem to indicate, for each community, the best neighbor they have, right? 
         */
        vector<double> neigh_weight;
        vector<unsigned int> neigh_pos;
        unsigned int neigh_last;

        vector<int> node_to_community; // community to which each node belongs

        struct count { 
            double in; /* number of arcs (i.e. self-loops) within the community */
            double tot_in; /* number of outcoming arcs of the community */
            double tot_out; /* number of incoming arcs of the community */
            double tot; /* number of arcs of the community */
        };

        vector<count> communities_arcs;

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
        ~Community(); 

        // initiliazes the partition with something else than all nodes alone
        void init_partition(string filename_part);

        // display the community of each node
        void display();

        // compute the set of neighboring communities of node
        // for each community, gives the number of arcs from node to comm
        void neigh_comm(unsigned int node);

        // compute the modularity of the current partition
        double modularity();

        /* FIXME: commented for the moment: is it useful??? */
        // displays the graph of communities as computed by one_level
        //void partition_to_graph();
        // displays the current partition (with communities renumbered from 0 to k-1)
        void display_partition();

        // generates the binary graph of communities as computed by one_level
        Graph* partition_to_graph();

        // compute communities of the graph for one level
        // return true if some nodes have been moved
        bool one_level();
        inline const Graph *get_graph() {
            return this->g; 
        }
        inline const vector<int>& get_node_to_community() const {
            return this->node_to_community; 
        }
        inline unsigned int get_size() const {
            return this->size; 
        }

        /*inline void remove(unsigned int node, unsigned int comm, double dnodecomm);

        // insert the node in comm with which it shares dnodecomm arcs
        inline void insert(unsigned int node, unsigned int comm, double dnodecomm);

        // compute the gain of modularity if node where inserted in comm
        // given that node has dnodecomm arcs to comm.  The formula is:
        // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
        // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
        // where In(comm)    = number of half-arcs strictly inside comm
        //       Tot(comm)   = number of half-arcs inside or outside comm (sum(degrees))
        //       d(node,com) = number of arcs from node to comm
        //       deg(node)   = node degree
        //       m           = number of arcs
        inline double modularity_gain(unsigned int node, unsigned int comm, double dnodecomm, double w_degree_out, double w_degree_in);*/

        // remove the node from its current community with which it has dnodecomm arcs
        friend void remove(Community &c, unsigned int node, unsigned int comm, double dnodecomm);

        // insert the node in comm with which it shares dnodecomm arcs
        friend void insert(Community &c, unsigned int node, unsigned int comm, double dnodecomm);

        // compute the gain of modularity if node where inserted in comm
        // given that node has dnodecomm arcs to comm.  The formula is:
        // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
        // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
        // where In(comm)    = number of half-arcs strictly inside comm
        //       Tot(comm)   = number of half-arcs inside or outside comm (sum(degrees))
        //       d(node,com) = number of arcs from node to comm
        //       deg(node)   = node degree
        //       m           = number of arcs
        friend double modularity_gain(const Community &c, unsigned int node, unsigned int comm, double dnodecomm, double w_degree_out, double w_degree_in);
};

/*inline double Community::modularity_gain(unsigned int node, unsigned int comm, double dnodecomm, double w_degree_out, double w_degree_in) {
    assert(node<size);
    double totc_out                 = communities_arcs[comm].tot_out;
    double totc_in                  = communities_arcs[comm].tot_in;
    double m                        = (g)->get_total_weight();

    return (dnodecomm / m - ((w_degree_out * totc_in + w_degree_in * totc_out) / (m*m)));
}

inline void Community::remove(unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<size);
    communities_arcs[comm].tot_out  -= (g)->weighted_out_degree(node);
    communities_arcs[comm].tot_in   -= (g)->weighted_in_degree(node);
    communities_arcs[comm].tot      -= (communities_arcs[comm].tot_out + communities_arcs[comm].tot_in);
    communities_arcs[comm].in       -= dnodecomm + (g)->count_selfloops(node);
    node_to_community[node]         = -1;
}

inline void Community::insert(unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<size);
    communities_arcs[comm].tot_out  += (g)->weighted_out_degree(node);
    communities_arcs[comm].tot_in   += (g)->weighted_in_degree(node);
    communities_arcs[comm].tot      += (communities_arcs[comm].tot_out + communities_arcs[comm].tot_in);
    communities_arcs[comm].in       += dnodecomm + (g)->count_selfloops(node);
    node_to_community[node]         = comm;
}*/
#endif // COMMUNITY_HPP
