/*! \file community.hpp
 *  \brief Header for class Community 
 *         Greedy algorithm optimizing directed modularity
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

#ifndef COMMUNITY_HPP
#define COMMUNITY_HPP

#include "graph.hpp"
#include <memory>

typedef struct count Count;

class Community {
    private:
        Graph* g; // network to compute communities for
        unsigned int size; // nummber of nodes in the network and size of all vectors

        vector<int> node_to_community; // community to which each node belongs

        struct count { 
            double in; /* number of arcs (i.e. self-loops) within the community */
            double tot_in; /* number of outcoming arcs of the community */
            double tot_out; /* number of incoming arcs of the community */
            double tot; /* number of arcs of the community */
            count() : in(0.), tot_in(0.), tot_out(0.), tot(0.) { }
        };
        vector< Count > communities_arcs;

        // number of pass for one level computation
        // if -1, compute as many pass as needed to increase modularity
        /* FIXME: this is _never_ used, must we keep it? */
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
        //void list_neighboring_communities(unsigned int node);

        // compute the modularity of the current partition
        double modularity();

        // displays the current partition (with communities renumbered from 0 to k-1)
        void display_partition();

        // generates the binary graph of communities as computed by one_level
        void partition_to_graph();

        // compute communities of the graph for one level
        // return true if some nodes have been moved
        bool one_level();

        // remove the node from its current community with which it has dnodecomm arcs
        friend void remove(Community&, unsigned int, unsigned int, double);

        // insert the node in comm with which it shares dnodecomm arcs
        friend void insert(Community&, unsigned int, unsigned int, double);

        // compute the gain of modularity if node where inserted in comm
        // given that node has dnodecomm arcs to comm.  The formula is:
        // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
        // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
        // where In(comm)    = number of half-arcs strictly inside comm
        //       Tot(comm)   = number of half-arcs inside or outside comm (sum(degrees))
        //       d(node,com) = number of arcs from node to comm
        //       deg(node)   = node degree
        //       m           = number of arcs
        friend double modularity_gain(const Community&, unsigned int, unsigned int, double);
        friend void list_neighboring_communities(const Community&, vector<double>&, vector<unsigned int>&, unsigned int, unsigned int&);

        inline const Graph *get_graph() {
            return this->g; 
        }
        inline const vector<int>& get_node_to_community() const {
            return this->node_to_community; 
        }
        inline unsigned int get_size() const {
            return this->size; 
        }
        inline unsigned int get_community(unsigned int node) const {
            return this->node_to_community[node];
        }

};

#endif // COMMUNITY_HPP
