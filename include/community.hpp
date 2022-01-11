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

typedef struct count Count;

/*! \class Community 
 * \brief   Class implementing the Directed Louvain algorithm. 
 *          Stores information regarding partial graphs that are computed during the algorithm  
 *          that can be displayed if needed (hierarchical community detection)
 */
class Community {
    private:
        Graph* g;                           /*!< A graph object to compute communities for */
        unsigned int size;                  /*!< Number of nodes in the graph */

        vector<int> node_to_community;      /*!< Community to which each node belongs */

        /** 
         * A structure containing information regarding the arcs within all communities 
         */
        struct count { 
            double in;      //<! Number of arcs (i.e. self-loops) within the community */
            double tot_in;  //!< Number of outcoming arcs from the community */
            double tot_out; //!< Number of incoming arcs to the community */
            double tot;     //!< Total number of arcs around the community */
            count() : in(0.), tot_in(0.), tot_out(0.), tot(0.) { }
        };
        vector< Count > communities_arcs;   /*!< A vector of Count structures with arcs information for all communities */

        double precision;                   /*!< A real number describing the minimum improvement on modularity to carry on computation */

    public:
        //! Constructor from edgelist format (initializes Graph object)
        /*! 
         * \param filename          the graph (edgelist format) needed for initializing Graph object attribute
         * \param weighted          boolean value indicating whether the graph is weighted
         * \param reproducibility   boolean value indicating whether to write the renumbered graph on hard drive (readable format)
         * \param renumbering       boolean value indicating whether the graph must be renumbered
         * \sa Graph()
         */
        Community (string filename, bool weighted, double precision, bool reproducibility, bool renumbering);
        //! Destructor
        ~Community(); 

        //! Member function initiliazing first partition with something different than identity
        /*!
         * \param filename the partition to be read
         */
        void init_partition(string filename);

        //! Member function displaying the community of each node
        void display();

        // Member function computing the directed modularity of the current partition
        /*!
         * \return the value of directed modularity for the current partition
         */
        double modularity();

        //! Member function displaying the current partition (with communities renumbered from 0 to k-1) on standard output
        void display_partition();

        //! Member function updating the graph to compute communities for 
        /*!
         * \sa one_level()
         */
        void partition_to_graph();

        //! Member function computing communities of the Graph attribute for one level
        /*!
         * \return true if some nodes have been moved 
         * \sa modularity_gain()
         */
        bool one_level();

        //! Friend method removing a node from its current community with which it has dnodecomm arcs
        /*! 
         * \param c the Community object 
         * \param node the node to remove from a community
         * \param comm the community to remove node from
         * \param dnodecomm the weighted degree of node within its community */
        friend void remove(Community &c, unsigned int node, unsigned int comm, double dnodecomm);

        //! Friend method inserting a node to a new community with which it has dnodecomm arcs
        /*! 
         * \param c the Community object 
         * \param node the node to insert within a community
         * \param comm the community to insert node in
         * \param dnodecomm the weighted degree of node within the community */
        friend void insert(Community &c, unsigned int node, unsigned int comm, double dnodecomm);

        // Friend method computing the gain of modularity if node is inserted into comm
        /*! 
         * Given that node has dnodecomm arcs to comm.  The formula is:
         * [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
         * [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
         * where In(comm)    = number of half-arcs strictly inside comm
         *       Tot(comm)   = number of half-arcs inside or outside comm (sum(degrees))
         *       d(node,com) = weights of arcs from node to comm
         *       deg(node)   = node degree
         *       m           = number of arcs
         *       
         * \param c the Community object
         * \param node the node to consider
         * \param comm the community to consider 
         * \param dnodecomm the weight of arcs from node to comm
         * \return the modularity gained from inserting node into comm (can be a negative value)
         */
        friend double modularity_gain(const Community &c, unsigned int node, unsigned int comm, double dnodecomm);

        //! Friend method computing the set of neighboring communities of a given node
        /*!
         * \param c the Community object
         * \param neighbor_weigh a vector containing, for each community, the total weight of arcs between node and comm 
         * \param neigh_pos a vector representing the communities that are neighbors from node (including its own)
         * \param node the node to consider
         * \return the number of communities neighboring node
         */ 
        friend unsigned int list_neighboring_communities(unsigned int node, const Community &c, vector<double> &neighbor_weight, vector<unsigned int> &neigh_pos):

        //! Getter for the graph to compute communities for
        inline const Graph *get_graph() {
            return this->g; 
        }
        //! Getter for the vector associating nodes to communities
        inline const vector<int>& get_node_to_community() const {
            return this->node_to_community; 
        }
        //! Getter for the size (i.e. number of communities)
        inline unsigned int get_size() const {
            return this->size; 
        }
        //! Getter for the community of a given node
        /*!
         * \param node the node considered
         * \return An integer between 0 and k-1 representing the community to which node belongs
         */
        inline unsigned int get_community(unsigned int node) const {
            return this->node_to_community[node];
        }

};

#endif // COMMUNITY_HPP
