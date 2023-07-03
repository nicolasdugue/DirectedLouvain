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

/** 
 * A structure containing information regarding the arcs within all communities 
 */
struct count { 
    double total_arcs_inside;       /*<! Number of arcs (i.e. self-loops) within the community */
    double total_incoming_arcs;     /*!< Number of outcoming arcs from the community */
    double total_outcoming_arcs;    /*!< Number of incoming arcs to the community */
    count() : total_arcs_inside(0.), total_incoming_arcs(0.), total_outcoming_arcs(0.) { }
};
typedef struct count Count;

/*! \class Community 
 * \brief   Class implementing the Directed Louvain algorithm. 
 *          Stores information regarding partial graphs that are computed during the algorithm  
 *          that can be displayed if needed (hierarchical community detection)
 */
class Community {
    private:
        Graph* g;                               /*!< A graph object to compute communities for */
        Graph* community_graph;                 /*!< A copy of the original graph object to compute communities for */
        unsigned int size;                      /*!< Number of nodes in the graph */
        double precision;                       /*!< A real number describing the minimum improvement on modularity to carry on computation */
        double gamma;                           /*!< Define the size of generated clusters. Higher gamma means smaller clusters */
        bool randomized;                            /*!< Use to shuffle or not the process order of nodes */

        vector<int> node_to_community;          /*!< Community to which each node belongs */
        vector < vector<int> > levels;          /*! Hierarchical community structure */
        vector< Count > communities_arcs;       /*!< A vector of Count structures with arcs information for all communities */

        void init_attributes();

        //! Private member function initiliazing first partition with something different than identity
        /*!
         * \param filename the partition to be read
         */
        void init_partition(string filename);

        //! Private member function updating the graph to compute communities for 
        /*!
         * \sa one_level()
         */
        void partition_to_graph();

        //! Private member function computing communities of the Graph attribute for one level
        /*!
         * The algorithm proceeds as follow: repeat main procedure while
         * + there is an improvement of modularity
         * + or there is an improvement of modularity greater than a given epsilon
         * + FIXME: the possibility to end after a given number of passes has been removed because  
         *          we never used it. Should we plug it back? (easy to do but...)
         * \param modularity double value containing the new modularity at the end
         * \return true if some nodes have been moved 
         * \sa modularity_gain()
         */
        bool one_level(double &modularity);
        Graph* egc_graph(unsigned int);
        Graph* linear_egc_graph(unsigned int);
    public:
        //! Constructor from edgelist format (initializes Graph object)
        /*! 
         * \param filename          the graph (edgelist format) needed for initializing Graph object attribute
         * \param precision         double value indicating the threshold for modularity improvement
         * \param gamma             double value indicating the size of generated clusters (higher gamma means smaller clusters)
         * \param reproducibility   boolean value indicating whether to write the renumbered graph on hard drive (readable format)
         * \param renumbering       boolean value indicating whether the graph must be renumbered
         * \param randomized        boolean value indicating whether vertices are considered in random order
         * \sa Graph()
         */
        Community (
            const string &filename,
            const double precision=0.0001,
            const double gamma=1,
            bool reproducibility=false,
            bool renumbering=false,
            bool randomized=true);
        //! Destructor
        ~Community(); 

        //! Getter for a specific level
        /*!
        * \param level the level to return
        * \return A map for node associate with their community
        */
        map<unsigned int, unsigned int> get_level(int level);

        //! Getter for the last level of the hierarchical structure
        /*!
        * \return A map for node associate with their community
        */
        map<unsigned int, unsigned int> get_last_level();

        //! Member function displaying the community of each node
        void display();

        //! Member function computing the directed modularity of the current partition
        /*!
         * \return the value of directed modularity for the current partition
         */
        double modularity();

        //! Member function displaying the current partition (with communities renumbered from 0 to k-1) on standard output
        void display_partition();

        //! Member function computing communities of the Graph attribute for one level
        /*!
         * \param verbose       boolean value indicating if verbose mode is activated
         * \param display_level integer value representing the level to display 
         * \param filename_part     an initial partition file (absolute or relative path)
         * \return the number of levels computed by the algorithm
         * The algorithm proceeds while the one_level() function returns true
         * FIXME: the possibility to end after a given number of passes has been removed because  
         *        we never used it. Should we plug it back? (easy to do but...)
         * \sa one_level(), modularity_gain()
         */
        int run(bool verbose, const int& display_level, const string& filename_part, bool egc);

        //! Member function printing a given hierarchical level on standard output
        /*!
         * \param level       the level to print (must be between 0 and levels.size()-1)
         */
        void print_level(int level);

        //! Friend method removing a node from its current community with which it has dnodecomm arcs
        /*! 
         * \param c the Community object 
         * \param node the node to remove from a community
         * \param comm the community to remove node from
         * \param dnodecomm the weighted degree of node within its community */
        friend void remove(Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree);

        //! Friend method inserting a node to a new community with which it has dnodecomm arcs
        /*! 
         * \param c the Community object 
         * \param node the node to insert within a community
         * \param comm the community to insert node in
         * \param dnodecomm the weighted degree of node within the community */
        friend void insert(Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree);

        // Friend method computing the gain of modularity if node is inserted into comm
        /*! 
         * The formula is:
         * d(node,comm)/m - [(dout(node)*In(comm) + din(node)*Out(comm)) / (m*m)]
         * where In(comm)       = number of incoming arcs to comm
         *       Out(comm)      = number of outcoming arcs from comm
         *       dout(node)     = weighted out-degree of node in whole graph
         *       din(node)      = weighted in-degree of node in whole graph
         *       d(node,comm)   = weights of arcs from node to comm
         *       m              = number of arcs
         * \param c the Community object
         * \param node the node to consider
         * \param comm the community to consider 
         * \param dnodecomm the weight of arcs from node to comm
         * \return the modularity gained from inserting node into comm (can be a negative value)
         */
        friend double gain_from_removal(const Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree);
        friend double gain_from_insertion(const Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree);

        //! Friend method computing the set of neighboring communities of a given node
        /*!
         * \param c the Community object
         * \param neighbor_weigh a vector containing, for each community, the total weight of arcs between node and comm 
         * \param neigh_pos a vector representing the communities that are neighbors from node (including its own)
         * \param node the node to consider
         * \param neighboring_communities a reference to the number of communities neighboring node 
         */ 
        friend void list_neighboring_communities(const unsigned int& node, const Community &c, vector<double> &neighbor_weight, vector<unsigned int> &neigh_pos, unsigned int &neighboring_communities);

        //! Getter for the graph to compute communities for
        inline const Graph *get_graph() {
            return this->g; 
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
        //! Getter for the hierarchical community structure
        /*!
         * \return A vector of vector of int containing each level of the hierarchical community structure
         */
        inline const vector< vector<int> > & get_hierarchy() const {
            return this->levels;
        }

};

#endif // COMMUNITY_HPP
