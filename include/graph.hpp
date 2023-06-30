/*! \file graph.hpp
 *  \brief Header for class Graph (under [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format)
 *         Base of the Directed Louvain community detection algorithm
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
 *
 */

#pragma GCC optimize("O3")

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

/*! \class Graph 
 * \brief   Class handling directed graphs. 
 *          Graphs are read from edgelist or binary formats and stored using [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format.
 *
 *          Nodes are renumbered unless stated otherwise, and the resulting correspondance is built.    
 *          Outcoming and incoming arcs are represented using cumulative degree sequences and a list of 
 *          both out- and in-neighbors and the corresponding weights.
 */
class Graph {
    private:
        friend class Community;                 /*!< Community is the main class of the algorithm, and is declared friend for convenience */
        bool weighted;                          /*!< A boolean value indicating whether the graph is weighted */

        unsigned int nodes;                     //!< Number of nodes of the graph
        unsigned int arcs;                      //!< Number of arcs of the graph
        double total_weight;                    //!< Total weight on the arcs of the graph

        vector<unsigned long> outdegrees;       /*!< A vector containing cumulative out-degrees for nodes 0 to nodes-1 */
        vector<unsigned int> outcoming_arcs;    /*!< A vector containing the out-neighbors of every node according to [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format */
        vector<double> outcoming_weights;       /*!< A vector containing the weights of outgoing arcs of every node according to [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format */ 
        vector<unsigned long> indegrees;        /*!< A vector containing cumulative in-degrees for nodes 0 to nodes-1 */
        vector<unsigned int> incoming_arcs;     /*!< A vector containing the in-neighbors of every node according to [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format */
        vector<double> incoming_weights;        /*!< A vector containing the weights of ingoing arcs of every node according to [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format */

        vector<unsigned long> correspondance;   /*!< A vector containing the original label of the input graph (if the graph is not renumbered this is identity) */

    public:
        //! Default constructor
        Graph();
        //! Constructor with arguments
        /*! 
         * \param filename          the file to read the graph from 
         * \param weighted          boolean value indicating whether the graph is weighted
         * \param reproducibility   boolean value indicating whether to write the renumbered graph on hard drive (readable format)
         * \param renumbering       boolean value indicating whether the graph must be renumbered
         * \param verbose           boolean value indicating whether to print information
         */
        Graph(string filename, bool reproducibility, bool renumbering=true, bool weighted=false, bool verbose=false); 
        //! Friend method to initialize all attributes 
        /*!
         * \param g the Graph object to initialize
         * \param LOUT adjacency list for outcoming arcs
         * \param LIN adjacency list for incoming arcs
         * \param verbose           boolean value indicating whether to print information
         * \sa Graph()
         */
        friend void init_attributes(Graph &g, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN, bool verbose);
        //! Copy constructor
        /*! 
         * \param g the Graph object to be copied
         */ 
        Graph (const Graph &g);

        //! Member function loading and initializing Graph object from binary file under [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format
        /*!
         * \param filename path (absolute or relative) to the ".bin" file
         * \param verbose           boolean value indicating whether to print information
         */
        void load(string filename, bool verbose=false); 
        //! Member function writing Graph object into binary file ".bin" under [CSR](https://github.com/nicolasdugue/DirectedLouvain/tree/c+%2B11#CSR) format
        /*! 
         * \param filename path (absolute or relative) to the ".bin" file
         */
        void write(string filename);
        //! Member function printing the Graph object in edgelist format on standard output
        void display() const;

        //! Member function returning the weight (if weighted) or number (else) of self loops of the node
        /*!
         * \param node the node to consider
         * \return the weight or number of self loops of node
         */
        double count_selfloops(unsigned int node);
        //! Member function returning the weighted out-degree of the node
        /*!
         * \param node the node to compute weighted out-degree for
         * \result weighted out-degree of node
         */ 
        double weighted_out_degree(unsigned int node);
        //! Member function returnin the weighted in-degree of the node
        /*!
         * \param node the node to compute weighted in-degree for
         * \result weighted in-degree of node
         */ 
        double weighted_in_degree(unsigned int node);

        //! Member function returning positions of the first out-neighbor of a given node 
        /*! 
         * If the out-degree of the node is 0, this actually returns the position of the first 
         * out-neighbor of node-1. However, this method is always used by looping on the out-degree of node. 
         * \param node the node to consider
         * \result the position of its out-neighbors and weights, according to cumulative out-degree sequence.
         */
        size_t out_neighbors(unsigned int node) const; 
        //! Member function returning the out-degree of a given node 
        /*! 
         * \param node the node to consider
         * \result the out-degree of the node according to cumulative out-degree sequence.
         */
        unsigned int out_degree(unsigned int node) const;
        //! Member function returning positions of the first in-neighbor of a given node 
        /*! 
         * If the in-degree of the node is 0, this actually returns the position of the first 
         * in-neighbor of node-1. However, this method is always used by looping on the out-degree of node. 
         * \param node the node to consider
         * \result the position of its in-neighbors and weights, according to cumulative in-degree sequence.
         */
        size_t in_neighbors(unsigned int node);
        //! Member function returning the in-degree of a given node 
        /*! 
         * \param node the node to consider
         * \result the in-degree of the node according to cumulative in-degree sequence.
         */
        unsigned int in_degree(unsigned int node);
        //! Member function returning the total degree of a given node 
        /*! 
         * \param node the node to consider
         * \result the sum of out- and in-degrees of the node according to cumulative in-degree sequence.
         */
        double weighted_degree(unsigned int node);

        unsigned int get_out_neighbor(size_t index) {
            assert(index < this->arcs);
            return this->outcoming_arcs[index];
        }

        unsigned int get_in_neighbor(size_t index) {
            assert(index < this->arcs);
            return this->incoming_arcs[index];
        }

        unsigned int get_weighted_out_neighbor(size_t index) {
            assert(index < this->arcs);
            return this->outcoming_weights[index];
        }

        unsigned int get_weighted_in_neighbor(size_t index) {
            assert(index < this->arcs);
            return this->incoming_weights[index];
        }

        bool is_weighted() {
            return this->weighted;
        }

        //! Getter for the number of nodes
        unsigned int get_nodes() const {
            return this->nodes; 
        } 
        //! Getter for the number of arcs
        unsigned int get_arcs() const { 
            return this->arcs; 
        }
        //! Getter for the total_weight of the graph
        double get_total_weight() const { 
            return this->total_weight; 
        }
        //! Getter for renumbering correspondance 
        const vector<unsigned long> &get_correspondance() { 
            return this->correspondance; 
        }

};

inline size_t Graph::out_neighbors(unsigned int node) const {
    assert(node<this->nodes);
    return (node==0 ? 0 : this->outdegrees[node-1]);
}

inline unsigned int Graph::out_degree(unsigned int node) const {
    assert(node<this->nodes);
    return (node==0 ? this->outdegrees[0] : this->outdegrees[node]-this->outdegrees[node-1]);
}

inline size_t Graph::in_neighbors(unsigned int node) {
    assert(node<this->nodes);
    return (node==0 ? 0 : this->indegrees[node-1]);
}

inline unsigned int Graph::in_degree(unsigned int node) {
    assert(node<this->nodes);
    return (node==0 ? this->indegrees[0] : this->indegrees[node]-this->indegrees[node-1]);
}

inline double Graph::weighted_degree(unsigned int node) {
    assert(node<this->nodes);
    return this->weighted_out_degree(node) + this->weighted_in_degree(node);
}

#endif // GRAPH_HPP
