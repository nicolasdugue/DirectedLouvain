// File: community.h
// -- community detection source file
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

#include "../include/community.hpp"
#include <climits>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

static unsigned int renumber_communities(const Community &c, vector< int > &renumber);

Community::Community(string in_filename, int weighted, int nbp, double minm, bool reproducibility, bool renumbering) {
    g = new Graph(in_filename, weighted, reproducibility, renumbering);
    size = g->nodes;
    neigh_weight.resize(size, -1);
    neigh_pos.resize(size);
    neigh_last = 0;

    node_to_community.resize(size); 
    communities_arcs.resize(size);

    for (unsigned int i = 0; i < size; ++i) {
        // i belongs to its own community
        node_to_community[i]      = i;
        /* the total number of edges inside the community corresponds to 
         * the number of self-loops after contraction
         */
        communities_arcs[i].in       = g->count_selfloops(i);
        communities_arcs[i].tot_out  = g->weighted_out_degree(i);
        communities_arcs[i].tot_in   = g->weighted_in_degree(i);
        communities_arcs[i].tot   = communities_arcs[i].tot_out + communities_arcs[i].tot_in;
    }

    nb_pass = nbp;
    min_modularity = minm;
}

Community::Community(Graph * gc, int nbp, double minm) {
    g = new Graph(*gc);
    size = g->nodes;

    neigh_weight.resize(size, -1);
    neigh_pos.resize(size);
    neigh_last = 0;

    node_to_community.resize(size); 
    communities_arcs.resize(size);

    for (unsigned int i = 0; i < size; ++i) {
        node_to_community[i]      = i; 
        communities_arcs[i].in       = g->count_selfloops(i);
        communities_arcs[i].tot_out  = g->weighted_out_degree(i);
        communities_arcs[i].tot_in   = g->weighted_in_degree(i);
        communities_arcs[i].tot   = communities_arcs[i].tot_out + communities_arcs[i].tot_in;
    }

    nb_pass = nbp;
    min_modularity = minm;
}

/* FIXME: this needs to be tested! */
void Community::init_partition(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in );

    // read partition
    while (!finput.eof()) {
        unsigned int node, comm;
        finput >> node >> comm;

        if (finput) {
            int old_comm = node_to_community[node];
            neigh_comm(node);

            remove(node, old_comm, neigh_weight[old_comm]);

            unsigned int i = 0;
            for(i = 0 ; i < size; ++i) {
                unsigned int best_comm = neigh_pos[i];
                double best_nbarcs = neigh_weight[neigh_pos[i]];
                if (best_comm == comm) {
                    insert(node, best_comm, best_nbarcs);
                    break;
                }
            }
            if (i == neigh_last)
                insert(node, comm, 0.f);
        }
    }
    finput.close();
}

void Community::display() {
    for (unsigned int i = 0; i < size; ++i)
        cerr << " " << g->correspondance[i] << "/" << node_to_community[i] << "/" 
             << communities_arcs[i].in << "/" << communities_arcs[i].tot;
    cerr << endl;
}

double Community::modularity() {
    double q = 0.;
    double m = g->total_weight;
    for (unsigned int i = 0; i < size; ++i) {
        if (communities_arcs[i].tot_in > 0 || communities_arcs[i].tot_out > 0) {
            double tot_out_var = communities_arcs[i].tot_out / m;
            double tot_in_var = communities_arcs[i].tot_in / m;
            q += communities_arcs[i].in / m - (tot_out_var * tot_in_var);
        }
    }

    return q;
}

void Community::neigh_comm(unsigned int node) {
    for (unsigned int i = 0; i < neigh_last; ++i)
        neigh_weight[neigh_pos[i]] = -1.f;

    // at this stage no neighboring community has to be visited
    neigh_last = 0;
    pair < size_t, size_t > p = g->out_neighbors(node);
    unsigned int deg = g->out_degree(node);

    // the first neighboring community of each node is its own
    neigh_pos[0] = node_to_community[node];
    neigh_weight[neigh_pos[0]] = 0;
    neigh_last = 1;

    for (unsigned int i = 0; i < deg; ++i) {
        // fetching neighbors of i, their community and the corresponding degrees
        unsigned int neigh = g->outcoming_arcs[p.first + i];
        int neigh_comm = node_to_community[neigh];
        double neigh_w = (g->outcoming_weights.size() == 0) ? 1. : g->outcoming_weights[p.second + i];

        if (neigh != node) {
            // if the community is discovered for the first time
            if (neigh_weight[neigh_comm] == -1) {
                neigh_weight[neigh_comm] = 0.f;
                neigh_pos[neigh_last++] = neigh_comm;
            }
            // the degree of i towards this community is updated
            neigh_weight[neigh_comm] += neigh_w;
        }
    }

    // we proceed similarly on in-neighbors
    pair < size_t, size_t > p_in = g->in_neighbors(node);
    unsigned int deg_in = g->in_degree(node);

    for (unsigned int i = 0; i < deg_in; ++i) {

        unsigned int neigh_in = g->incoming_arcs[p_in.first + i];
        int neigh_comm_in = node_to_community[neigh_in];
        double neigh_w_in = (g->incoming_weights.size() == 0) ? 1. : g->incoming_weights[p_in.second + i];

        if (neigh_in != node) {
            if (neigh_weight[neigh_comm_in] == -1) {
                neigh_weight[neigh_comm_in] = 0.;
                neigh_pos[neigh_last++] = neigh_comm_in;
            }
            neigh_weight[neigh_comm_in] += neigh_w_in;
        }
    }
}

/*void Community::partition_to_graph() {
    vector < int > renumber(size, -1);
    renumber_communities(*this, renumber);

    for (unsigned int i = 0; i < size; ++i) {
        pair < size_t, size_t > p = g->out_neighbors(i);

        unsigned int deg = g->out_degree(i);
        for (unsigned int j = 0; j < deg; ++j) {
            unsigned int neigh = g->outcoming_arcs[p.first + j];
            cout << renumber[node_to_community[i]] << " " << renumber[node_to_community[neigh]] << endl;
        }
    }
}*/

void Community::display_partition() {
    vector < int > renumber(size, -1);
    renumber_communities(*this, renumber);

    for (unsigned int i = 0; i < size; ++i)
        cout << g->correspondance[i] << " " << renumber[node_to_community[i]] << endl;
}

Graph *Community::partition_to_graph() {
    // renumber communities
    vector < int > renumber(size, -1);
    unsigned int f = renumber_communities(*this, renumber);

    // compute communities
    vector < vector < int > > comm_nodes(f);
    for (unsigned int node = 0; node < size; ++node) 
        comm_nodes[renumber[node_to_community[node]]].push_back(node);

    // compute weighted graph
    Graph *g2 = new Graph();
    g2->nodes = comm_nodes.size();
    /* FIXME: replaces the ugly trick size() != 0 to check whether the graph is weighted! */
    g2->weighted = true;
    for(unsigned int i = 0; i < g2->nodes; ++i)
        g2->correspondance.push_back(i);

    g2->outdegrees.resize(comm_nodes.size());
    g2->indegrees.resize(comm_nodes.size());

    double neigh_weight;

    size_t comm_deg = comm_nodes.size();
    for (size_t comm = 0; comm < comm_deg; ++comm) {
        map < int, double > m_out, m_in;

        size_t comm_size = comm_nodes[comm].size();
        for (unsigned int node = 0; node < comm_size; ++node) {
            // we first deal with out-neighbors communities
            pair < size_t, size_t > p = g->out_neighbors(comm_nodes[comm][node]);
            unsigned int deg = g->out_degree(comm_nodes[comm][node]);
            for (unsigned int i = 0; i < deg; ++i) {
                unsigned int neigh = g->outcoming_arcs[p.first + i];
                unsigned int neigh_comm = renumber[node_to_community[neigh]];
                neigh_weight = (g->outcoming_weights.size() == 0) ? 1. : g->outcoming_weights[p.second + i];

                auto it_out = m_out.find(neigh_comm);
                if (it_out == m_out.end())
                    m_out.insert(make_pair(neigh_comm, neigh_weight));
                else
                    it_out -> second += neigh_weight;
            }

            // same thing for in-neighbors communities
            pair < size_t, size_t > p_in = g->in_neighbors(comm_nodes[comm][node]);
            deg = g->in_degree(comm_nodes[comm][node]);
            for (unsigned int i = 0; i < deg; ++i) {
                unsigned int neigh = g->incoming_arcs[p_in.first + i];
                unsigned int nc = renumber[node_to_community[neigh]];
                neigh_weight = (g->incoming_weights.size() == 0) ? 1.f : g->incoming_weights[p_in.second + i];

                auto it_in = m_in.find(nc);
                if (it_in == m_in.end())
                    m_in.insert(make_pair(nc, neigh_weight));
                else
                    it_in -> second += neigh_weight;
            }
        }

        g2->outdegrees[comm] = (comm == 0) ? m_out.size() : g2->outdegrees[comm - 1] + m_out.size();
        g2->arcs += m_out.size();

        for (auto it_out = m_out.begin(); it_out != m_out.end(); ++it_out) {
            g2->total_weight += it_out -> second;
            g2->outcoming_arcs.push_back(it_out -> first);
            g2->outcoming_weights.push_back(it_out -> second);
        }

        g2->indegrees[comm] = (comm == 0) ? m_in.size() : g2->indegrees[comm - 1] + m_in.size();
        //g2->arcs += m_in.size();

        for (auto it_in = m_in.begin(); it_in != m_in.end(); ++it_in) {
            g2->incoming_arcs.push_back(it_in -> first);
            g2->incoming_weights.push_back(it_in -> second);
        }
    }

    return g2;
}

bool
Community::one_level() {
    bool improvement = false;
    int nb_moves = 0;
    int nb_pass_done = 0;
    double new_mod = modularity();
    double cur_mod = new_mod;

    vector < int > random_order(size);
    for (unsigned int i = 0; i < size; ++i)
        random_order[i] = i;

    /*unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(random_order.begin(), random_order.end(), std::default_random_engine(seed));*/

    // repeat while
    //   there is an improvement of modularity
    //   or there is an improvement of modularity greater than a given epsilon
    //   or a predefined number of pass have been done
    do {
        cur_mod = new_mod;
        nb_moves = 0;
        nb_pass_done++;

        // for each node: remove the node from its community and insert it in the best community
        for (unsigned int node_tmp = 0; node_tmp < size; ++node_tmp) {
            int node = random_order[node_tmp];
            int node_comm = node_to_community[node];
            double w_degree_out = g->weighted_out_degree(node);
            double w_degree_in = g->weighted_in_degree(node);

            // computation of all neighboring communities of current node
            neigh_comm(node);
            // remove node from its current community
            remove(node, node_comm, neigh_weight[node_comm]);

            // compute the nearest community for node
            // default choice for future insertion is the former community
            int best_comm = node_comm;
            double best_nbarcs = 0.;
            double best_increase = 0.;
            for (unsigned int i = 0; i < neigh_last; ++i) {
                double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree_out, w_degree_in);
                if (increase > best_increase) {
                    best_comm = neigh_pos[i];
                    best_nbarcs = neigh_weight[neigh_pos[i]];
                    best_increase = increase;
                }
            }

            // insert node in the nearest community
            insert(node, best_comm, best_nbarcs);

            if (best_comm != node_comm)
                ++nb_moves;
        }

        new_mod = modularity();

        if (nb_moves > 0)
            improvement = true;

    } while (nb_moves > 0 && new_mod - cur_mod > min_modularity);

    return improvement;
}

static unsigned int renumber_communities(const Community &c, vector< int > &renumber) {
    for (unsigned int node = 0; node < c.get_size(); ++node) 
        renumber[c.get_node_to_community()[node]]++;
    
    unsigned int f = 0;
    for (unsigned int i = 0; i < c.get_size(); ++i)
        if (renumber[i] != -1)
            renumber[i] = f++;

    return f;
}

