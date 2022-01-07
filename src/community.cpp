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
    this->g              = new Graph(in_filename, weighted, reproducibility, renumbering);
    this->size           = g->nodes;
    this->neigh_last     = 0;
    this->nb_pass        = nbp;
    this->min_modularity = minm;

    this->neighbor_weight.resize(size, -1.f);
    this->neigh_pos.resize(size);
    this->node_to_community.resize(size); 
    this->communities_arcs.resize(size);

    for (unsigned int i = 0; i < size; ++i) {
        /* the total number of edges inside the community corresponds to 
         * the number of self-loops after contraction
         */
        this->communities_arcs[i].in      = g->count_selfloops(i);
        this->communities_arcs[i].tot_out = g->weighted_out_degree(i);
        this->communities_arcs[i].tot_in  = g->weighted_in_degree(i);
        this->communities_arcs[i].tot     = this->communities_arcs[i].tot_out + this->communities_arcs[i].tot_in;
        // i belongs to its own community
        this->node_to_community[i]        = i;
    }
}

Community::Community(Graph * gc, int nbp, double minm) {
    this->g              = new Graph(*gc);
    this->size           = g->nodes;
    this->nb_pass        = nbp;
    this->min_modularity = minm;
    this->neigh_last     = 0;

    this->neighbor_weight.resize(size, -1.f);
    this->neigh_pos.resize(size);
    this->node_to_community.resize(size); 
    this->communities_arcs.resize(size);

    for (unsigned int i = 0; i < size; ++i) {
        this->communities_arcs[i].in      = g->count_selfloops(i);
        this->communities_arcs[i].tot_out = g->weighted_out_degree(i);
        this->communities_arcs[i].tot_in  = g->weighted_in_degree(i);
        this->communities_arcs[i].tot     = this->communities_arcs[i].tot_out + this->communities_arcs[i].tot_in;
        this->node_to_community[i]        = i; 
    }
}

Community::~Community() {
    delete this->g;
}

/* FIXME: this needs to be tested! 
 * weird: uses neighbor_weight but everyone is at -1 and never changes? */
void Community::init_partition(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in);

    // read partition
    unsigned int node, comm;
    while (finput >> node >> comm) {
        int old_comm = this->node_to_community[node];
        this->neigh_comm(node);

        remove(*this, node, old_comm, neighbor_weight[old_comm]);

        unsigned int best_comm  = 0; 
        double best_nbarcs      = 0.;
        unsigned int i;

        for(i = 0 ; i < size; ++i) {
            best_comm  = neigh_pos[i];
            best_nbarcs      = neighbor_weight[neigh_pos[i]];
            if (best_comm == comm) {
                insert(*this, node, best_comm, best_nbarcs);
                break;
            }
        }

        if (i == neigh_last)
            insert(*this, node, comm, 0.f);
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
    double m = g->get_total_weight();
    for (unsigned int i = 0; i < size; ++i) {
        if (this->communities_arcs[i].tot_in > 0 || this->communities_arcs[i].tot_out > 0) {
            double tot_out_var = this->communities_arcs[i].tot_out / m;
            double tot_in_var = this->communities_arcs[i].tot_in / m;
            q += this->communities_arcs[i].in / m - (tot_out_var * tot_in_var);
        }
    }

    return q;
}

/* FIXME: static function returning neighbor_weight? */
void Community::neigh_comm(unsigned int node) {
    for (unsigned int i = 0; i < neigh_last; ++i)
        neighbor_weight[neigh_pos[i]] = -1.f;

    // at this stage no neighboring community has to be visited
    neigh_last = 0;
    pair< size_t, size_t > p = g->out_neighbors(node);
    unsigned int deg = g->out_degree(node);

    // the first neighboring community of each node is its own
    neigh_pos[0] = node_to_community[node];
    neighbor_weight[neigh_pos[0]] = 0;
    neigh_last = 1;

    for (unsigned int i = 0; i < deg; ++i) {
        // fetching neighbors of i, their community and the corresponding degrees
        unsigned int neigh = g->outcoming_arcs[p.first + i];
        int neigh_comm = node_to_community[neigh];
        double neigh_w = (g->weighted) ? g->outcoming_weights[p.second + i] : 1.f;

        if (neigh != node) {
            // if the community is discovered for the first time
            if (neighbor_weight[neigh_comm] == -1.f) {
                neighbor_weight[neigh_comm] = 0.f;
                neigh_pos[neigh_last++] = neigh_comm;
            }
            // the degree of i towards this community is updated
            neighbor_weight[neigh_comm] += neigh_w;
        }
    }

    // we proceed similarly on in-neighbors
    /* FIXME: don't we count twice the same thing ?! 
     * or do we really need to know w_out+w_in ?*
     */
    pair< size_t, size_t > p_in = g->in_neighbors(node);
    unsigned int deg_in = g->in_degree(node);

    for (unsigned int i = 0; i < deg_in; ++i) {
        unsigned int neigh_in = g->incoming_arcs[p_in.first + i];
        int neigh_comm_in = node_to_community[neigh_in];
        double neigh_w_in = (g->weighted) ? g->incoming_weights[p_in.second + i] : 1.f;

        if (neigh_in != node) {
            if (neighbor_weight[neigh_comm_in] == -1) {
                neighbor_weight[neigh_comm_in] = 0.;
                neigh_pos[neigh_last++] = neigh_comm_in;
            }
            neighbor_weight[neigh_comm_in] += neigh_w_in;
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
        cout << (this->g)->correspondance[i] << " " << renumber[this->node_to_community[i]] << endl;
}

Graph *Community::partition_to_graph() {
    // renumber communities
    vector < int > renumber(size, -1);
    unsigned int f = renumber_communities(*this, renumber);

    // compute communities
    vector < vector < int > > comm_nodes(f);
    for (unsigned int node = 0; node < size; ++node) 
        comm_nodes[renumber[this->node_to_community[node]]].push_back(node);

    // compute weighted graph
    Graph *g2 = new Graph();
    g2->nodes = comm_nodes.size();
    g2->weighted = true;
    for(unsigned int i = 0; i < g2->nodes; ++i)
        g2->correspondance.push_back(i);

    g2->outdegrees.resize(g2->nodes);
    g2->indegrees.resize(g2->nodes);

    unsigned int neigh, neigh_comm;
    double neigh_weight;

    for (size_t comm = 0; comm < g2->nodes; ++comm) {
        map < int, double > m_out, m_in;
        size_t comm_size = comm_nodes[comm].size();
        for (unsigned int node = 0; node < comm_size; ++node) {
            // we first deal with out-neighbors communities
            pair< size_t, size_t > p = (this->g)->out_neighbors(comm_nodes[comm][node]);
            unsigned int deg = (this->g)->out_degree(comm_nodes[comm][node]);
            for (unsigned int i = 0; i < deg; ++i) {
                neigh = (this->g)->outcoming_arcs[p.first + i];
                neigh_comm = renumber[this->node_to_community[neigh]];
                neigh_weight = ((this->g)->weighted) ? (this->g)->outcoming_weights[p.second + i] : 1.f;

                auto it_out = m_out.find(neigh_comm);
                if (it_out == m_out.end())
                    m_out.insert(make_pair(neigh_comm, neigh_weight));
                else
                    it_out -> second += neigh_weight;
            }

            // same thing for in-neighbors communities
            pair< size_t, size_t > p_in = (this->g)->in_neighbors(comm_nodes[comm][node]);
            deg = (this->g)->in_degree(comm_nodes[comm][node]);
            for (unsigned int i = 0; i < deg; ++i) {
                neigh = (this->g)->incoming_arcs[p_in.first + i];
                neigh_comm = renumber[this->node_to_community[neigh]];
                neigh_weight = ((this->g)->weighted) ? (this->g)->incoming_weights[p_in.second + i] : 1.f;

                auto it_in = m_in.find(neigh_comm);
                if (it_in == m_in.end())
                    m_in.insert(make_pair(neigh_comm, neigh_weight));
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

        for (auto it_in = m_in.begin(); it_in != m_in.end(); ++it_in) {
            g2->incoming_arcs.push_back(it_in -> first);
            g2->incoming_weights.push_back(it_in -> second);
        }
    }

    return g2;
}

bool Community::one_level() {
    bool improvement = false;
    int nb_moves = 0;
    int nb_pass_done = 0;
    double new_mod = modularity();
    double cur_mod = new_mod;

    vector < int > random_order(size);
    for (unsigned int i = 0; i < size; ++i)
        random_order[i] = i;

    /* FIXME: uncomment for final release 
     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
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

            // computation of all neighboring communities of current node
            neigh_comm(node);

            // remove node from its current community
            remove(*this, node, node_comm, neighbor_weight[node_comm]);

            // compute the nearest community for node
            // default choice for future insertion is the former community
            int best_comm = node_comm;
            double best_nbarcs = 0.;
            double best_increase = 0.;
            for (unsigned int i = 0; i < neigh_last; ++i) {
                double increase = modularity_gain(*this, node, neigh_pos[i]);
                if (increase > best_increase) {
                    best_comm = neigh_pos[i];
                    best_nbarcs = neighbor_weight[neigh_pos[i]];
                    best_increase = increase;
                }
            }

            // insert node in the nearest community
            insert(*this, node, best_comm, best_nbarcs);

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
    size_t size = c.get_size();
    for (unsigned int node = 0; node < size; ++node) 
        renumber[c.get_node_to_community()[node]]++;
    
    unsigned int f = 0;
    for (unsigned int i = 0; i < size; ++i)
        if (renumber[i] != -1)
            renumber[i] = f++;

    return f;
}

double modularity_gain(const Community &c, unsigned int node, unsigned int comm) {
    assert(node<c.size);
    double dnodecomm = c.neighbor_weight[comm];
    double weighted_out_degree = (c.g)->weighted_out_degree(node);
    double weighted_in_degree = (c.g)->weighted_in_degree(node);
    double totc_out                 = c.communities_arcs[comm].tot_out;
    double totc_in                  = c.communities_arcs[comm].tot_in;
    double m                        = (c.g)->get_total_weight();

    return (dnodecomm / m - ((weighted_out_degree * totc_in + weighted_in_degree * totc_out) / (m*m)));
}

void remove(Community &c, unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<c.size);
    c.communities_arcs[comm].tot_out  -= (c.g)->weighted_out_degree(node);
    c.communities_arcs[comm].tot_in   -= (c.g)->weighted_in_degree(node);
    c.communities_arcs[comm].tot      -= (c.communities_arcs[comm].tot_out + c.communities_arcs[comm].tot_in);
    c.communities_arcs[comm].in       -= dnodecomm + (c.g)->count_selfloops(node);
    c.node_to_community[node]         = -1;
}

void insert(Community &c, unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<c.size);
    c.communities_arcs[comm].tot_out  += (c.g)->weighted_out_degree(node);
    c.communities_arcs[comm].tot_in   += (c.g)->weighted_in_degree(node);
    c.communities_arcs[comm].tot      += (c.communities_arcs[comm].tot_out + c.communities_arcs[comm].tot_in);
    c.communities_arcs[comm].in       += dnodecomm + (c.g)->count_selfloops(node);
    c.node_to_community[node]         = comm;
}

