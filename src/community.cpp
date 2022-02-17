#include "../include/community.hpp"
#include <climits>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

// Static function renumbering communities from 0 to k-1 (returns k)
static unsigned int renumber_communities(const Community &c, vector< int > &renumber);

Community::Community(const string& in_filename, bool weighted, const double& minm, bool reproducibility, bool renumbering) {
    this->g              = new Graph(in_filename, weighted, reproducibility, renumbering);
    this->size           = g->nodes;
    this->precision = minm;

    this->node_to_community.resize(this->size); 
    this->communities_arcs.resize(this->size);

    for (unsigned int i = 0; i < this->size; ++i) {
        // i belongs to its own community
        this->node_to_community[i]                      = i;
        // the total number of edges inside the community corresponds to 
        // the number of self-loops of node i (after at least one step)
        this->communities_arcs[i].total_arcs_inside     = g->count_selfloops(i);
        this->communities_arcs[i].total_outcoming_arcs  = g->weighted_out_degree(i);
        this->communities_arcs[i].total_incoming_arcs   = g->weighted_in_degree(i);
    }
}

Community::~Community() {
    delete this->g;
}


/* FIXME: this needs to be tested! */
void Community::init_partition(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in);
    assert(finput.rdstate() == ios::goodbit);

    unsigned int node, comm;
    while (finput >> node >> comm) {
        vector<unsigned int> positions_neighboring_communities(this->size);
        vector<double> neighbor_weight(this->size, -1);

        int old_comm = this->node_to_community[node];
        unsigned int neighboring_communities = 0;
        list_neighboring_communities(node, *this, neighbor_weight, positions_neighboring_communities, neighboring_communities);

        remove(*this, node, old_comm, neighbor_weight[old_comm], (this->g)->weighted_out_degree(node), (this->g)->weighted_in_degree(node));

        unsigned int best_community  = 0; 
        double best_nbarcs      = 0.;
        unsigned int i;

        for(i = 0 ; i < size; ++i) {
            best_community   = positions_neighboring_communities[i];
            best_nbarcs = neighbor_weight[positions_neighboring_communities[i]];
            if (best_community == comm) {
                insert(*this, node, best_community, best_nbarcs, (this->g)->weighted_out_degree(node), (this->g)->weighted_in_degree(node));
                break;
            }
        }

        if (i == neighboring_communities)
            insert(*this, node, comm, 0.f, (this->g)->weighted_out_degree(node), (this->g)->weighted_in_degree(node));
    }
    finput.close();
}

void Community::display() {
    for (unsigned int i = 0; i < size; ++i)
        cout << " " << g->correspondance[i] << "/" << node_to_community[i] << "/" 
             << this->communities_arcs[i].total_arcs_inside << "/" 
             << this->communities_arcs[i].total_outcoming_arcs << "/" << this->communities_arcs[i].total_incoming_arcs;
    cout << endl;
}

double Community::modularity() {
    double q = 0.;
    double m = g->get_total_weight();
    for (unsigned int i = 0; i < size; ++i) {
        if (this->communities_arcs[i].total_incoming_arcs > 0 || this->communities_arcs[i].total_outcoming_arcs > 0) {
            double total_outcoming_arcs_var     = this->communities_arcs[i].total_outcoming_arcs / m;
            double total_incoming_arcs_var      = this->communities_arcs[i].total_incoming_arcs / m;
            q                                   += this->communities_arcs[i].total_arcs_inside / m - (total_outcoming_arcs_var * total_incoming_arcs_var);
        }
    }
    return q;
}

void Community::display_partition() {
    vector < int > renumber(size, -1);
    renumber_communities(*this, renumber);
    // Marking the beginning of new level
    cout << -1 << " " << -1 << endl;
    for (unsigned int i = 0; i < size; ++i)
        cout << (this->g)->correspondance[i] << " " << renumber[this->node_to_community[i]] << endl;
}

void Community::partition_to_graph() {
    // Renumbering communities
    vector<int> renumber(size, -1);
    unsigned int f = renumber_communities(*this, renumber);

    // Computing communities (k lists of nodes)
    vector <vector<int>> comm_nodes(f);
    for (unsigned int node = 0; node < size; ++node) 
        comm_nodes[renumber[this->node_to_community[node]]].push_back(node);

    // Computing contracted weighted graph
    Graph *g2 = new Graph();
    g2->nodes = comm_nodes.size();
    g2->weighted = true;

    // Correspondance is set to identity since the graph is renumbered from 0 to k-1 (communities)
    for(unsigned int i = 0; i < g2->nodes; ++i)
        g2->correspondance.push_back(i);

    g2->outdegrees.resize(g2->nodes);
    g2->indegrees.resize(g2->nodes);

    unsigned int out_neighbor, out_neighboring_community;
    double out_weight;
    unsigned int in_neighbor, in_neighboring_community;
    double in_weight;

    // Computing arcs between communities
    for (size_t comm = 0; comm < g2->nodes; ++comm) {
        map <int,double> m_out, m_in;
        size_t comm_size = comm_nodes[comm].size();
        for (unsigned int node = 0; node < comm_size; ++node) {
            // Out-neighbors
            size_t p = (this->g)->out_neighbors(comm_nodes[comm][node]);
            unsigned int deg = (this->g)->out_degree(comm_nodes[comm][node]);
            // Looking for communities of every out-neighbor of node and then storing/updating weighted out-degrees
            for (unsigned int i = 0; i < deg; ++i) {
                out_neighbor = (this->g)->outcoming_arcs[p + i];
                out_neighboring_community = renumber[this->node_to_community[out_neighbor]];
                out_weight = ((this->g)->weighted) ? (this->g)->outcoming_weights[p + i] : 1.f;

                auto it_out = m_out.find(out_neighboring_community);
                if (it_out == m_out.end())
                    m_out.insert(make_pair(out_neighboring_community, out_weight));
                else
                    it_out -> second += out_weight;
            }

            // In-neighbors
            size_t p_in = (this->g)->in_neighbors(comm_nodes[comm][node]);
            deg = (this->g)->in_degree(comm_nodes[comm][node]);
            // Looking for communities of every in-neighbor of node and then storing/updating weighted in-degrees
            for (unsigned int i = 0; i < deg; ++i) {
                in_neighbor = (this->g)->incoming_arcs[p_in + i];
                in_neighboring_community = renumber[this->node_to_community[in_neighbor]];
                in_weight = ((this->g)->weighted) ? (this->g)->incoming_weights[p_in + i] : 1.f;

                auto it_in = m_in.find(in_neighboring_community);
                if (it_in == m_in.end())
                    m_in.insert(make_pair(in_neighboring_community, in_weight));
                else
                    it_in -> second += in_weight;
            }
        }

        // Building outcoming and incoming arcs according to previously computed weights
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

    // Updating graph attribute with computed graph g2
    delete this->g;
    this->g = g2;

    // Updating other attributes according to constructed graph g
    this->size           = this->g->nodes;
    this->node_to_community.resize(this->size); 

    for (unsigned int i = 0; i < this->size; ++i) {
        this->node_to_community[i]                      = i; 
        this->communities_arcs[i].total_arcs_inside     = this->g->count_selfloops(i);
        this->communities_arcs[i].total_outcoming_arcs  = this->g->weighted_out_degree(i);
        this->communities_arcs[i].total_incoming_arcs   = this->g->weighted_in_degree(i);
    }
}

bool Community::one_level(double &modularity) {
    int nb_moves = 0;
    bool improvement = false;
    double current_modularity = this->modularity();
    double delta;

    // Order in which to proceed nodes of the graph...
    vector < int > random_order(size);
    for (unsigned int i = 0; i < size; ++i)
        random_order[i] = i;

    // ... randomized: (Directed) Louvain's algorithm is not deterministic
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(random_order.begin(), random_order.end(), std::default_random_engine(seed));

    // Vectors containing weights and positions of neighboring communities
    vector<double> neighbor_weight(this->size,-1);
    vector<unsigned int> positions_neighboring_communities(size);
    // Every node neighbors its own community
    unsigned int neighboring_communities = 1;
    double total_increase = 0.;

    do {
        nb_moves = 0;
        total_increase = 0.;

        // For each node: remove it from its community and insert it in the best community (if any)
        for (unsigned int node_tmp = 0; node_tmp < size; ++node_tmp) {
            int node = random_order[node_tmp];
            int node_community = this->node_to_community[node];
            double weighted_out_degree  = (this->g)->weighted_out_degree(node);
            double weighted_in_degree   = (this->g)->weighted_in_degree(node);
            double self_loops           = (this->g)->count_selfloops(node);

            // Computating all neighboring communities of current node (the number of such communities is stored in neighboring_communities)
            list_neighboring_communities(node, *this, neighbor_weight, positions_neighboring_communities, neighboring_communities);

            // Gain from removing node from its current community
            //start_mod = this->modularity();
            double removal = gain_from_removal(*this, node, node_community, neighbor_weight[node_community], weighted_out_degree, weighted_in_degree);
            remove(*this, node, node_community, neighbor_weight[node_community]+self_loops, weighted_out_degree, weighted_in_degree);

            // Default choice for future insertion is the former community
            int best_community = node_community;
            double best_nbarcs = neighbor_weight[node_community];
            double best_increase = 0.;

            // Computing modularity gain for all neighboring communities
            for (unsigned int i = 0; i < neighboring_communities; ++i) {
                // (Gain from) inserting note to neighboring community
                double insertion = gain_from_insertion(*this, node, positions_neighboring_communities[i], neighbor_weight[positions_neighboring_communities[i]], weighted_out_degree, weighted_in_degree);
                double increase = insertion+removal;
                if(increase > best_increase) {
                    best_community = positions_neighboring_communities[i];
                    best_nbarcs = neighbor_weight[best_community];
                    best_increase = increase;
                }
            }
            // Inserting node in the nearest community
            insert(*this, node, best_community, best_nbarcs+self_loops, weighted_out_degree, weighted_in_degree);

            // If a move was made then we do one more step
            if (best_community != node_community) {
                total_increase+=best_increase;
                improvement = true;
                ++nb_moves;
            }
        }
        
        // Computing the difference between the two modularities
        delta = (current_modularity+total_increase) - current_modularity;
        current_modularity = delta + current_modularity;

    } while (nb_moves > 0 && delta > precision);

    modularity = current_modularity;
    return improvement;
}

void Community::run(bool verbose, const int& display_level, const string& filename_part) {
    int level = 0;
    double mod = this->modularity();

    bool improvement = true;
    if (filename_part != "")
        this->init_partition(filename_part);
    do {
        const Graph *community_graph = this->get_graph();
        if (verbose) {
            cerr << "level " << level << ":\n";
            cerr << "  network size: " <<
                community_graph->get_nodes() << " nodes, " <<
                community_graph->get_arcs() << " arcs, " <<
                community_graph->get_total_weight() << " weight." << endl;
        }

        // Directed Louvain: main procedure
        double new_mod = 0;
        improvement = this->one_level(new_mod);
        if (++level == display_level)
            community_graph->display();
        if (display_level == -1)
            this->display_partition();
        // Updating the graph to computer hierarchical structurer
        this->partition_to_graph();
        if (verbose)
            cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

        mod = new_mod;
        // Doing at least one more computation if partition is provided
        if (filename_part != "" && level == 1) 
            improvement = true;
    } while (improvement);
}

// Friend and static functions are defered to a different file for readability 
#include "community_friend_static.cpp"
