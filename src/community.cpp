#include "../include/community.hpp"
#include <climits>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

// Static function renumbering communities from 0 to k-1 (returns k)
static unsigned int renumber_communities(const Community &c, vector< int > &renumber);
static void update_levels(const Community &c, vector< vector<int> > &levels, int level);

Community::Community(const string& in_filename, const double precision, const double gamma, bool reproducibility, bool renumbering, bool randomized) {
    this->g                 = new Graph(in_filename, reproducibility, renumbering);
    this->precision         = precision;
    this->gamma             = gamma;
    this->randomized        = randomized;
    init_attributes();
}

void Community::init_attributes() {
    this->size              = g->nodes;
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
    delete this->community_graph;
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
        cout << " " << (this->g)->correspondance[i] << "/" << node_to_community[i] << "/" 
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
        cout << (this->community_graph)->correspondance[i] << " " << renumber[this->node_to_community[i]] << endl;
}

map<unsigned int, unsigned int> Community::get_level(int level){
    assert(level >= 0 && level < (int)this->levels.size());
    vector < int > n2c(this->g->nodes);
    map<unsigned int, unsigned int> lvl;

    for (unsigned int i = 0; i < this->g->nodes; i++)
        n2c[i] = i;

    for (int l = 0; l < level; l++)
        for (unsigned int node = 0; node < this->g->nodes; node++)
            n2c[node] = this->levels[l][n2c[node]];

    for (unsigned int node = 0; node < this->g->nodes; node++)
        lvl[(this->g)->correspondance[node]] = n2c[node];
    return lvl;
}

map<unsigned int, unsigned int> Community::get_last_level(){
    return get_level(levels.size()-1);
}

void Community::print_level(int level) {
    assert(level >= 0 && level < (int)this->levels.size());
    vector < int > n2c(this->g->nodes);

    for (unsigned int i = 0; i < this->g->nodes; i++)
        n2c[i] = i;

    for (int l = 0; l < level; l++)
        for (unsigned int node = 0; node < this->g->nodes; node++)
            n2c[node] = this->levels[l][n2c[node]];

    for (unsigned int node = 0; node < this->g->nodes; node++) 
        cout << (this->g)->correspondance[node] << " " << n2c[node] << endl;
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
            size_t p = (this->community_graph)->out_neighbors(comm_nodes[comm][node]);
            unsigned int deg = (this->community_graph)->out_degree(comm_nodes[comm][node]);
            // Looking for communities of every out-neighbor of node and then storing/updating weighted out-degrees
            for (unsigned int i = 0; i < deg; ++i) {
                out_neighbor = (this->community_graph)->outcoming_arcs[p + i];
                out_neighboring_community = renumber[this->node_to_community[out_neighbor]];
                out_weight = ((this->community_graph)->weighted) ? (this->community_graph)->outcoming_weights[p + i] : 1.f;

                auto it_out = m_out.find(out_neighboring_community);
                if (it_out == m_out.end())
                    m_out.insert(make_pair(out_neighboring_community, out_weight));
                else
                    it_out -> second += out_weight;
            }

            // In-neighbors
            size_t p_in = (this->community_graph)->in_neighbors(comm_nodes[comm][node]);
            deg = (this->community_graph)->in_degree(comm_nodes[comm][node]);
            // Looking for communities of every in-neighbor of node and then storing/updating weighted in-degrees
            for (unsigned int i = 0; i < deg; ++i) {
                in_neighbor = (this->community_graph)->incoming_arcs[p_in + i];
                in_neighboring_community = renumber[this->node_to_community[in_neighbor]];
                in_weight = ((this->community_graph)->weighted) ? (this->community_graph)->incoming_weights[p_in + i] : 1.f;

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
    delete this->community_graph;
    this->community_graph = g2;

    // Updating other attributes according to constructed graph g
    this->size           = this->community_graph->nodes;
    this->node_to_community.resize(this->size); 

    for (unsigned int i = 0; i < this->size; ++i) {
        this->node_to_community[i]                      = i; 
        this->communities_arcs[i].total_arcs_inside     = this->community_graph->count_selfloops(i);
        this->communities_arcs[i].total_outcoming_arcs  = this->community_graph->weighted_out_degree(i);
        this->communities_arcs[i].total_incoming_arcs   = this->community_graph->weighted_in_degree(i);
    }
}

bool Community::one_level(double &modularity) {
    int nb_moves = 0;
    bool improvement = false;
    double current_modularity = this->modularity();

    // Order in which to proceed nodes of the graph...
    vector < int > random_order(size);
    for (unsigned int i = 0; i < size; ++i)
        random_order[i] = i;

    // ... randomized: (Directed) Louvain's algorithm is not deterministic
    if(randomized) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle(random_order.begin(), random_order.end(), std::default_random_engine(seed));
    }

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
            double weighted_out_degree  = (this->community_graph)->weighted_out_degree(node);
            double weighted_in_degree   = (this->community_graph)->weighted_in_degree(node);
            double self_loops           = (this->community_graph)->count_selfloops(node);

            // Computating all neighboring communities of current node (the number of such communities is stored in neighboring_communities)
            list_neighboring_communities(node, *this, neighbor_weight, positions_neighboring_communities, neighboring_communities);

            // Gain from removing node from its current community
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
                //double increase = insertion;
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
        current_modularity += total_increase;

    } while (nb_moves > 0 && total_increase > precision);

    modularity = this->modularity();
    return improvement;
}

int Community::run(bool verbose, const int& display_level, const string& filename_part, bool egc, unsigned int nb_runs=1) {
    int level = 0;

    if(verbose) {
        cerr << "Directed Louvain runs on the following parameters: " << endl; 
        cerr << "graph with: " << 
                this->g->get_nodes() << " nodes, " <<
                this->g->get_arcs() << " arcs, " <<
                this->g->get_total_weight() << " weight." << endl;
        cerr << "precision is set to: " << this->precision << endl;
        cerr << "resolution parameter is set to: " << this->gamma << endl;
    }

    bool improvement = true;
    if (filename_part != "")
        this->init_partition(filename_part);

    if(egc) {
        cerr << "computing Ensemble Graph" << endl; 
        this->community_graph = this->undirected_egc_graph(nb_runs);
    }

    else
        this->community_graph   = new Graph(*(this->g));

    this->init_attributes();
    double mod = this->modularity();
    cerr << mod << endl;
    do {
        if (verbose) {
            cerr << "level " << level << ":\n";
            cerr << "  network size: " <<
                this->community_graph->get_nodes() << " nodes, " <<
                this->community_graph->get_arcs() << " arcs, " <<
                this->community_graph->get_total_weight() << " weight." << endl;
        }

        // Directed Louvain: main procedure
        double new_mod = 0;
        improvement = this->one_level(new_mod);
        // Maintaining levels
        levels.resize(++level);
        update_levels(*this, levels, level-1);
        if (level == display_level || display_level == -1)
            this->display_partition();
        // Updating the graph to computer hierarchical structure
        this->partition_to_graph();
        if (verbose)
            cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

        mod = new_mod;
        // Doing at least one more computation if partition is provided
        if (filename_part != "" && level == 1) 
            improvement = true;
    } while (improvement);
    if (display_level == -2)
        print_level(levels.size()-1);
    return level;
}

Graph* Community::undirected_egc_graph(unsigned int nb_runs) {
    Graph *tmp = new Graph(*(this->g));

    //Ensemble-Clustering-for-Graphs
    cerr << "computing votes ensemble..." << endl;
    map<pair<unsigned int, unsigned int>, unsigned int> weights;
    for (unsigned int i = 0; i < nb_runs; ++i) {
        // Directed Louvain: main procedure
        this->community_graph = new Graph(*(this->g));
        this->init_attributes();

        double new_mod = 0;
        this->one_level(new_mod);
        cerr << "modularity at step " << i << ":" << this->modularity() << endl;
        for (unsigned int origin = 0; origin < this->g->nodes; origin++) {
            for (unsigned int j = 0; j < this->g->out_degree(origin); ++j) {
                size_t p = this->g->out_neighbors(origin);
                int destination = this->g->outcoming_arcs[p+j];
                if(this->node_to_community[origin] == this->node_to_community[destination]){
                    pair<unsigned int, unsigned int> to_add(origin, destination);
                    weights[to_add]++;
                }
            }
        }
        delete this->community_graph;
    }
    cerr << "done" << endl;

    cerr << "computing core indices" << endl;

    vector<unsigned int> bin;
    vector<unsigned int> vert(this->g->nodes);
    vector<unsigned int> pos(this->g->nodes);
    
    /* Computing max degree */
    vector<unsigned int> cores(this->g->nodes);

    /* We do not consider weights for this algorithm, is this correct? */
    for(size_t i = 0; i < this->g->nodes; ++i)
        cores[i] = this->g->out_degree(i)+this->g->in_degree(i);

    unsigned int max_degree = *(max_element(cores.begin(), cores.end()));

    /* Computing degree histogram */
    bin.resize(max_degree+1,0);
    for(size_t i = 0; i < this->g->nodes; ++i)
        bin[cores[i]] += 1;

    /* Start pointers */
    unsigned int j = 0;
    for(unsigned int i = 0; i < max_degree+1; ++i) {
        unsigned int k = bin[i];
        bin[i] = j;
        j += k;
    }

    /* Sorting in vert (and corrupting bin) */
    for(size_t i = 0; i < this->g->nodes; ++i) {
        pos[i] = bin[cores[i]];
        vert[pos[i]] = i;
        bin[cores[i]] += 1;
    }

    /* Correcting bin */
    for(size_t i = max_degree; i > 0; --i)
        bin[i] = bin[i-1];

    /* Main algorithm */
    for(size_t i = 0; i < this->g->nodes; ++i) {
        unsigned int v = vert[i];
        size_t p = this->g->out_neighbors(v);
        for(size_t j = 0; j < this->g->out_degree(v); ++j) {
            unsigned int u = this->g->outcoming_arcs[p+j];
            if(cores[u] > cores[v]) {
                unsigned int du = cores[u];
                unsigned int pu = pos[u];
                unsigned int pw = bin[du];
                unsigned int w = vert[pw];
                if(u != w) {
                    pos[u] = pw; 
                    vert[pu] = w;
                    pw = bin[du];
                    w = vert[pw];
                }
                bin[du] += 1;
                cores[u] -= 1;
            }
        }
    } 

    cerr << "computing ensemble graph" << endl;
    ofstream foutput;
    foutput.open("graph/EGC.txt", fstream::out | fstream::binary);
    // Computing weighted graph from previous steps
    double min_weight = .05;
    tmp->outcoming_weights.assign(tmp->outcoming_weights.size(), min_weight);
    tmp->incoming_weights.assign(tmp->incoming_weights.size(), min_weight);

    /* Out-cores */
    for (unsigned int node = 0; node < tmp->nodes; ++node) {
        size_t p = tmp->out_neighbors(node);
        for (unsigned int i = 0; i < tmp->out_degree(node); ++i) {
            // If the arc is not in the map or has one core value false we assign min_weight to it
            if(weights.count(make_pair(node,tmp->outcoming_arcs[p + i]))>0) 
                if(cores[i]>1 && cores[tmp->outcoming_arcs[p + i]] > 1)  {
                    tmp->outcoming_weights[p + i] = min_weight + (1-min_weight)*(double)weights[make_pair(node,tmp->outcoming_arcs[p + i])]/nb_runs;
                    unsigned int pos_in_neighbor = tmp->in_neighbors(tmp->outcoming_arcs[p+i]);
                    for(size_t j= 0; j < tmp->in_degree(tmp->outcoming_arcs[p+i]); ++j) {
                        if(tmp->incoming_arcs[pos_in_neighbor+j]==node) {
                            tmp->incoming_weights[pos_in_neighbor+j] = min_weight + ((1-min_weight)*((double)weights[make_pair(node,tmp->outcoming_arcs[p + i])]/nb_runs));
                        }
                    }
                }
            foutput << g->correspondance[node] << " " << g->correspondance[tmp->outcoming_arcs[p + i]] << " " << 
            tmp->outcoming_weights[p+i] << endl;
        }
    }

    cerr << 
    accumulate(tmp->outcoming_weights.begin(), tmp->outcoming_weights.end(), decltype(tmp->outcoming_weights)::value_type(0)) << " " <<  
    accumulate(tmp->incoming_weights.begin(), tmp->incoming_weights.end(), decltype(tmp->incoming_weights)::value_type(0)) << 
    endl;
    
    tmp->total_weight =  accumulate(tmp->outcoming_weights.begin(), tmp->outcoming_weights.end(), decltype(tmp->outcoming_weights)::value_type(0)); 

    foutput.close();
    
    return tmp;
}

/* Inspired by https://github.com/igraph/igraph/blob/7632007bdfd837bfc68d5087ce6f34f1e9139385/src/centrality/coreness.c#L31 
 * Extracted from https://arxiv.org/pdf/cs/0310049.pdf Algorithm 1
 */
Graph* Community::linear_egc_graph(unsigned int nb_runs) {
    Graph *tmp = new Graph(*(this->g));

    //Ensemble-Clustering-for-Graphs
    cerr << "computing votes ensemble..." << endl;
    map<pair<unsigned int, unsigned int>, unsigned int> weights;
    for (unsigned int i = 0; i < nb_runs; ++i) {
        // Directed Louvain: main procedure
        this->community_graph = new Graph(*(this->g));
        this->init_attributes();
        double new_mod = 0;
        this->one_level(new_mod);
        for (unsigned int origin = 0; origin < this->g->nodes; origin++) {
            for (unsigned int j = 0; j < this->g->out_degree(origin); ++j) {
                size_t p = this->g->out_neighbors(origin);
                int destination = this->g->outcoming_arcs[p+j];
                if(this->node_to_community[origin] == this->node_to_community[destination]){
                    pair<unsigned int, unsigned int> to_add(origin, destination);
                    weights[to_add]++;
                }
            }
        }
        delete this->community_graph;
    }
    cerr << "done" << endl;

    cerr << "computing core indices" << endl;

    vector<unsigned int> bin;
    vector<unsigned int> vert(this->g->nodes);
    vector<unsigned int> pos(this->g->nodes);
    
    /* Computing max degree */
    vector<unsigned int> out_cores(this->g->nodes);

    /* We do not consider weights for this algorithm, is this correct? */
    for(size_t i = 0; i < this->g->nodes; ++i)
        out_cores[i] = this->g->out_degree(i);

    unsigned int out_max_degree = *(max_element(out_cores.begin(), out_cores.end()));

    /* Computing degree histogram */
    bin.resize(out_max_degree+1,0);
    for(size_t i = 0; i < this->g->nodes; ++i)
        bin[out_cores[i]] += 1;

    /* Start pointers */
    unsigned int j = 0;
    for(unsigned int i = 0; i < out_max_degree+1; ++i) {
        unsigned int k = bin[i];
        bin[i] = j;
        j += k;
    }

    /* Sorting in vert (and corrupting bin) */
    for(size_t i = 0; i < this->g->nodes; ++i) {
        pos[i] = bin[out_cores[i]];
        vert[pos[i]] = i;
        bin[out_cores[i]] += 1;
    }

    /* Correcting bin */
    for(size_t i = out_max_degree; i > 0; --i)
        bin[i] = bin[i-1];

    /* Main algorithm */
    for(size_t i = 0; i < this->g->nodes; ++i) {
        unsigned int v = vert[i];
        size_t p = this->g->out_neighbors(v);
        for(size_t j = 0; j < this->g->out_degree(v); ++j) {
            unsigned int u = this->g->outcoming_arcs[p+j];
            if(out_cores[u] > out_cores[v]) {
                unsigned int du = out_cores[u];
                unsigned int pu = pos[u];
                unsigned int pw = bin[du];
                unsigned int w = vert[pw];
                if(u != w) {
                    pos[u] = pw; 
                    vert[pu] = w;
                    pw = bin[du];
                    w = vert[pw];
                }
                bin[du] += 1;
                out_cores[u] -= 1;
            }
        }
    } 

    /* Computing max degree */
    vector<unsigned int> in_cores(this->g->nodes);
    vert.clear();
    pos.clear();
    bin.clear();

    /* We do not consider weights for this algorithm, is this correct? */
    for(size_t i = 0; i < this->g->nodes; ++i)
        in_cores[i] = this->g->in_degree(i);

    unsigned int in_max_degree = *(max_element(in_cores.begin(), in_cores.end()));

    /* Computing degree histogram */
    bin.resize(in_max_degree+1,0);
    for(size_t i = 0; i < this->g->nodes; ++i)
        bin[in_cores[i]] += 1;

    /* Start pointers */
    j = 0;
    for(unsigned int i = 0; i < in_max_degree+1; ++i) {
        unsigned int k = bin[i];
        bin[i] = j;
        j += k;
    }

    /* Sorting in vert (and corrupting bin) */
    for(size_t i = 0; i < this->g->nodes; ++i) {
        pos[i] = bin[in_cores[i]];
        vert[pos[i]] = i;
        bin[in_cores[i]] += 1;
    }

    /* Correcting bin */
    for(size_t i = in_max_degree; i > 0; --i)
        bin[i] = bin[i-1];

    /* Main algorithm */
    for(size_t i = 0; i < this->g->nodes; ++i) {
        unsigned int v = vert[i];
        size_t p = this->g->in_neighbors(v);
        for(size_t j = 0; j < this->g->in_degree(v); ++j) {
            unsigned int u = this->g->incoming_arcs[p+j];
            if(in_cores[u] > in_cores[v]) {
                unsigned int du = in_cores[u];
                unsigned int pu = pos[u];
                unsigned int pw = bin[du];
                unsigned int w = vert[pw];
                if(u != w) {
                    pos[u] = pw; 
                    vert[pu] = w;
                    pw = bin[du];
                    w = vert[pw];
                }
                bin[du] += 1;
                in_cores[u] -= 1;
            }
        }
    } 
    cerr << "done" << endl;

    cerr << "computing ensemble graph" << endl;
    ofstream foutput;
    foutput.open("graph/EGC.txt", fstream::out | fstream::binary);
    // Computing weighted graph from previous steps
    double min_weight = .05;
    tmp->outcoming_weights.assign(tmp->outcoming_weights.size(), min_weight);
    tmp->incoming_weights.assign(tmp->incoming_weights.size(), min_weight);

    /* Out-cores */
    for (unsigned int node = 0; node < tmp->nodes; ++node) {
        size_t p = tmp->out_neighbors(node);
        for (unsigned int i = 0; i < tmp->out_degree(node); ++i) {
            // If the arc is not in the map or has one core value false we assign min_weight to it
            if(weights.count(make_pair(node,tmp->outcoming_arcs[p + i]))>0) 
                if((out_cores[i]>1 && out_cores[tmp->outcoming_arcs[p + i]] > 1) || (in_cores[i]>1 && in_cores[tmp->outcoming_arcs[p + i]] > 1))  {
                    tmp->outcoming_weights[p + i] = min_weight + (1-min_weight)*(double)weights[make_pair(node,tmp->outcoming_arcs[p + i])]/nb_runs;
                    unsigned int pos_in_neighbor = tmp->in_neighbors(tmp->outcoming_arcs[p+i]);
                    for(size_t j= 0; j < tmp->in_degree(tmp->outcoming_arcs[p+i]); ++j) {
                        if(tmp->incoming_arcs[pos_in_neighbor+j]==node) {
                            tmp->incoming_weights[pos_in_neighbor+j] = min_weight + ((1-min_weight)*((double)weights[make_pair(node,tmp->outcoming_arcs[p + i])]/nb_runs));
                        }
                    }
                }
            foutput << g->correspondance[node] << " " << g->correspondance[tmp->outcoming_arcs[p + i]] << " " << 
            tmp->outcoming_weights[p+i] << endl;
        }
    }

    cerr << 
    accumulate(tmp->outcoming_weights.begin(), tmp->outcoming_weights.end(), decltype(tmp->outcoming_weights)::value_type(0)) << " " <<  
    accumulate(tmp->incoming_weights.begin(), tmp->incoming_weights.end(), decltype(tmp->incoming_weights)::value_type(0)) << 
    endl;
    
    tmp->total_weight =  accumulate(tmp->outcoming_weights.begin(), tmp->outcoming_weights.end(), decltype(tmp->outcoming_weights)::value_type(0)); 

    foutput.close();
    
    return tmp;
}

Graph* Community::egc_graph(unsigned int nb_runs) {
    Graph *tmp = new Graph(*(this->g));

    //Ensemble-Clustering-for-Graphs
    cerr << "computing votes ensemble..." << endl;
    map<pair<unsigned int, unsigned int>, unsigned int> weights;
    for (unsigned int i = 0; i < nb_runs; ++i) {
        // Directed Louvain: main procedure
        this->community_graph = new Graph(*(this->g));
        this->init_attributes();
        double new_mod = 0;
        this->one_level(new_mod);
        for (unsigned int origin = 0; origin < this->g->nodes; origin++) {
            for (unsigned int j = 0; j < this->g->out_degree(origin); ++j) {
                size_t p = this->g->out_neighbors(origin);
                int destination = this->g->outcoming_arcs[p+j];
                if(this->node_to_community[origin] == this->node_to_community[destination]){
                    pair<unsigned int, unsigned int> to_add(origin, destination);
                    weights[to_add]++;
                }
            }
        }
        delete this->community_graph;
    }
    cerr << "done" << endl;
    
    cerr << "computing 2-core" << endl;
    // 2-core
    vector<bool> core(size, true);
    Graph* core_graph = new Graph(*(this->g));
    bool marked = false;
    do {
        marked = false;
        for (unsigned int i = 0; i < core.size(); ++i) {
            if (core[i]){
                unsigned int outdegree = core_graph->out_degree(i);
                unsigned int indegree = core_graph->in_degree(i);
                unsigned int degree = outdegree + indegree;
                if (degree < 2){
                    core[i] = false;
                    marked = true;
                    if(outdegree==1) {
                        unsigned int pos_out_neighbor = core_graph->out_neighbors(i);
                        unsigned int in_neighbor = core_graph->outcoming_arcs[pos_out_neighbor];
                        unsigned int pos_in_neighbor = core_graph->in_neighbors(in_neighbor);
                        vector<unsigned int>::iterator it;
                        it = find(core_graph->incoming_arcs.begin()+pos_in_neighbor, core_graph->incoming_arcs.begin()+pos_in_neighbor + core_graph->in_degree(in_neighbor), i);
                        if(it!=core_graph->incoming_arcs.end()) {
                            core_graph->incoming_arcs.erase(it);
                            core_graph->arcs--;
                        }
                        for (unsigned int j = i; j < core_graph->outdegrees.size(); ++j)
                            core_graph->outdegrees[j]--;
                        core_graph->outcoming_arcs.erase(core_graph->outcoming_arcs.begin()+pos_out_neighbor);
                        for (unsigned int j = in_neighbor; j < core_graph->indegrees.size(); ++j)
                            core_graph->indegrees[j]--;
                    }
                    if(indegree == 1){
                        unsigned int pos_in_neighbor = core_graph->in_neighbors(i);
                        unsigned int out_neighbor = core_graph->incoming_arcs[pos_in_neighbor];
                        unsigned int pos_out_neighbor = core_graph->out_neighbors(out_neighbor);
                        vector<unsigned int>::iterator it;
                        it = find(core_graph->outcoming_arcs.begin()+pos_out_neighbor, core_graph->outcoming_arcs.begin()+pos_out_neighbor + core_graph->out_degree(out_neighbor), i);
                        if(it!=core_graph->outcoming_arcs.end()) {
                            core_graph->outcoming_arcs.erase(it);
                            core_graph->arcs--;
                        }
                        for (unsigned int j = i; j < core_graph->indegrees.size(); ++j)
                            core_graph->indegrees[j]--;
                        core_graph->incoming_arcs.erase(core_graph->incoming_arcs.begin()+pos_in_neighbor);
                        for (unsigned int j = out_neighbor; j < core_graph->outdegrees.size(); ++j)
                            core_graph->outdegrees[j]--;
                    }
                }
            }
        }
    } while(marked);

    cerr << "done" << endl;
    cerr << "computing ensemble graph" << endl;
    ofstream foutput;
    foutput.open("graph/EGC.txt", fstream::out | fstream::binary);
    // Computing weighted graph from previous steps
    double min_weight = .05;
    tmp->outcoming_weights.assign(tmp->outcoming_weights.size(), min_weight);
    tmp->incoming_weights.assign(tmp->incoming_weights.size(), min_weight);
    for (unsigned int node = 0; node < tmp->nodes; ++node) {
        size_t p = tmp->out_neighbors(node);
        for (unsigned int i = 0; i < tmp->out_degree(node); ++i) {
            // If the arc is not in the map or has one core value false we assign min_weight to it
            if(weights.count(make_pair(node,tmp->outcoming_arcs[p + i]))>0) 
                if(true) {
                //if(core[i] && core[tmp->outcoming_arcs[p + i]]) {
                    tmp->outcoming_weights[p + i] = min_weight + ((1-min_weight)*((double)weights[make_pair(node,tmp->outcoming_arcs[p + i])]/50));
                    unsigned int pos_in_neighbor = tmp->in_neighbors(tmp->outcoming_arcs[p+i]);
                    for(size_t j= 0; j < tmp->in_degree(tmp->outcoming_arcs[p+i]); ++j) {
                        if(tmp->incoming_arcs[pos_in_neighbor+j]==node) {
                            tmp->incoming_weights[pos_in_neighbor+j] = min_weight + ((1-min_weight)*((double)weights[make_pair(node,tmp->outcoming_arcs[p + i])]/50));
                        }
                    }
                }
            foutput << g->correspondance[node] << " " << g->correspondance[tmp->outcoming_arcs[p + i]] << " " << 
            tmp->outcoming_weights[p+i] << endl;
        }
    }

    cerr << 
    accumulate(tmp->outcoming_weights.begin(), tmp->outcoming_weights.end(), decltype(tmp->outcoming_weights)::value_type(0)) << " " <<  
    accumulate(tmp->incoming_weights.begin(), tmp->incoming_weights.end(), decltype(tmp->incoming_weights)::value_type(0)) << 
    endl;
    
    tmp->total_weight =  accumulate(tmp->outcoming_weights.begin(), tmp->outcoming_weights.end(), decltype(tmp->outcoming_weights)::value_type(0)); 

    foutput.close();
    
    return tmp;
}
// Friend and static functions are defered to a different file for readability 
#include "community_friend_static.cpp"
