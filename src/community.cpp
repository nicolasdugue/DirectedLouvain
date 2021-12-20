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

#include "../include/community.h"

using namespace std;

Community::Community(string in_filename, int type, int nbp, double minm) {
    cerr << "Reading graph" << endl;
    g = new Graph(in_filename, type);
    cerr << "Graph read" << endl;
    size = ( * g).nb_nodes;
    neigh_weight.resize(size, -1);
    neigh_pos.resize(size);
    neigh_last = 0;

    n2c.resize(size); 
    in .resize(size);
    tot.resize(size);
    tot_out.resize(size);
    tot_in.resize(size);

    for (int i = 0; i < size; i++) {
        // i belongs to its own community
        n2c[i] = i;
        // computing weighted directed degree of i
        tot_out[i] = ( * g).out_weighted_degree(i);
        tot_in[i] = ( * g).in_weighted_degree(i);
        tot[i] = tot_out[i] + tot_in[i];
        /* the total number of edges inside the community corresponds to 
         * the number of self-loops after contraction
         */
        in [i] = ( * g).nb_selfloops(i);
    }

    nb_pass = nbp;
    min_modularity = minm;
}

Community::Community(string bin_filename, string filename_w, int type, int nbp, double minm) {
    cerr << "Reading graph" << endl;
    g = new Graph(bin_filename, filename_w, type);
    cerr << "Graph read" << endl;
    size = ( * g).nb_nodes;
    neigh_weight.resize(size, -1);
    neigh_pos.resize(size);
    neigh_last = 0;

    n2c.resize(size); 
    in .resize(size);
    tot.resize(size);
    tot_out.resize(size);
    tot_in.resize(size);

    for (int i = 0; i < size; i++) {
        // i belongs to its own community
        n2c[i] = i;
        // computing weighted directed degree of i
        tot_out[i] = ( * g).out_weighted_degree(i);
        tot_in[i] = ( * g).in_weighted_degree(i);
        tot[i] = tot_out[i] + tot_in[i];
        /* the total number of edges inside the community corresponds to 
         * the number of self-loops after contraction
         */
        in [i] = ( * g).nb_selfloops(i);
    }

    nb_pass = nbp;
    min_modularity = minm;
}

Community::Community(Graph * gc, int nbp, double minm) {
    g = new Graph(*gc);
    size = ( * g).nb_nodes;

    neigh_weight.resize(size, -1);
    neigh_pos.resize(size);
    neigh_last = 0;

    n2c.resize(size); in .resize(size);
    tot_out.resize(size);
    tot_in.resize(size);
    tot.resize(size);

    for (int i = 0; i < size; i++) {
        n2c[i] = i; in [i] = ( * g).nb_selfloops(i);
        tot_out[i] = ( * g).out_weighted_degree(i);
        tot_in[i] = ( * g).in_weighted_degree(i);
        tot[i] = tot_out[i] + tot_in[i];
    }

    nb_pass = nbp;
    min_modularity = minm;
}

void
Community::init_partition(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in );

    // read partition
    while (!finput.eof()) {
        float node, comm;
        finput >> node >> comm;

        if (finput) {
            int old_comm = n2c[node];
            neigh_comm(node);

            remove(node, old_comm, neigh_weight[old_comm]);

            float i = 0;
            for (i = 0; i < neigh_last; i++) {
                float best_comm = neigh_pos[i];
                double best_nblinks = neigh_weight[neigh_pos[i]];
                if (best_comm == comm) {
                    insert(node, best_comm, best_nblinks);
                    break;
                }
            }
            if (i == neigh_last)
                insert(node, comm, 0);
        }
    }
    finput.close();
}

void
Community::display() {
    for (int i = 0; i < size; i++)
        cerr << " " << i << "/" << n2c[i] << "/" << in [i] << "/" << tot[i];
    cerr << endl;
}

double
Community::modularity() {
    double q = 0.;
    double m = (double)( * g).total_weight;
    for (int i = 0; i < size; i++) {
        if (tot_in[i] > 0 || tot_out[i] > 0) {
            double tot_out_var, tot_in_var;
            tot_out_var = (double) tot_out[i] / m;
            tot_in_var = (double) tot_in[i] / m;
            q += (double) in [i] / m - (tot_out_var * tot_in_var);
        }
    }

    return q;
}

void
Community::neigh_comm(float node) {
    for (float i = 0; i < neigh_last; i++)
        neigh_weight[neigh_pos[i]] = -1;
    // at this stage no neighboring community has to be visited
    neigh_last = 0;

    pair < size_t, size_t > p = ( * g).neighbors(node);

    float deg = ( * g).nb_neighbors_out(node);

    // the first neighboring community of each node is its own
    neigh_pos[0] = n2c[node];
    neigh_weight[neigh_pos[0]] = 0;
    neigh_last = 1;

    for (float i = 0; i < deg; i++) {
        // fetching neighbors of i, their community and the corresponding degrees
        float neigh = g->links[p.first + i];
        float neigh_comm = n2c[neigh];
        double neigh_w = (( * g).weights.size() == 0) ? 1. : g->weights[p.second + i];

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
    pair < size_t, size_t > p_in = ( * g).in_neighbors(node);

    float deg_in = ( * g).nb_neighbors_in(node);

    for (float i = 0; i < deg_in; i++) {

        float neigh_in = g->links_in[p_in.first + i];
        float neigh_comm_in = n2c[neigh_in];
        double neigh_w_in = (( * g).weights_in.size() == 0) ? 1. : g->weights_in[p_in.second + i];

        if (neigh_in != node) {
            if (neigh_weight[neigh_comm_in] == -1) {
                neigh_weight[neigh_comm_in] = 0.;
                neigh_pos[neigh_last++] = neigh_comm_in;
            }
            neigh_weight[neigh_comm_in] += neigh_w_in;
        }
    }
}

void
Community::partition2graph() {
    vector < int > renumber(size, -1);
    for (int node = 0; node < size; node++) {
        renumber[n2c[node]]++;
    }

    int final = 0;
    for (int i = 0; i < size; i++)
        if (renumber[i] != -1)
            renumber[i] = final++;

    for (int i = 0; i < size; i++) {
        pair < size_t, size_t > p = ( * g).neighbors(i);

        int deg = ( * g).nb_neighbors_out(i);
        for (int j = 0; j < deg; j++) {
            int neigh = g->links[p.first + j];
            cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
        }
    }
}

void
Community::display_partition() {
    vector < int > renumber(size, -1);
    for (int node = 0; node < size; node++) {
        renumber[n2c[node]]++;
    }

    int final = 0;
    for (int i = 0; i < size; i++)
        if (renumber[i] != -1)
            renumber[i] = final++;

    for (int i = 0; i < size; i++)
        cout << i << " " << renumber[n2c[i]] << endl;
}

Graph *
    Community::partition2graph_binary() {
        // renumber communities
        vector < int > renumber(size, -1);
        for (int node = 0; node < size; node++) {
            renumber[n2c[node]]++;
        }

        // give a number to every community
        int final = 0;
        for (int i = 0; i < size; i++)
            if (renumber[i] != -1)
                renumber[i] = final++;

        // compute communities
        vector < vector < int > > comm_nodes(final);
        for (int node = 0; node < size; node++) {
            comm_nodes[renumber[n2c[node]]].push_back(node);
        }

        // compute weighted graph
        Graph *g2 = new Graph();
        (*g2).nb_nodes = comm_nodes.size();

        (*g2).degrees_out.resize(comm_nodes.size());
        (*g2).degrees_in.resize(comm_nodes.size());

        double neigh_weight;

        int comm_deg = comm_nodes.size();
        for (int comm = 0; comm < comm_deg; comm++) {
            map < int, double > m_out, m_in;
            map < int, double > ::iterator it_out, it_in;

            int comm_size = comm_nodes[comm].size();
            for (int node = 0; node < comm_size; node++) {
                // we first deal with out-neighbors communities
                pair < size_t, size_t > p = ( * g).neighbors(comm_nodes[comm][node]);
                int deg = ( * g).nb_neighbors_out(comm_nodes[comm][node]);
                for (int i = 0; i < deg; i++) {
                    int neigh = g->links[p.first + i];
                    int neigh_comm = renumber[n2c[neigh]];
                    neigh_weight = (( * g).weights.size() == 0) ? 1. : g->weights[p.second + i];

                    it_out = m_out.find(neigh_comm);
                    if (it_out == m_out.end())
                        m_out.insert(make_pair(neigh_comm, neigh_weight));
                    else
                        it_out -> second += neigh_weight;
                }

                // same thing for in-neighbors communities
                pair < size_t, size_t > p_in = ( * g).in_neighbors(comm_nodes[comm][node]);
                deg = ( * g).nb_neighbors_in(comm_nodes[comm][node]);
                for (int i = 0; i < deg; i++) {
                    int neigh = g->links_in[p_in.first + i];
                    int neigh_comm = renumber[n2c[neigh]];
                    neigh_weight = (( * g).weights_in.size() == 0) ? 1.f : g->weights_in[p_in.second + i];

                    it_in = m_in.find(neigh_comm);
                    if (it_in == m_in.end())
                        m_in.insert(make_pair(neigh_comm, neigh_weight));
                    else
                        it_in -> second += neigh_weight;
                }
            }

            ( * g2).degrees_out[comm] = (comm == 0) ? m_out.size() : ( * g2).degrees_out[comm - 1] + m_out.size();
            ( * g2).nb_links_out += m_out.size();

            for (it_out = m_out.begin(); it_out != m_out.end(); it_out++) {
                ( * g2).total_weight += it_out -> second;
                ( * g2).links.push_back(it_out -> first);
                ( * g2).weights.push_back(it_out -> second);
            }

            ( * g2).degrees_in[comm] = (comm == 0) ? m_in.size() : ( * g2).degrees_in[comm - 1] + m_in.size();
            ( * g2).nb_links_in += m_in.size();

            for (it_in = m_in.begin(); it_in != m_in.end(); it_in++) {
                ( * g2).links_in.push_back(it_in -> first);
                ( * g2).weights_in.push_back(it_in -> second);
            }
        }

        return g2;
    }

bool
Community::one_level() {
    bool improvement = false;
    int nb_moves;
    int nb_pass_done = 0;
    double new_mod = modularity();
    double cur_mod = new_mod;

    vector < int > random_order(size);
    for (int i = 0; i < size; i++)
        random_order[i] = i;
    /* FIXME: commented for the sake of reproducibility 
     * TODO: random_shuffle
     */
    /*for (int i = 0; i < size - 1; i++) {
        int rand_pos = rand() % (size - i) + i;
        int tmp = random_order[i];
        random_order[i] = random_order[rand_pos];
        random_order[rand_pos] = tmp;
    }*/

    // repeat while
    //   there is an improvement of modularity
    //   or there is an improvement of modularity greater than a given epsilon
    //   or a predefined number of pass have been done
    do {
        cur_mod = new_mod;
        nb_moves = 0;
        nb_pass_done++;

        // for each node: remove the node from its community and insert it in the best community
        for (int node_tmp = 0; node_tmp < size; node_tmp++) {
            int node = random_order[node_tmp];
            int node_comm = n2c[node];
            double w_degree_out = ( * g).out_weighted_degree(node);
            double w_degree_in = ( * g).in_weighted_degree(node);

            // computation of all neighboring communities of current node
            neigh_comm(node);
            // remove node from its current community
            remove(node, node_comm, neigh_weight[node_comm]);

            // compute the nearest community for node
            // default choice for future insertion is the former community
            int best_comm = node_comm;
            double best_nblinks = 0.;
            double best_increase = 0.;
            for (float i = 0; i < neigh_last; i++) {
                double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree_out, w_degree_in);
                if (increase > best_increase) {
                    best_comm = neigh_pos[i];
                    best_nblinks = neigh_weight[neigh_pos[i]];
                    best_increase = increase;
                }
            }

            // insert node in the nearest community
            insert(node, best_comm, best_nblinks);

            if (best_comm != node_comm)
                nb_moves++;

            //if(node_tmp % 5000000 ==0) cerr << "one less to go" << endl;
        }

        new_mod = modularity();

        if (nb_moves > 0)
            improvement = true;

    } while (nb_moves > 0 && new_mod - cur_mod > min_modularity);

    return improvement;
}
