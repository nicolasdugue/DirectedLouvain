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

double modularity_gain(const Community &c, unsigned int node, unsigned int comm, double dnodecomm) {
    assert(node<c.size);
    double weighted_out_degree  = (c.g)->weighted_out_degree(node);
    double weighted_in_degree   = (c.g)->weighted_in_degree(node);
    double totc_out             = c.communities_arcs[comm].tot_out;
    double totc_in              = c.communities_arcs[comm].tot_in;
    double m                    = (c.g)->get_total_weight();

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

void list_neighboring_communities(unsigned int node, const Community &c, vector<double> &neighbor_weight, vector<unsigned int> &neigh_pos, unsigned int &neighboring_communities) {
    for(unsigned int i = 0; i < neighboring_communities; ++i)
        /* FIXME: clear then resize? */
        neighbor_weight[neigh_pos[i]] = -1.f;

    size_t p = (c.g)->out_neighbors(node);
    unsigned int deg = (c.g)->out_degree(node);

    // the first neighboring community of each node is its own
    neigh_pos[0] = c.get_community(node);
    neighbor_weight[neigh_pos[0]] = 0;
    neighboring_communities = 1;

    for (unsigned int i = 0; i < deg; ++i) {
        // fetching neighbors of i, their community and the corresponding degrees
        unsigned int neigh = (c.g)->get_out_neighbor(p + i);
        int neigh_comm = c.get_community(neigh);
        double neigh_w = ((c.g)->is_weighted()) ? (c.g)->get_weighted_out_neighbor(p + i) : 1.f;

        if (neigh != node) {
            // if the community is discovered for the first time
            if (neighbor_weight[neigh_comm] == -1.f) {
                neighbor_weight[neigh_comm] = 0.f;
                neigh_pos[neighboring_communities++] = neigh_comm;
            }
            // the degree of i towards this community is updated
            neighbor_weight[neigh_comm] += neigh_w;
        }
    }

    // we proceed similarly on in-neighbors
    /* FIXME: don't we count twice the same thing ?! 
     * or do we really need to know w_out+w_in ?*
     * FIXME: same question for neigh_pos: encompasses the weights 
     * of out-comm then in-comm?
     */
    size_t p_in = (c.g)->in_neighbors(node);
    unsigned int deg_in = (c.g)->in_degree(node);

    for (unsigned int i = 0; i < deg_in; ++i) {
        unsigned int neigh_in = (c.g)->get_in_neighbor(p_in + i);
        int neigh_comm_in = c.get_community(neigh_in);
        double neigh_w_in = ((c.g)->is_weighted()) ? (c.g)->get_weighted_in_neighbor(p_in + i) : 1.f;

        if (neigh_in != node) {
            if (neighbor_weight[neigh_comm_in] == -1) {
                neighbor_weight[neigh_comm_in] = 0.;
                neigh_pos[neighboring_communities++] = neigh_comm_in;
            }
            neighbor_weight[neigh_comm_in] += neigh_w_in;
        }
    }
}

