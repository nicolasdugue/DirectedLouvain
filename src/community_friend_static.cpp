static unsigned int renumber_communities(const Community &c, vector< int > &renumber) {
    size_t size = c.get_size();
    for (unsigned int node = 0; node < size; ++node) 
        ++renumber[c.get_community(node)];
    
    unsigned int f = 0;
    for (unsigned int i = 0; i < size; ++i)
        if (renumber[i] != -1)
            renumber[i] = f++;
    return f;
}

static void update_levels(const Community &c, vector< vector<int> > &levels, int level) {
    vector < int > renumber(c.get_size(), -1);
    renumber_communities(c, renumber);
    for (unsigned int i = 0; i < c.get_size(); ++i)
        levels[level].push_back(renumber[c.get_community(i)]);
}

// Function updating the total number of arcs going from, to and inside a community after the removal of a node
void remove(Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree) {
    assert(node<c.size);
    c.communities_arcs[comm].total_outcoming_arcs   -= weighted_out_degree;
    c.communities_arcs[comm].total_incoming_arcs    -= weighted_in_degree;
    c.communities_arcs[comm].total_arcs_inside      -= dnodecomm;
    c.node_to_community[node]                       = -1;
}

// Function updating the total number of arcs going from, to and inside a community after the insertion of a node
void insert(Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree) {
    assert(node<c.size);
    c.communities_arcs[comm].total_outcoming_arcs   += weighted_out_degree;
    c.communities_arcs[comm].total_incoming_arcs    += weighted_in_degree;
    c.communities_arcs[comm].total_arcs_inside      += dnodecomm;
    c.node_to_community[node]                       = comm;
}

// Function computing the modularity gain from inserting a node within a given community
double gain_from_removal(const Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree) {
    assert(node<c.size);
    double totc_out             = c.communities_arcs[comm].total_outcoming_arcs;
    double totc_in              = c.communities_arcs[comm].total_incoming_arcs;
    double m                    = (c.community_graph)->get_total_weight();
    return ((-dnodecomm / m) + c.gamma * ((weighted_out_degree * totc_in + weighted_in_degree * totc_out)) / (m*m));
}

// Function computing the modularity gain from inserting a node within a given community
double gain_from_insertion(const Community &c, const unsigned int& node, const int& comm, const double& dnodecomm, const double& weighted_out_degree, const double& weighted_in_degree) {
    assert(node<c.size);
    double totc_out             = c.communities_arcs[comm].total_outcoming_arcs + weighted_out_degree;
    double totc_in              = c.communities_arcs[comm].total_incoming_arcs + weighted_in_degree;
    double m                    = (c.community_graph)->get_total_weight();
    return ((dnodecomm / m) - c.gamma * ((weighted_out_degree * totc_in + weighted_in_degree * totc_out)) / (m*m));
}

// Function computing the weights and positions of communities neighboring a given node (including its own) 
void list_neighboring_communities(const unsigned int& node, const Community &c, vector<double> &neighbor_weight, vector<unsigned int> &positions_neighboring_communities, unsigned int &neighboring_communities) {
    // Cleaning the previously computed neighbors
    for(unsigned int i = 0; i < neighboring_communities; ++i)
        neighbor_weight[positions_neighboring_communities[i]] = -1.f;

    size_t p = (c.community_graph)->out_neighbors(node);
    unsigned int deg = (c.community_graph)->out_degree(node);

    // The first neighboring community of each node is its own
    positions_neighboring_communities[0] = c.get_community(node);
    neighbor_weight[positions_neighboring_communities[0]] = 0;
    neighboring_communities = 1;

    for (unsigned int i = 0; i < deg; ++i) {
        // Fetching neighbors of i, their community and the corresponding degrees
        unsigned int neigh = (c.community_graph)->get_out_neighbor(p + i);
        int neigh_comm = c.get_community(neigh);
        double neigh_w = ((c.community_graph)->is_weighted()) ? (c.community_graph)->get_weighted_out_neighbor(p + i) : 1.f;

        if (neigh != node) {
            // If the community is discovered for the first time
            if (neighbor_weight[neigh_comm] == -1.f) {
                neighbor_weight[neigh_comm] = 0.f;
                positions_neighboring_communities[neighboring_communities++] = neigh_comm;
            }
            // The degree of i toward this community is updated
            neighbor_weight[neigh_comm] += neigh_w;
        }
    }

    // Proceeding similarly on in-neighbors
    size_t p_in = (c.community_graph)->in_neighbors(node);
    unsigned int deg_in = (c.community_graph)->in_degree(node);

    for (unsigned int i = 0; i < deg_in; ++i) {
        unsigned int neigh_in = (c.community_graph)->get_in_neighbor(p_in + i);
        int neigh_comm_in = c.get_community(neigh_in);
        double neigh_w_in = ((c.community_graph)->is_weighted()) ? (c.community_graph)->get_weighted_in_neighbor(p_in + i) : 1.f;

        if (neigh_in != node) {
            if (neighbor_weight[neigh_comm_in] == -1.f) {
                neighbor_weight[neigh_comm_in] = 0.f;
                positions_neighboring_communities[neighboring_communities++] = neigh_comm_in;
            }
            neighbor_weight[neigh_comm_in] += neigh_w_in;
        }
    }
}

