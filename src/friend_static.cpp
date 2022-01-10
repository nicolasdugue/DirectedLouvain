static void add_to_map(unsigned int node, unsigned int &cpt, vector<unsigned long> &correspondance, vector<int> &corres, map<unsigned long, unsigned int> &corres_big_ids, bool renumbering) {
    if(renumbering) {
        if (node < MAP_LIMIT) {
            if (corres[node] == -1) {
                corres[node] = cpt++;
                    correspondance.push_back(node);
            }
        } else {
            if(corres_big_ids.find(node) == corres_big_ids.end()) {
                corres_big_ids[node] = cpt++;
                    correspondance.push_back(node);
            }
        }
    }
    else if (cpt < node)
        cpt = node;
}

static inline unsigned int get_mapped_node(unsigned long int node, const vector<int>& corres, const map< unsigned long, unsigned int> &corres_big_ids) {
    unsigned int mapped_node;
    if (node < MAP_LIMIT)
        mapped_node = corres[node];
    else
        mapped_node = corres_big_ids.at(node);
    return mapped_node;
}

static unsigned int build_map(string filename, vector<unsigned long> &correspondance, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN, bool weighted, bool renumbering, bool reproducibility) {

    //Creates the correspondance table
    vector<int> corres(MAP_LIMIT,-1);
    //Creates the specific table for huge ints that have to be stored as long long int
    map < unsigned long, unsigned int > corres_big_ids;
    if(renumbering)
        cerr << "renumbering graph..." << endl;
    ofstream foutput;
    if(reproducibility) {
        size_t name_size = filename.size();
        string name = filename.substr(0,name_size-4);
        string extension = filename.substr(name_size-4,name_size);
        string tmp=name+"_renum";
        tmp+=extension;
        foutput.open(tmp, fstream::out | fstream::binary);
    }

    ifstream finput;
    finput.open(filename, fstream:: in );
    unsigned int cpt = 0;
    double weight = 1.f;
    if (finput) {
        unsigned int src, dest;

        while (finput >> src >> dest) {
            if (weighted)
                finput >> weight;

            add_to_map(src, cpt, correspondance, corres, corres_big_ids, renumbering);
            add_to_map(dest, cpt, correspondance, corres, corres_big_ids, renumbering);
        }
    }

    LOUT.resize(cpt);
    LIN.resize(cpt);

    weight = 1.f;
    unsigned int src, dest, map_src, map_dest;

    finput.clear();
    finput.seekg(0);
    
    while (finput >> src >> dest) {
            weight = 1.f;
            if (weighted)
                finput >> weight;

            map_src = get_mapped_node(src, corres, corres_big_ids);
            map_dest = get_mapped_node(dest, corres, corres_big_ids);

            LOUT[map_src].push_back(make_pair(map_dest, weight));
            if(map_src!=map_dest)
                LIN[map_dest].push_back(make_pair(map_src, weight));

            if(reproducibility) {
                foutput << map_src << " " << map_dest;
                if (weighted)
                    foutput << " " << weight;
                foutput << endl;
            }
        }

        if(reproducibility)
            foutput.close();

    /* If the graph is already renumbered the correspondance must be identity */
    if(!renumbering) {
        /* Number of nodes is cpt+1 */
        ++cpt;
        for(unsigned int i = 0; i < cpt; ++i)
            correspondance.push_back(i);
    }

    if(reproducibility)
        foutput.close();

    return cpt;
}


void init_attributes(Graph &g, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN) {
    cerr << "number of nodes: " << g.nodes << endl;

    g.outdegrees.resize(g.nodes);
    unsigned long int tot = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        tot += (unsigned long int) LOUT[i].size();
        g.outdegrees[i] = tot;
    }
    g.arcs = g.outdegrees[g.nodes - 1];
    g.outcoming_arcs.resize(g.arcs);

    if(g.weighted) 
        g.outcoming_weights.resize(g.arcs);
    else 
        g.outcoming_weights.resize(0);

    unsigned long int total_LOUT = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        for (auto edge : LOUT[i]) {
            g.outcoming_arcs[total_LOUT]=edge.first;
            if(g.weighted)
                g.outcoming_weights[total_LOUT]=edge.second;
            ++total_LOUT;
        }
    }

    // Release memory
    /*for (size_t i = 0; i < LOUT.size(); ++i) {
        LOUT[i].clear();
        vector < pair < unsigned int, double > > ().swap(LOUT[i]);
    }

    LOUT.clear();
    vector < vector < pair < unsigned int, double > > > ().swap(LOUT);*/

    g.indegrees.resize(g.nodes);
    tot = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        tot += (unsigned long int) LIN[i].size();
        g.indegrees[i] = tot;
    }
    g.incoming_arcs.resize(g.arcs);

    if(g.weighted) 
        g.incoming_weights.resize(g.arcs);
    else 
        g.incoming_weights.resize(0);

    unsigned long int total_LIN = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        for (auto edge : LIN[i]) {
            g.incoming_arcs[total_LIN]=edge.first;
            if(g.weighted)
                g.incoming_weights[total_LIN]=edge.second;
            ++total_LIN;
        }
    }

    // Release memory
    /*for (size_t i = 0; i < LIN.size(); ++i) {
        LIN[i].clear();
        vector < pair < unsigned int, double > > ().swap(LIN[i]);
    }

    LIN.clear();
    vector < vector < pair < unsigned int, double > > > ().swap(LIN);*/

    cerr << "number of arcs: ";
    cerr << g.arcs << endl;

    // Compute total weight
    g.total_weight = 0.;

    double &total_weight = g.total_weight;
    for (unsigned int i = 0; i < g.nodes; ++i) {
        total_weight += g.weighted_out_degree(i);
    }

    cerr << "total weight: " << g.total_weight << endl;
}
