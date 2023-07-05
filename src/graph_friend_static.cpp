#include<sstream>

// This method adds node to the small or large int correspondance "map"
// A reference to a counter is given, which will be the number of nodes in the end
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

// This function returns the renumbered identifier of a given node
static inline unsigned int get_mapped_node(unsigned long int node, const vector<int>& corres, const map< unsigned long, unsigned int> &corres_big_ids, bool renumbering) {
    unsigned int mapped_node = node;
    if(renumbering) {
        if (node < MAP_LIMIT)
            mapped_node = corres[node];
        else
            mapped_node = corres_big_ids.at(node);
    }
    return mapped_node;
}

// This function builds the map for renumbering the input graph and returns cpt (the number of nodes)
// If reproducibility is set to true, the renumbered graph is written into a file under edgelist format: src dest (weight)
static unsigned int build_map(Graph &g, string filename, vector<unsigned long> &correspondance, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN, bool renumbering, bool reproducibility, bool verbose) {

    vector<int> corres(MAP_LIMIT,-1);
    map < unsigned long, unsigned int > corres_big_ids;
    if(renumbering && verbose)
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
    assert(finput.rdstate() == ios::goodbit);

    unsigned int cpt = 0;

    // Read the graph file to generate a map of node
    string line;
    while (getline(finput, line)) {
        unsigned int src = 0; 
        unsigned int dest = 0; 
        size_t number_of_tokens = 0;
       
        char* line_to_split = new char[line.size()+1];
        strcpy(line_to_split, line.c_str());

        char* token = strtok(line_to_split, " \t");

        while(token != NULL && number_of_tokens<2) {
            if(number_of_tokens==0)
                src = atoi(token);
            if(number_of_tokens==1)
                dest = atoi(token);
            number_of_tokens++;
            token = strtok(NULL, " ");
        }

        assert(number_of_tokens == 2);

        delete[] line_to_split;

        add_to_map(src, cpt, correspondance, corres, corres_big_ids, renumbering);
        add_to_map(dest, cpt, correspondance, corres, corres_big_ids, renumbering);
    }


    // If the graph is already renumbered the correspondance must be identity
    if(!renumbering) {
        // Number of nodes in that case is cpt+1
        ++cpt;
        for(unsigned int i = 0; i < cpt; ++i)
            correspondance.push_back(i);
    }

    LOUT.resize(cpt);
    LIN.resize(cpt);

    finput.clear();
    finput.seekg(0);
   
    while (getline(finput, line)) {
        unsigned int src = 0;
        unsigned int dest = 0; 
        unsigned int map_src, map_dest;
        double weight = 1.f;
        size_t number_of_tokens = 0;
       
        char* line_to_split = new char[line.size()+1];
        strcpy(line_to_split, line.c_str());

        char* token = strtok(line_to_split, " \t");

        while(token != NULL && number_of_tokens<2) {
            if(number_of_tokens==0)
                src = atoi(token);
            if(number_of_tokens==1)
                dest = atoi(token);
            number_of_tokens++;
            token = strtok(NULL, " ");
        }

        if (token !=NULL) {
            
            if(!g.is_weighted()) 
                g.set_weighted(true);

            weight = atof(token);
        }

        delete[] line_to_split;

        map_src = get_mapped_node(src, corres, corres_big_ids, renumbering);
        map_dest = get_mapped_node(dest, corres, corres_big_ids, renumbering);

        LOUT[map_src].push_back(make_pair(map_dest, weight));
        LIN[map_dest].push_back(make_pair(map_src, weight));

        if(reproducibility) {
            foutput << map_src << " " << map_dest;
            if (g.is_weighted())
                foutput << " " << weight;
            foutput << endl;
        }
    }

    if(reproducibility)
        foutput.close();

    return cpt;
}

// This method initializes all attributs of Graph object g from out- and in- adjacency lists LOUT and LIN
void init_attributes(Graph &g, vector<vector<pair<unsigned int,double> > > &LOUT, vector<vector<pair<unsigned int,double> > > &LIN, bool verbose) {
    if(verbose) {
        cerr << "initializing graph..." << endl;
        cerr << "number of nodes: " << g.nodes << endl;
    }

    // Building cumulative out-degree sequence
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

    // Stocking out-neighbors and weights (if any)
    unsigned long int total_LOUT = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        for (auto edge : LOUT[i]) {
            g.outcoming_arcs[total_LOUT]=edge.first;
            if(g.weighted)
                g.outcoming_weights[total_LOUT]=edge.second;
            ++total_LOUT;
        }
    }

    // Releasing memory for efficiency purposes
    for (size_t i = 0; i < LOUT.size(); ++i) {
        LOUT[i].clear();
        vector < pair < unsigned int, double > > ().swap(LOUT[i]);
    }

    LOUT.clear();
    vector < vector < pair < unsigned int, double > > > ().swap(LOUT);

    // Building cumulative in-degree sequence
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

    // Stocking in-neighbors and weights (if any)
    unsigned long int total_LIN = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        for (auto edge : LIN[i]) {
            g.incoming_arcs[total_LIN]=edge.first;
            if(g.weighted)
                g.incoming_weights[total_LIN]=edge.second;
            ++total_LIN;
        }
    }

    // Releasing memory for efficiency purposes
    for (size_t i = 0; i < LIN.size(); ++i) {
        LIN[i].clear();
        vector < pair < unsigned int, double > > ().swap(LIN[i]);
    }

    LIN.clear();
    vector < vector < pair < unsigned int, double > > > ().swap(LIN);

    if(verbose)
        cerr << "number of arcs: " << g.arcs << endl;

    // Computing the total weight of the graph
    g.total_weight = 0.;

    for (unsigned int i = 0; i < g.nodes; ++i) {
        g.total_weight += g.weighted_out_degree(i);
    }

    if(verbose)
        cerr << "total weight: " << g.total_weight << endl << "done." << endl;
}
