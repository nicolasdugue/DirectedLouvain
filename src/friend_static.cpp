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

    /* Based on https://stackoverflow.com/questions/15115943/what-is-the-best-efficient-way-to-read-millions-of-integers-separated-by-lines-f */
    FILE *f = fopen(filename.c_str(), "r");

    unsigned int cpt = 0;
    int hasnum = 0;
    int num = 0;
    int bytes = 0;
    char buffer[8192];
    int eof = 0;
    char *p = NULL;
    while(!eof) {
        fread(buffer, 1, sizeof(buffer), f);
        p = buffer;
        bytes = 8192;
        while(bytes > 0) {
            if (*p == 26) {
                eof = 1;
                break;
            }
            if (*p >= '0' &&  *p <= '9') {
                hasnum = 1;
                num *= 10;
                num += *p-'0';
                ++p;
                --bytes;
            }
            else if (*p == ' ') {
                if (hasnum) 
                    add_to_map(num, cpt, correspondance, corres, corres_big_ids, renumbering);
                num = 0;
                ++p;
                --bytes;
                hasnum = 0;
            }
            else if (*p == '\n') {
                /* skipping weight for now */
                if(!weighted)
                    if (hasnum) 
                        add_to_map(num, cpt, correspondance, corres, corres_big_ids, renumbering);
                num = 0;
                ++p;
                --bytes;
                hasnum = 0;
            }
            else {
                cerr << "Error reading graph." << endl;
                exit(1);
            }
        }
        memset(buffer, 26, sizeof(buffer));  // To detect end of files. 
    }
    p = NULL;
    fclose(f);

    /* Reading the file again to avoid too many resize for LOUT and LIN */
    f = fopen(filename.c_str(), "r");

    LOUT.resize(cpt);
    LIN.resize(cpt);

    double weight = 1.f;
    unsigned int map_src, map_dest;

    size_t i = 0;
    bytes = 0;
    hasnum = 0;
    num = 0;
    unsigned long arc[2];
    eof = 0;
    while(!eof) {
        fread(buffer, 1, sizeof(buffer), f);
        p = buffer;
        bytes = 8192;
        while(bytes > 0) {
            if (*p == 26) {
                eof = 1;
                break;
            }
            /* FIXME: deal with double weights */
            if (*p >= '0' &&  *p <= '9') {
                hasnum = 1;
                num *= 10;
                num += *p-'0';
                ++p;
                --bytes;
            }
            else if (*p == ' ') {
                if (hasnum) 
                    arc[i++] = num;
                num = 0;
                ++p;
                --bytes;
                hasnum = 0;
            }
            else if(*p == '\n') {
                /* FIXME: need to deal with double weights */
                /* src and dest have been read */
                if(weighted) 
                    weight = num;
                /* otherwise we need to read dest */
                else {
                    if (hasnum) 
                        arc[i++] = num;
                }
                num = 0;
                i = 0;
                ++p;
                --bytes;
                hasnum = 0;
                if(renumbering) {
                    map_src = get_mapped_node(arc[0], corres, corres_big_ids);
                    map_dest = get_mapped_node(arc[1], corres, corres_big_ids);
                }
                else {
                    map_src = arc[0];
                    map_dest = arc[1];
                }

                LOUT[map_src].push_back(make_pair(map_dest, weight));
                if(map_src!=map_dest) 
                    LIN[map_dest].push_back(make_pair(map_src, weight));

                /* FIXME: find a faster way */
                if(reproducibility) {
                    foutput << map_src << " " << map_dest;
                    if (weighted)
                        foutput << " " << weight;
                    foutput << endl;
                }
            }
            else {
                cout << "Error..." << endl;
                exit(1);
            }
        }
        memset(buffer, 26, sizeof(buffer));  // To detect end of files. 
    }
    p = NULL;
    fclose(f);

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
