// File: graph_binary.cpp
// -- graph handling source
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

#include <sys/mman.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>
#include "../include/graph_binary.hpp"

/* FIXME: rename */
const unsigned int nodes = 10000000;
static unsigned int build_map(string filename, vector<ULI> &correspondance, vector<int> &corres, map<ULI, unsigned int> &corres_big_ids, int type, bool renumbering);

Graph::Graph() {
    nb_nodes = 0;
    nb_links_out = 0;
    nb_links_in = 0;
    total_weight = 0;
}

Graph::Graph(string in_filename, short type, bool reproducibility, bool renumbering) {
    vector<vector<pair<unsigned int,float> > > LOUT;
    vector<vector<pair<unsigned int,float> > > LIN;

    this->type = type;

    float weight = 1.f;

    string extension = in_filename.substr(in_filename.size()-4,in_filename.size());

    if(extension!=".bin") {
        correspondance.resize(0);

        //Creates the correspondance table
        vector<int> corres(nodes,-1);
        //Creates the specific table for huge ints that have to be stored as long long int
        map < ULI, unsigned int > corres_big_ids;
        this->nb_nodes = build_map(in_filename, this->correspondance, corres, corres_big_ids, this->type, renumbering);

        /* FIXME: since the graph is renumbered, we may avoid resizing so much */
        LOUT.resize(this->correspondance.size());
        LIN.resize(this->correspondance.size());

        ifstream finput;
        finput.open(in_filename, fstream:: in );
        ofstream foutput;
        string name = in_filename.substr(0,in_filename.size()-4);
        unsigned int src, dest, map_src, map_dest;
        if(reproducibility) {
            /* TODO: split on string to add "_renum" before extension */
            string tmp=name+"_renum";
            tmp+=extension;
            foutput.open(tmp, fstream::out | fstream::binary);
        }
        while (finput >> src >> dest) {
            weight = 1.f;
            if (type == WEIGHTED)
                finput >> weight;

            if (src < nodes) 
                map_src = corres[src];
            else 
                map_src = corres_big_ids[src];

            if (dest < nodes) 
                map_dest = corres[dest];
            else 
                map_dest = corres_big_ids[dest];

            LOUT[map_src].push_back(make_pair(map_dest, weight));
            LIN[map_dest].push_back(make_pair(map_src, weight));

            if(reproducibility) {
                foutput << map_src << " " << map_dest;
                if (type == WEIGHTED)
                    foutput << " " << weight;
                foutput << endl;
            }
        }

        if(reproducibility) 
            foutput.close();

        finput.close();
        init_attributes(*this, LOUT, LIN);

        if(reproducibility) 
            this->display_binary(name+".bin");
        
    }
    else
        this->load_from_binary(in_filename);
}

Graph::Graph(const Graph &g) {
    this->type = g.type; 

    this->links = g.links;
    this->links_in = g.links_in;
    this->degrees_out = g.degrees_out;
    this->degrees_in = g.degrees_in;
    this->weights = g.weights;
    this->weights_in = g.weights_in;
    this->correspondance = g.correspondance;

    this->nb_nodes = g.nb_nodes;
    this->nb_links_out = g.nb_links_out;
    this->nb_links_in = g.nb_links_in;
    this->total_weight = g.total_weight;
}


/* Reading from binary */
void Graph::load_from_binary(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in | fstream::binary);
    weights.resize(0);
    weights_in.resize(0);

    cerr << "number of nodes: ";
    // Read number of nodes on 4 bytes
    finput.read((char * ) & nb_nodes, sizeof(unsigned int));
    assert(finput.rdstate() == ios::goodbit);
    cerr << nb_nodes << endl;

    // Read cumulative out degree sequence: 8 bytes for each node
    cerr << "total degrees out: ";
    degrees_out.resize(nb_nodes);
    finput.read((char * ) & degrees_out[0], nb_nodes * sizeof(unsigned long));
    cerr << degrees_out[nb_nodes - 1] << endl;

    /* FIXME: really the same thing? */
    this->total_weight = degrees_out[nb_nodes-1];

    // Read links: 4 bytes for each link
    nb_links_out = degrees_out[nb_nodes - 1];
    links.resize(nb_links_out);
    finput.read((char * )( & links[0]), (long) nb_links_out * sizeof(unsigned int));

    if (this->type == WEIGHTED) {
        weights.resize(nb_links_out);
        finput.read((char * ) & weights[0], nb_links_out * sizeof(float));
    }

    // Read cumulative in degree sequence: 8 bytes for each node
    cerr << "total degrees in: ";
    degrees_in.resize(nb_nodes);
    finput.read((char * ) & degrees_in[0], nb_nodes * sizeof(unsigned long));
    cerr << degrees_in[nb_nodes - 1] << endl;

    // Read links_in: 4 bytes for each link
    nb_links_in = degrees_in[nb_nodes - 1];
    links_in.resize(nb_links_in);
    finput.read((char * )( & links_in[0]), (long) nb_links_in * sizeof(unsigned int));

    if (type == WEIGHTED) {
        weights_in.resize(nb_links_in);
        finput.read((char * ) & weights_in[0], nb_links_in * sizeof(float));
    }

    // Read correspondance of labels
    correspondance.resize(nb_nodes);
    finput.read((char * )( & correspondance[0]), nb_nodes * sizeof(ULI));

    /* FIXME: why do we compute once again the total weighted out degree? 
    for (unsigned int i = 0; i < nb_nodes; ++i) {
        total_weight += out_weighted_degree(i);
    }*/
}

void
Graph::display() const {
    for (unsigned int node = 0; node < nb_nodes; node++) {
        pair < size_t, size_t > p = neighbors(node);
        cout << this->correspondance[node] << ":";
        for (unsigned int i = 0; i < nb_neighbors_out(node); ++i) {
            if (true) {
                if (weights.size() != 0)
                    cout << " (" << this->correspondance[links[p.first + i]] << " " << weights[p.second + i] << ")";
                else
                    cout << " " << this->correspondance[links[p.first+i]];
            }
        }
        cout << endl;
    }
}

void
Graph::display_binary(string outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);

    foutput.write((char * )( & nb_nodes), sizeof(unsigned int));
    foutput.write((char * )( & degrees_out[0]), sizeof(unsigned long) * nb_nodes);
    foutput.write((char * )( & links[0]), sizeof(unsigned int) * nb_links_out);
    if(this->type==WEIGHTED)
        foutput.write((char * )( & weights[0]), sizeof(float) * nb_links_out);
    foutput.write((char * )( & degrees_in[0]), sizeof(unsigned long) * nb_nodes);
    foutput.write((char * )( & links_in[0]), sizeof(unsigned int) * nb_links_in);
    if(this->type==WEIGHTED)
        foutput.write((char * )( & weights_in[0]), sizeof(float) * nb_links_in);
    
    for (auto c : this->correspondance) {
        foutput.write((char * )( & c), sizeof(ULI));
    }

    foutput.close();
}

static unsigned int build_map(string filename, vector<ULI> &correspondance, vector<int> &corres, map<ULI, unsigned int> &corres_big_ids, int type, bool renumbering) {

    if(renumbering)
        cerr << "renumbering graph..." << endl;
    ifstream finput;
    finput.open(filename, fstream:: in );
    unsigned int cpt = 0;
    float weight = 1.f;
    if (finput) {
        unsigned int src, dest;

        while (finput >> src >> dest) {

            if (type == WEIGHTED)
                finput >> weight;
            //If src is a long that can be stored as an int
            if (src < nodes) {
                if (corres[src] == -1) {
                    corres[src] = cpt++;
                    if(renumbering)
                        correspondance.push_back(src);
                }
            } else {
                if(corres_big_ids.find(src) == corres_big_ids.end()) {
                    corres_big_ids[src] = cpt++;
                    if(renumbering)
                        correspondance.push_back(src);
                }
            }

            if (dest < nodes) {
                if (corres[dest] == -1) {
                    corres[dest] = cpt++;
                    if(renumbering)
                        correspondance.push_back(dest);
                }
            } else {
                if(corres_big_ids.find(dest) == corres_big_ids.end()) {
                    corres_big_ids[dest] = cpt++;
                    if(renumbering)
                        correspondance.push_back(dest);
                }
            }
        }
    }

    /* If the graph is already renumbered the correspondance must be identity */
    if(!renumbering) {
        for(unsigned int i = 0; i < cpt; ++i)
            correspondance.push_back(i);
    }
    else 
        cerr << "done." << endl;

    finput.close();
    return cpt;
}

void init_attributes(Graph &g, vector<vector<pair<unsigned int,float> > > &LOUT, vector<vector<pair<unsigned int,float> > > &LIN) {
    cerr << "number of nodes: " << g.nb_nodes << endl;

    cerr << "degrees out: ";
    g.degrees_out.resize(g.nb_nodes);
    unsigned long int tot = 0;
    for (size_t i = 0; i < g.nb_nodes; ++i) {
        tot += (unsigned long int) LOUT[i].size();
        g.degrees_out[i] = tot;
    }
    g.nb_links_out = g.degrees_out[g.nb_nodes - 1];
    g.links.resize(g.nb_links_out);
    cerr << g.nb_links_out << endl;

    if(g.type==WEIGHTED) {
        g.weights.resize(g.nb_links_out);
    }
    else {
        g.weights.resize(0);
    }

    unsigned long int total_LOUT = 0;
    for (size_t i = 0; i < g.nb_nodes; ++i) {
        for (auto edge : LOUT[i]) {
            g.links[total_LOUT]=edge.first;
            if(g.type==WEIGHTED)
                g.weights[total_LOUT]=edge.second;
            ++total_LOUT;
        }
    }

    // Release memory
    for (size_t i = 0; i < LOUT.size(); ++i) {
        LOUT[i].clear();
        vector < pair < unsigned int, float > > ().swap(LOUT[i]);
    }

    LOUT.clear();
    vector < vector < pair < unsigned int, float > > > ().swap(LOUT);

    cerr << "degrees in: ";
    g.degrees_in.resize(g.nb_nodes);
    tot = 0;
    for (size_t i = 0; i < g.nb_nodes; ++i) {
        tot += (unsigned long int) LIN[i].size();
        g.degrees_in[i] = tot;
    }
    g.nb_links_in = g.degrees_in[g.nb_nodes - 1];
    g.links_in.resize(g.nb_links_in);
    cerr << g.nb_links_in << endl;

    if(g.type==WEIGHTED) {
        g.weights_in.resize(g.nb_links_in);
    }
    else {
        g.weights_in.resize(0);
    }

    unsigned long int total_LIN = 0;
    for (size_t i = 0; i < g.nb_nodes; ++i) {
        for (auto edge : LIN[i]) {
            g.links_in[total_LIN]=edge.first;
            if(g.type==WEIGHTED)
                g.weights_in[total_LIN]=edge.second;
            ++total_LIN;
        }
    }

    // Compute total weight
    g.total_weight = g.degrees_out[g.nb_nodes-1];
    cerr << "total weight: " << g.total_weight << endl;
}

