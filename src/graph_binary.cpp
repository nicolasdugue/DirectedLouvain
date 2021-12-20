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

Graph::Graph() {
    nb_nodes = 0;
    nb_links_out = 0;
    nb_links_in = 0;
    total_weight = 0;
}

/* FIXME: rename */
const unsigned int nodes = 1000000;

static void build_map(string filename, vector<ULLI> &correspondance, vector<unsigned int> &corres, map<ULLI, unsigned int> &corres_big_ids, int type);
static void init_attributes(Graph &g, vector<vector<pair<unsigned int,float> > > &LOUT, vector<vector<pair<unsigned int,float> > > &LIN, int type);

Graph::Graph(string in_filename, int type) {
    vector<vector<pair<unsigned int,float> > > LOUT;
    vector<vector<pair<unsigned int,float> > > LIN;

    /* FIXME: this has nothing to do in the constructor */
    float weight = 1.f;
    ofstream foutput;
    string tmp = "";
    /* TODO: split on string to add "_renum" before extension */
    tmp+=in_filename;
    tmp+="_renum";

    cerr << "Renumerotation begins..." << endl;

    correspondance.resize(0);

    //Creates the correspondance table
    vector<unsigned int> corres(nodes,0);
    //Creates the specific table for huge ints that have to be stored as long long int
    map < ULLI, unsigned int > corres_big_ids;

    build_map(in_filename, this->correspondance, corres, corres_big_ids, type);

    cerr << "Renumerotation ends..." << endl;

    // Then we build the graph reading the new graph
    // Out links first
    ifstream finput;
    finput.open(in_filename, fstream:: in );
    weight = 1.f;
    /* FIXME: since the graph is renumbered, we may avoid resizing so much */
    LOUT.resize(this->correspondance.size());
    LIN.resize(this->correspondance.size());

    unsigned int src, dest, map_src, map_dest;
    while (finput >> src >> dest) {
        if (src < nodes) 
            map_src = corres[src] - 1;
        else {
            map_src = corres_big_ids[src]-1;
            //map_src = it_src -> second - 1;
        }

        if (dest < nodes) 
            map_dest = corres[dest] - 1;
        else {
            map_dest = corres_big_ids[dest]-1;
            //map_dest = it_dest -> second - 1;
        }

        if (type == WEIGHTED)
            finput >> weight;

        LOUT[map_src].push_back(make_pair(map_dest, weight));
        LIN[map_dest].push_back(make_pair(map_src, weight));
    }

    cerr << "done." << endl;
    finput.close();

    init_attributes(*this, LIN, LOUT, type);
    /* FIXME: needed only if reproducibility option is chosen */
    //foutput.open(tmp, fstream::out | fstream::binary);
    //foutput.close();
    /*if(finput) {
        unsigned int src, dest;
        while (finput >> src >> dest) {
            unsigned int map_src, map_dest;
            if (src < nodes) 
                map_src = corres[src] - 1;
            else {
                auto it_src = corres_big_ids.find(src);
                map_src = it_src -> second - 1;
            }

            if (dest < nodes) 
                map_dest = corres[dest] - 1;
            else {
                auto it_dest = corres_big_ids.find(dest);
                map_dest = it_dest -> second - 1;
            }

            foutput << map_src << " " << map_dest;
            if (type == WEIGHTED)
                foutput << " " << weight;
            foutput << endl;
        }
    }*/

}

/* Reading from binary */
Graph::Graph(string filename, string filename_w, int type) {
    ifstream finput;
    finput.open(filename, fstream:: in | fstream::binary);

    cerr << "number of nodes" << endl;
    // Read number of nodes on 4 bytes
    finput.read((char * ) & nb_nodes, sizeof(size_t));
    assert(finput.rdstate() == ios::goodbit);
    cerr << "done: " << nb_nodes << endl;

    // Read cumulative out degree sequence: 8 bytes for each node
    cerr << "degrees out" << endl;
    degrees_out.resize(nb_nodes);
    finput.read((char * ) & degrees_out[0], nb_nodes * sizeof(unsigned long int));
    cerr << "done : " << degrees_out[nb_nodes - 1] << endl;

    // Read links: 4 bytes for each link
    nb_links_out = degrees_out[nb_nodes - 1];
    links.resize(nb_links_out);
    finput.read((char * )( & links[0]), (long) nb_links_out * sizeof(unsigned int));

    // Read cumulative in degree sequence: 8 bytes for each node
    cerr << "degrees in" << endl;
    degrees_in.resize(nb_nodes);
    finput.read((char * ) & degrees_in[0], nb_nodes * sizeof(unsigned long int));
    cerr << "done : " << degrees_in[nb_nodes - 1] << endl;

    // Read links_in: 4 bytes for each link
    nb_links_in = degrees_in[nb_nodes - 1];
    links_in.resize(nb_links_in);
    finput.read((char * )( & links_in[0]), (long) nb_links_in * sizeof(unsigned int));

    // Read correspondance of labels
    cerr << "correspondance" << endl;
    correspondance.resize(nb_nodes);
    finput.read((char * )( & correspondance[0]), nb_nodes * sizeof(unsigned long long int));

    // if weighted, read weights: 4 bytes for each link (each link is counted twice)
    weights.resize(0);
    weights_in.resize(0);

    total_weight = 0;
    if (type == WEIGHTED) {
        cerr << "Weights reading" << endl;
        ifstream finput_w;
        finput_w.open(filename_w, fstream:: in | fstream::binary);
        weights.resize(nb_links_out);
        finput_w.read((char * ) & weights[0], nb_links_out * sizeof(float));
        weights_in.resize(nb_links_in);
        finput_w.read((char * ) & weights_in[0], nb_links_in * sizeof(float));
        cerr << "Done" << endl;
    }

    // Compute total weight
    for (unsigned int i = 0; i < nb_nodes; i++) {
        total_weight += out_weighted_degree(i);
    }
}

void
Graph::display() {
    for (unsigned int node = 0; node < nb_nodes; node++) {
        pair < size_t, size_t > p = neighbors(node);
        cout << node << ":";
        for (unsigned int i = 0; i < nb_neighbors_out(node); i++) {
            if (true) {
                if (weights.size() != 0)
                    cout << " (" << links[p.first + i] << " " << weights[p.second + i] << ")";
                else
                    cout << " " << links[p.first+i];
            }
        }
        cout << endl;
    }
}

void
Graph::display_binary(char * outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);

    foutput.write((char * )( & nb_nodes), 4);
    foutput.write((char * )( & degrees_out[0]), 4 * nb_nodes);
    foutput.write((char * )( & links[0]), 8 * nb_links_out);
    foutput.write((char * )( & degrees_in[0]), 4 * nb_nodes);
    foutput.write((char * )( & links_in[0]), 8 * nb_links_in);
}

static void build_map(string filename, vector<ULLI> &correspondance, vector<unsigned int> &corres, map<ULLI, unsigned int> &corres_big_ids, int type) {

    ifstream finput;
    finput.open(filename, fstream:: in );
    /* FIXME: ugly trick, this starts at 1 to say "if corres[node] == 0 then it has not be assigned yet" */
    ULLI cpt = 1;
    float weight = 1.f;
    if (finput) {
        unsigned int src, dest;

        /* We first do the renumerotation and build the correspondance */
        while (finput >> src >> dest) {

            if (type == WEIGHTED)
                finput >> weight;
            //If src is a long that can be stored as an int
            if (src < nodes) {
                if (corres[src] == 0) {
                    corres[src] = cpt++;
                    correspondance.push_back(src);
                }
            } else {
                corres_big_ids[src] = cpt++;
                correspondance.push_back(src);
            }

            if (dest < nodes) {
                if (corres[dest] == 0) {
                    corres[dest] = cpt++;
                    correspondance.push_back(dest);
                }
            } else {
                corres_big_ids[dest] = cpt++;
                correspondance.push_back(dest);
            }
        }
    }
    finput.close();
}

static void init_attributes(Graph &g, vector<vector<pair<unsigned int,float> > > &LOUT, vector<vector<pair<unsigned int,float> > > &LIN, int type) {
    g.nb_nodes = LOUT.size();
    cerr << "number of nodes: " << g.nb_nodes << endl;

    cerr << "degrees out:";
    g.degrees_out.resize(g.nb_nodes);
    unsigned long int tot = 0;
    for (size_t i = 0; i < g.nb_nodes; i++) {
        tot += (unsigned long int) LOUT[i].size();
        g.degrees_out[i] = tot;
    }
    g.nb_links_out = g.degrees_out[g.nb_nodes - 1];
    g.links.resize(g.nb_links_out);
    cerr << g.nb_links_out << endl;

    unsigned long int total_LOUT = 0;
    for (size_t i = 0; i < g.nb_nodes; i++) {
        for (auto edge : LOUT[i]) {
            g.links[total_LOUT]=edge.first;
            ++total_LOUT;
        }
    }

    // Release memory
    for (size_t i = 0; i < LOUT.size(); i++) {
        LOUT[i].clear();
        vector < pair < unsigned int, float > > ().swap(LOUT[i]);
    }

    LOUT.clear();
    vector < vector < pair < unsigned int, float > > > ().swap(LOUT);

    cerr << "degrees in:";
    // Read cumulative in degree sequence: 8 bytes for each node
    g.degrees_in.resize(g.nb_nodes);
    tot = 0;
    for (size_t i = 0; i < g.nb_nodes; i++) {
        tot += (unsigned long int) LIN[i].size();
        g.degrees_in[i] = tot;
    }
    g.nb_links_in = g.degrees_in[g.nb_nodes - 1];
    g.links_in.resize(g.nb_links_in);
    cerr << g.nb_links_in << endl;


    unsigned long int total_LIN = 0;
    for (size_t i = 0; i < g.nb_nodes; i++) {
        for (auto edge : LIN[i]) {
            g.links_in[total_LIN]=edge.first;
            ++total_LIN;
        }
    }

    // if weighted, read weights: 4 bytes for each link (each link is counted twice)
    if(type==WEIGHTED) {
        g.weights.resize(g.nb_links_out);
        g.weights_in.resize(g.nb_links_in);

        for (size_t i = 0; i < g.nb_nodes; i++) {
            for (auto edge : LOUT[i]) 
                g.weights[i]=edge.second;
        }

        for (size_t i = 0; i < g.nb_nodes; i++) {
            for (auto edge : LIN[i]) 
                g.weights_in[i]=edge.second;
        }
    }
    else {
        g.weights.resize(0);
        g.weights_in.resize(0);
    }

    // Compute total weight
    g.total_weight = 0.;
    for (unsigned int i = 0; i < g.nb_nodes; i++) {
        g.total_weight += g.out_weighted_degree(i);
    }
    cerr << "total weight: " << g.total_weight << endl;
}

