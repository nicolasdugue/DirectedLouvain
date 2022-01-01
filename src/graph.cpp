// File: graph.cpp
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
#include "../include/graph.hpp"

/* FIXME: rename */
const unsigned int MAP_LIMIT = 10000000;
static unsigned int build_map(string filename, vector<ULI> &correspondance, vector<int> &corres, map<ULI, unsigned int> &corres_big_ids, int type, bool renumbering);

Graph::Graph() {
    this->nodes = 0;
    this->arcs = 0;
    this->total_weight = 0;
}

Graph::Graph(string in_filename, short type, bool reproducibility, bool renumbering) {
    vector<vector<pair<unsigned int,double> > > LOUT;
    vector<vector<pair<unsigned int,double> > > LIN;

    this->type = type;

    double weight = 1.f;

    string extension = in_filename.substr(in_filename.size()-4,in_filename.size());

    if(extension!=".bin") {
        /* FIXME: handle wrong filename: add exceptions! */
        correspondance.resize(0);

        //Creates the correspondance table
        vector<int> corres(MAP_LIMIT,-1);
        //Creates the specific table for huge ints that have to be stored as long long int
        map < ULI, unsigned int > corres_big_ids;
        this->nodes = build_map(in_filename, this->correspondance, corres, corres_big_ids, this->type, renumbering);

        LOUT.resize(this->correspondance.size());
        LIN.resize(this->correspondance.size());

        ifstream finput;
        finput.open(in_filename, fstream:: in );
        ofstream foutput;
        string name = in_filename.substr(0,in_filename.size()-4);
        unsigned int src, dest, map_src, map_dest;
        /* FIXME: move into build_map to avoid double reading */
        if(reproducibility) {
            string tmp=name+"_renum";
            tmp+=extension;
            foutput.open(tmp, fstream::out | fstream::binary);
        }
        cerr << "initializing graph..." << endl;
        while (finput >> src >> dest) {
            weight = 1.f;
            if (type == WEIGHTED)
                finput >> weight;

            if (src < MAP_LIMIT) 
                map_src = corres[src];
            else 
                map_src = corres_big_ids[src];

            if (dest < MAP_LIMIT) 
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
        /* FIXME: actually writing directly into file then reading from it seems faster */
        init_attributes(*this, LOUT, LIN);

        if(reproducibility) 
            this->write(name+".bin");
        cerr << "done." << endl;
        
    }
    else
        this->load(in_filename);
}

Graph::Graph(const Graph &g) {
    this->type = g.type; 

    this->nodes = g.nodes;
    this->arcs = g.arcs;
    this->total_weight = g.total_weight;

    this->outcoming_arcs = g.outcoming_arcs;
    this->incoming_arcs = g.incoming_arcs;
    this->outdegrees = g.outdegrees;
    this->indegrees = g.indegrees;
    this->outcoming_weights = g.outcoming_weights;
    this->incoming_weights = g.incoming_weights;
    this->correspondance = g.correspondance;

}

/* FIXME: do not use vector if reproducibility, but instead write _on the fly_ and then read */
void
Graph::write(string outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);

    foutput.write((char * )( & this->nodes), sizeof(unsigned int));
    foutput.write((char * )( & outdegrees[0]), sizeof(unsigned long) * this->nodes);
    foutput.write((char * )( & outcoming_arcs[0]), sizeof(unsigned int) * arcs);
    if(this->type==WEIGHTED)
        foutput.write((char * )( & outcoming_weights[0]), sizeof(double) * arcs);
    foutput.write((char * )( & indegrees[0]), sizeof(unsigned long) * this->nodes);
    foutput.write((char * )( & incoming_arcs[0]), sizeof(unsigned int) * arcs);
    if(this->type==WEIGHTED)
        foutput.write((char * )( & incoming_weights[0]), sizeof(double) * arcs);
    
    for (auto c : this->correspondance) {
        foutput.write((char * )( & c), sizeof(ULI));
    }

    foutput.close();
}

/* Reading from binary */
void Graph::load(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in | fstream::binary);
    outcoming_weights.resize(0);
    incoming_weights.resize(0);

    cerr << "number of nodes: ";
    finput.read((char * ) & this->nodes, sizeof(unsigned int));
    assert(finput.rdstate() == ios::goodbit);
    cerr << this->nodes << endl;

    outdegrees.resize(this->nodes);
    finput.read((char * ) & outdegrees[0], this->nodes * sizeof(unsigned long));

    arcs = outdegrees[this->nodes - 1];
    outcoming_arcs.resize(arcs);
    finput.read((char * )( & outcoming_arcs[0]), (long) arcs * sizeof(unsigned int));

    if (this->type == WEIGHTED) {
        outcoming_weights.resize(arcs);
        finput.read((char * ) & outcoming_weights[0], arcs * sizeof(double));
    }

    indegrees.resize(this->nodes);
    finput.read((char * ) & indegrees[0], this->nodes * sizeof(unsigned long));
    
    cerr << "number of arcs: ";
    cerr << indegrees[this->nodes - 1] << endl;

    incoming_arcs.resize(arcs);
    finput.read((char * )( & incoming_arcs[0]), (long) arcs * sizeof(unsigned int));

    if (type == WEIGHTED) {
        incoming_weights.resize(arcs);
        finput.read((char * ) & incoming_weights[0], arcs * sizeof(double));
    }

    correspondance.resize(this->nodes);
    finput.read((char * )( & correspondance[0]), this->nodes * sizeof(ULI));

    this->total_weight = 0.f;
    for (unsigned int i = 0; i < this->nodes; ++i) {
        this->total_weight += weighted_out_degree(i);
    }
    finput.close();
}

void Graph::display() const {
    for (unsigned int node = 0; node < this->nodes; node++) {
        pair < size_t, size_t > p = neighbors(node);
        cout << this->correspondance[node] << ":";
        for (unsigned int i = 0; i < out_degree(node); ++i) {
            if (true) {
                if (outcoming_weights.size() != 0)
                    cout << " (" << this->correspondance[outcoming_arcs[p.first + i]] << " " << outcoming_weights[p.second + i] << ")";
                else
                    cout << " " << this->correspondance[outcoming_arcs[p.first+i]];
            }
        }
        cout << endl;
    }
}

static unsigned int build_map(string filename, vector<ULI> &correspondance, vector<int> &corres, map<ULI, unsigned int> &corres_big_ids, int type, bool renumbering) {

    if(renumbering)
        cerr << "renumbering graph..." << endl;
    ifstream finput;
    finput.open(filename, fstream:: in );
    unsigned int cpt = 0;
    double weight = 1.f;
    if (finput) {
        unsigned int src, dest;

        while (finput >> src >> dest) {

            if (type == WEIGHTED)
                finput >> weight;
            //If src is a long that can be stored as an int
            if (src < MAP_LIMIT) {
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

            if (dest < MAP_LIMIT) {
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

/* FIXME: is this really faster than writing into/reading from a binary file?
 * FIXME: if this is not the case, then remove the reproducibility attribute !
 */
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

    if(g.type==WEIGHTED) 
        g.outcoming_weights.resize(g.arcs);
    else 
        g.outcoming_weights.resize(0);

    unsigned long int total_LOUT = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        for (auto edge : LOUT[i]) {
            g.outcoming_arcs[total_LOUT]=edge.first;
            if(g.type==WEIGHTED)
                g.outcoming_weights[total_LOUT]=edge.second;
            ++total_LOUT;
        }
    }

    // Release memory
    for (size_t i = 0; i < LOUT.size(); ++i) {
        LOUT[i].clear();
        vector < pair < unsigned int, double > > ().swap(LOUT[i]);
    }

    LOUT.clear();
    vector < vector < pair < unsigned int, double > > > ().swap(LOUT);

    g.indegrees.resize(g.nodes);
    tot = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        tot += (unsigned long int) LIN[i].size();
        g.indegrees[i] = tot;
    }
    g.incoming_arcs.resize(g.arcs);

    if(g.type==WEIGHTED) 
        g.incoming_weights.resize(g.arcs);
    else 
        g.incoming_weights.resize(0);

    unsigned long int total_LIN = 0;
    for (size_t i = 0; i < g.nodes; ++i) {
        for (auto edge : LIN[i]) {
            g.incoming_arcs[total_LIN]=edge.first;
            if(g.type==WEIGHTED)
                g.incoming_weights[total_LIN]=edge.second;
            ++total_LIN;
        }
    }

    cerr << "number of arcs: ";
    cerr << g.arcs << endl;

    // Compute total weight
    g.total_weight = 0.;
    for (unsigned int i = 0; i < g.nodes; ++i) {
        g.total_weight += g.weighted_out_degree(i);
    }

    cerr << "total weight: " << g.total_weight << endl;
}

