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

#include <string>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>

#include <sys/mman.h>
#include "../include/graph.hpp"

const unsigned int MAP_LIMIT = 5000000;
static unsigned int build_map(string, vector<unsigned long>&, vector<vector<pair<unsigned int,double> > >&, vector<vector<pair<unsigned int,double> > >&, bool, bool, bool);

Graph::Graph() {
    this->nodes         = 0;
    this->arcs          = 0;
    this->total_weight  = 0;
    this->weighted      = false;

    this->outcoming_arcs.resize(0);
    this->incoming_arcs.resize(0);
    this->outcoming_weights.resize(0);
    this->incoming_weights.resize(0);
    this->outdegrees.resize(0);
    this->indegrees.resize(0);
}

Graph::Graph(string in_filename, bool weighted, bool reproducibility, bool renumbering) {
    vector<vector<pair<unsigned int,double> > > LOUT;
    vector<vector<pair<unsigned int,double> > > LIN;

    this->weighted = weighted;
    string extension = in_filename.substr(in_filename.size()-4,in_filename.size());

    /* FIXME: handle wrong filename: add exceptions! */
    if(extension!=".bin") {
        this->correspondance.resize(0);
        this->nodes = build_map(in_filename, this->correspondance, LOUT, LIN, this->weighted, renumbering, reproducibility);

        cerr << "initializing graph..." << endl;
        init_attributes(*this, LOUT, LIN);
        cerr << "done." << endl;
        string name = in_filename.substr(0,in_filename.size()-4);
        /* TODO: modify README to take into account that this is now done in every case */
        this->write(name+".bin");
    }
    else
        this->load(in_filename);
}

Graph::Graph(const Graph &g) {
    this->weighted          = g.weighted; 

    this->nodes             = g.nodes;
    this->arcs              = g.arcs;
    this->total_weight      = g.total_weight;

    this->outcoming_arcs    = g.outcoming_arcs;
    this->incoming_arcs     = g.incoming_arcs;
    this->outdegrees        = g.outdegrees;
    this->indegrees         = g.indegrees;
    this->outcoming_weights = g.outcoming_weights;
    this->incoming_weights  = g.incoming_weights;
    this->correspondance    = g.correspondance;
}

/* FIXME: do not use vector if reproducibility, but instead write _on the fly_ and then read */
void
Graph::write(string outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);

    foutput.write((char*)( & this->nodes), sizeof(unsigned int));
    foutput.write((char*)( & this->outdegrees[0]), sizeof(unsigned long) * this->nodes);
    foutput.write((char*)( & this->outcoming_arcs[0]), sizeof(unsigned int) * this->arcs);
    if(this->weighted)
        foutput.write((char*)( & this->outcoming_weights[0]), sizeof(double) * this->arcs);
    foutput.write((char*)( & this->indegrees[0]), sizeof(unsigned long) * this->nodes);
    foutput.write((char*)( & this->incoming_arcs[0]), sizeof(unsigned int) * this->arcs);
    if(this->weighted)
        foutput.write((char*)( & this->incoming_weights[0]), sizeof(double) * this->arcs);
    foutput.write((char*)( & this->correspondance[0]), sizeof(unsigned long) * this->nodes);
    foutput.close();
}

/* Reading from binary */
void Graph::load(string filename) {
    ifstream finput;
    finput.open(filename, fstream:: in | fstream::binary);
    this->outcoming_weights.resize(0);
    this->incoming_weights.resize(0);

    cerr << "number of nodes: ";
    finput.read((char*) & this->nodes, sizeof(unsigned int));
    assert(finput.rdstate() == ios::goodbit);
    cerr << this->nodes << endl;

    this->outdegrees.resize(this->nodes);
    finput.read((char*) & this->outdegrees[0], this->nodes * sizeof(unsigned long));

    this->arcs = outdegrees[this->nodes - 1];
    this->outcoming_arcs.resize(this->arcs);
    finput.read((char*)( & this->outcoming_arcs[0]), this->arcs * sizeof(unsigned int));

    if (this->weighted) {
        this->outcoming_weights.resize(arcs);
        finput.read((char*) & this->outcoming_weights[0], this->arcs * sizeof(double));
    }

    this->indegrees.resize(this->nodes);
    finput.read((char*) & this->indegrees[0], this->nodes * sizeof(unsigned long));
    
    cerr << "number of arcs: ";
    cerr << this->indegrees[this->nodes - 1] << endl;

    this->incoming_arcs.resize(this->arcs);
    finput.read((char*)( & this->incoming_arcs[0]), this->arcs * sizeof(unsigned int));

    if (this->weighted) {
        this->incoming_weights.resize(this->arcs);
        finput.read((char*) & this->incoming_weights[0], this->arcs * sizeof(double));
    }

    this->correspondance.resize(this->nodes);
    finput.read((char*)( & this->correspondance[0]), this->nodes * sizeof(unsigned long));

    this->total_weight = 0.f;
    for (unsigned int i = 0; i < this->nodes; ++i) {
        this->total_weight += this->weighted_out_degree(i);
    }

    cerr << "total weight: " << this->total_weight << endl;
    finput.close();
}

void Graph::display() const {
    for (unsigned int node = 0; node < this->nodes; ++node) {
        size_t p = this->out_neighbors(node);
        cout << this->correspondance[node] << ":";
        for (unsigned int i = 0; i < out_degree(node); ++i) {
            if (this->weighted)
                cout << " (" << this->outcoming_arcs[p + i] << " " << this->outcoming_weights[p + i] << ")";
            else
                cout << " " << this->outcoming_arcs[p+i];
        }
        cout << endl;
    }
}

double Graph::count_selfloops(unsigned int node) {
    assert(node<this->nodes);
    size_t p = this->out_neighbors(node);
    for (unsigned int i=0 ; i < this->out_degree(node) ; ++i) {
        if (this->outcoming_arcs[p+i]==node) {
            if (this->weighted)
                return this->outcoming_weights[p+i];
            else 
                return 1.;
        }
    }

    return 0.;
}

double Graph::weighted_out_degree(unsigned int node) {
    assert(node<this->nodes);
    if (!this->weighted)
        return this->out_degree(node);
    else {
        size_t p = this->out_neighbors(node);
        double res = 0;
        for (unsigned int i=0 ; i < this->out_degree(node) ; ++i) 
            res += this->outcoming_weights[p+i];
        return res;
    }
}

double Graph::weighted_in_degree(unsigned int node) {
    assert(node<this->nodes);
    if (!this->weighted)
        return this->in_degree(node);
    else {
        size_t p = this->in_neighbors(node);
        double res = 0;
        for (unsigned int i=0 ; i < this->in_degree(node) ; ++i) 
            res += this->incoming_weights[p+i];
        return res;
    }
}

/* Friend and static functions are defered to a different file for readability */
#include "friend_static.cpp"
