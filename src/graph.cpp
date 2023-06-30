#include <string>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>

#include <sys/mman.h>
#include "../include/graph.hpp"

const unsigned int MAP_LIMIT = 5000000;
static unsigned int build_map(const Graph &g, string, vector<unsigned long>&, vector<vector<pair<unsigned int,double> > >&, vector<vector<pair<unsigned int,double> > >&, bool, bool, bool, bool);

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

Graph::Graph(string filename, bool reproducibility, bool renumbering, bool weighted, bool verbose) {
    vector<vector<pair<unsigned int,double> > > LOUT;
    vector<vector<pair<unsigned int,double> > > LIN;

    this->weighted = weighted;
    string extension = filename.substr(filename.size()-4,filename.size());

    if(extension!=".bin") {
        this->correspondance.resize(0);
        this->nodes = build_map(this, filename, this->correspondance, LOUT, LIN, this->weighted, renumbering, reproducibility, verbose);

        init_attributes(*this, LOUT, LIN, verbose);
        string name = filename.substr(0,filename.size()-4);
        this->write(name+".bin");
    }
    else
        this->load(filename, verbose);
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

void
Graph::write(string outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);
    assert(foutput.rdstate() == ios::goodbit);

    // Writing all attributes of the graph in binary file, out- first then in-
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

void Graph::load(string filename, bool verbose) {
    ifstream finput;
    finput.open(filename, fstream:: in | fstream::binary);
    assert(finput.rdstate() == ios::goodbit);
    this->outcoming_weights.resize(0);
    this->incoming_weights.resize(0);

    finput.read((char*) & this->nodes, sizeof(unsigned int));
    if(verbose) 
        cerr << "number of nodes: " << this->nodes << endl;

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
  
    if(verbose)  
        cerr << "number of arcs:" << this->indegrees[this->nodes - 1] << endl;

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

    if(verbose)
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

// Friend and static functions are defered to a different file for readability 
#include "graph_friend_static.cpp"
