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
        auto p = this->out_neighbors(node);
        cout << this->correspondance[node] << ":";
        for (unsigned int i = 0; i < out_degree(node); ++i) {
            if (this->weighted)
                cout << " (" << this->correspondance[this->outcoming_arcs[p.first + i]] << " " << this->outcoming_weights[p.second + i] << ")";
            else
                cout << " " << this->correspondance[this->outcoming_arcs[p.first+i]];
        }
        cout << endl;
    }
}

double Graph::count_selfloops(unsigned int node) {
    assert(node<this->nodes);
    auto p = this->out_neighbors(node);
    for (unsigned int i=0 ; i < this->out_degree(node) ; ++i) {
        if (this->outcoming_arcs[p.first+i]==node) {
            if (this->weighted)
                return this->outcoming_weights[p.second+i];
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
        auto p = this->out_neighbors(node);
        double res = 0;
        for (unsigned int i=0 ; i < this->out_degree(node) ; ++i) 
            res += this->outcoming_weights[p.second+i];
        return res;
    }
}

double Graph::weighted_in_degree(unsigned int node) {
    assert(node<this->nodes);
    if (!this->weighted)
        return this->in_degree(node);
    else {
        auto p = this->in_neighbors(node);
        double res = 0;
        for (unsigned int i=0 ; i < this->in_degree(node) ; ++i) 
            res += this->incoming_weights[p.second+i];
        return res;
    }
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
    for (size_t i = 0; i < LIN.size(); ++i) {
        LIN[i].clear();
        vector < pair < unsigned int, double > > ().swap(LIN[i]);
    }

    LIN.clear();
    vector < vector < pair < unsigned int, double > > > ().swap(LIN);

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
    int f = open(filename.c_str(), O_RDONLY);
    off_t size = lseek(f, 0, SEEK_END);
    char *buffer = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, f, 0);

    unsigned int cpt = 0;
    int hasnum = 0;
    int num = 0;
    off_t bytes = size;
    char *p = buffer;
    while(bytes > 0) {
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
            /* FIXME: try inserting here in LOUT and LIN using insert and iterator */
        }
        else {
            cerr << "Error reading graph." << endl;
            exit(1);
        }
    }
    close(f);
    munmap(buffer,size);

    /* If the graph is already renumbered the correspondance must be identity */
    if(!renumbering) {
        /* Number of nodes is cpt+1 */
        ++cpt;
        for(unsigned int i = 0; i < cpt; ++i)
            correspondance.push_back(i);
    }

    /* Reading the file again to avoid too many resize for LOUT and LIN */
    f = open(filename.c_str(), O_RDONLY);
    size = lseek(f, 0, SEEK_END);
    buffer = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, f, 0);

    LOUT.resize(cpt);
    LIN.resize(cpt);

    double weight = 1.f;
    unsigned int map_src, map_dest;

    size_t i = 0;
    bytes = size;
    hasnum = 0;
    num = 0;
    unsigned long arc[2];
    p = buffer;
    while(bytes > 0) {
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
    close(f);
    munmap(buffer,size);

    if(reproducibility)
        foutput.close();

    return cpt;
}

