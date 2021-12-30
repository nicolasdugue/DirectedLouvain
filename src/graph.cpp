// File: graph.cpp
// -- simple graph handling source file
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

#include "../include/graph.hpp"
#include <climits>

using namespace std;

/* FIXME: rename */
const unsigned int nodes = 1000000;

static void build_map(string filename, vector<ULLI> &correspondance, vector<unsigned int> &corres, map<ULLI, unsigned int> &corres_big_ids, int type) {

    ifstream finput;
    finput.open(filename, fstream:: in );
    /* FIXME: ugly trick, this starts at 1 to say "if corres[node] == 0 then it has not be assigned yet" */
    ULLI cpt = 1;
    float weight = 1.f;

    cerr << "Start reading" << endl;
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
                if(corres_big_ids.find(src) == corres_big_ids.end()) {
                    corres_big_ids[src] = cpt++;
                    correspondance.push_back(src);
                }
            }

            if (dest < nodes) {
                if (corres[dest] == 0) {
                    corres[dest] = cpt++;
                    correspondance.push_back(dest);
                }
            } else {
                if(corres_big_ids.find(dest) == corres_big_ids.end()) {
                    corres_big_ids[dest] = cpt++;
                    correspondance.push_back(dest);
                }
            }
        }
    }
    finput.close();
}

Graph::Graph(string in_filename, string filename, string filename_w, int type, bool do_renumber) {
    ifstream finput;

    /* FIXME: this has nothing to do in the constructor */
    if (do_renumber) {

        unsigned int nb_links = 0;
        float weight = 1.f;
        ofstream foutput;
        string tmp = "";
        /* TODO: split on string to add "_renum" before extension */
        tmp+=in_filename;
        tmp+="_renum";

        finput.open(in_filename, fstream:: in );
        foutput.open(tmp, fstream::out | fstream::binary);

        cerr << "Renumerotation begins..." << endl;

        correspondance.resize(0);

        //Creates the correspondance table
        vector<unsigned int> corres(nodes,0);
        //Creates the specific table for huge ints that have to be stored as long long int
        map < ULLI, unsigned int > corres_big_ids;

        build_map(in_filename, this->correspondance, corres, corres_big_ids, type);

        if(finput) {
            unsigned int src, dest;
            while (finput >> src >> dest) {
                unsigned int pos_src, pos_dest;
                if (src < nodes) 
                    pos_src = corres[src] - 1;
                else {
                    auto it_src = corres_big_ids.find(src);
                    pos_src = it_src -> second - 1;
                }

                if (dest < nodes) 
                    pos_dest = corres[dest] - 1;
                else {
                    auto it_dest = corres_big_ids.find(dest);
                    pos_dest = it_dest -> second - 1;
                }

                nb_links++;
                // TODO: count nb_links while parsing for the maximum, and pourcent the progression
                foutput << pos_src << " " << pos_dest;
                if (type == WEIGHTED)
                    foutput << " " << weight;
                foutput << endl;
            }
        }

        cerr << "Renumerotation ends..." << endl;
        cerr << "Building the graph... links out" << endl;
        foutput.close();
        finput.close();

        // Then we build the graph reading the new graph
        // Out links first
        finput.open(tmp, fstream:: in );
        weight = 1.f;
        /* FIXME: since the graph is renumbered, we may avoid resizing so much */
        links_out.resize(this->correspondance.size());
        links_in.resize(this->correspondance.size());

        unsigned int src, dest;
        while (finput >> src >> dest) {
            if (type == WEIGHTED)
                finput >> weight;

            links_out[src].push_back(make_pair(dest, weight));
            links_in[dest].push_back(make_pair(src, weight));

            nb_links++;
        }

        cerr << "done." << endl;
        finput.close();

    } else {

        long unsigned int src, dest;
        int nb_links = 0;
        float weight = 1.f;

        finput.open(in_filename, fstream:: in );
        while (finput >> src >> dest) {
            if (type == WEIGHTED)
                finput >> weight;

            if (links_out.size() <= max(src, dest) + 1) 
                links_out.resize(max(src, dest) + 1);
            if (links_in.size() <= max(src, dest) + 1) 
                links_in.resize(max(src, dest) + 1);

            links_out[src].push_back(make_pair(dest, weight));
            links_in[dest].push_back(make_pair(src, weight));

            nb_links++;
        }

        cout << "done." << endl;
        finput.close();

    }
    this->display_binary(filename, filename_w, type, do_renumber);

}

Graph::Graph(const Graph &g) {
    this->links_out = g.links_out;
    this->links_in = g.links_in;
    this->correspondance = g.correspondance;
}

/* This procedure allows for multigraphs: if the input file contains 1 4 3 and 1 4 6 
 * then this procedure will replace (4,3) and (4,6) in the adjacency list of 1 by (4,9)
 */
void
Graph::clean(int type) {
    for (unsigned int i = 0; i < links_out.size(); i++) {
        map < ULLI, float > m;

        for (size_t j = 0; j < links_out[i].size(); j++) {
            auto it = m.find(links_out[i][j].first);
            if (it == m.end())
                m.insert(make_pair(links_out[i][j].first, links_out[i][j].second));
            else if (type == WEIGHTED)
                it -> second += links_out[i][j].second;
        }

        vector < pair < unsigned int, float > > v;
        for (auto it = m.begin(); it != m.end(); it++)
            v.push_back( * it);
        links_out[i].clear();
        links_out[i] = v;
    }
}

/*void
Graph::display(int type) {
    for (unsigned int i = 0; i < links_out.size(); i++) {
        for (unsigned int j = 0; j < links_out[i].size(); j++) {
            int dest = links_out[i][j].first;
            float weight = links_out[i][j].second;
            if (type == WEIGHTED)
                cout << i << " " << dest << " " << weight << endl;
            else
                cout << i << " " << dest << endl;
        }
    }
}*/

void
Graph::display_binary(string filename, string filename_w, int type, bool do_renumber) {
    ofstream foutput;
    foutput.open(filename, fstream::out | fstream::binary);
    /* FIXME: this belongs to the constructor since links_out is released 
     * once written (a better way would be to have a function that writes 
     * by appending (file opened before) and releases memory at the same time */
    size_t s = links_out.size();
    /* /!\FIXME/!\: this is not the number of nodes if not renumbered! */
    // outputs weights in a separate file
    if (type == WEIGHTED) {
        ofstream foutput_w;
        foutput_w.open(filename_w, fstream::out | fstream::binary);
        for (size_t i = 0; i < s; i++) {
            for (auto edge : links_out[i]) {
                float weight = edge.second;
                foutput_w.write((char * )( & weight), sizeof(float));
            }
        }
        s = links_in.size();
        for (size_t i = 0; i < s; i++) {
            for (auto edge : links_in[i]) {
                float weight = edge.second;
                foutput_w.write((char * )( & weight), sizeof(float));
            }
        }
        foutput_w.close();
    }
    cerr << "done" << endl;

    ofstream fbin;
    fbin.open(filename, fstream::out | fstream::binary);
    cerr << "number of nodes : " << s << endl;
    cerr << "writing in binary file..." << endl;
    // outputs number of nodes
    fbin.write((char * )( & s), sizeof(size_t));

    // outputs cumulative degree sequence
    /* Contient uniquement les degres sortants en oriente */
    long tot = 0;
    for (size_t i = 0; i < s; i++) {
        tot += (unsigned long int) links_out[i].size();
        fbin.write((char * )( & tot), sizeof(unsigned long int));
    }

    // outputs links_out
    for (size_t i = 0; i < s; i++) {
        for (auto edge : links_out[i]) {
            unsigned int dest = edge.first;
            fbin.write((char * )( & dest), sizeof(unsigned int));
        }
    }
    cerr << "done." << endl;
    cerr << "releasing memory..." << endl;
    // Release memory
    for (size_t i = 0; i < links_out.size(); i++) {
        links_out[i].clear();
        vector < pair < unsigned int, float > > ().swap(links_out[i]);
    }

    links_out.clear();
    vector < vector < pair < unsigned int, float > > > ().swap(links_out);

    // Writing information --- only in-degrees
    long tot_in = 0;
    for (size_t i = 0; i < s; i++) {
        tot_in += (unsigned long int) links_in[i].size();
        fbin.write((char * )( & tot_in), sizeof(unsigned long int));
    }

    // outputs links_in
    for (size_t i = 0; i < s; i++) {
        for (auto edge : links_in[i]) {
            unsigned int dest = edge.first;
            fbin.write((char * )( & dest), sizeof(unsigned int));
        }
    }
    cerr << "releasing memory..." << endl;
    // Release memory
    for (size_t i = 0; i < links_in.size(); i++) {
        links_in[i].clear();
        vector < pair < unsigned int, float > > ().swap(links_in[i]);
    }

    links_in.clear();
    vector < vector < pair < unsigned int, float > > > ().swap(links_in);

    // outputs correspondance
    if (do_renumber) {
        for (auto c : this->correspondance) {
            fbin.write((char * )( & c), sizeof(ULLI));
        }
    }

    fbin.close();
}
