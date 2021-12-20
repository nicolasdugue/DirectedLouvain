// File: graph.h
// -- simple graph handling header file
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

#ifndef GRAPH_HPP
#define GRAPH_HPP

#define WEIGHTED   0
#define UNWEIGHTED 1

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>

typedef unsigned long long int ULLI;

using namespace std;

class Graph {
    private:
        /* TODO: ULLI if not renumbered ? */
        vector<vector<pair<unsigned int,float> > > links_out;
        vector<vector<pair<unsigned int,float> > > links_in;
        /* FIXME: is it possible to have a map with adjustable value types ?
         * vector<void*> and use necessary memory? */
        vector<ULLI> correspondance;

    public:
        Graph (string in_filename, string filename, string filename_w, int type, bool do_renumber);
        Graph (const Graph& );
        ~Graph() {}

        /* FIXME: is it used somewhere? */
        void clean(int type);
        void display(int type);
        void display_binary(string filename, string filename_w, int type, bool do_renumber);
};

#endif // GRAPH_H
