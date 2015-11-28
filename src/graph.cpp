// File: graph.cpp
// -- simple graph handling source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include<cstring>
#include "../include/graph.h"

using namespace std;

Graph::Graph(char *in_filename, char *filename, char* filename_w, int type, bool do_renumber) {
  ifstream finput;
  
  if(do_renumber) {
	map<unsigned long long int, unsigned int>::iterator it_src;
	map<unsigned long long int, unsigned int>::iterator it_dest;	
  
  unsigned int cpt = 1;

        
  unsigned int nb_links=0;
  float weight = 1.f;
  ofstream foutput;
	char* tmp = new char[strlen(in_filename) + 6];
	strcat(tmp, in_filename);
	strcat(tmp, "_renum");
  foutput.open(tmp, fstream::out | fstream::binary);

		cout << "Renumerotation begins..." << endl;
		
		correspondance.resize(0);

        /* Creating a bunch of maps... NOT ANYMORE :( */
	unsigned int* corres = new unsigned int[numeric_limits<unsigned int>::max()];
	memset(corres, 0, numeric_limits<unsigned int>::max()*sizeof(unsigned int));
	map<unsigned long long int , unsigned int > corres_big_ids;

        finput.open(in_filename,fstream::in);
	unsigned int src_prec, dest_prec;
	src_prec = 1;
	dest_prec = 1;

        if(finput) {

	unsigned long long int src, dest;
	unsigned int pos_src, pos_dest;

	/* We first do the renumerotation and build the correspondance */
        while(finput >> src >> dest) {

                        if(type == WEIGHTED)
                                finput >> weight;

			if(src < numeric_limits<unsigned int>::max()) {
				if(corres[src] == 0) {

                                        corres[src] = cpt;

                                correspondance.resize(correspondance.size()+1);
                                correspondance[cpt-1] = src;
                                cpt++;
			}}

			else {

				if(corres_big_ids.find(src) == corres_big_ids.end()) {

				corres_big_ids.insert(make_pair(src, cpt));

				correspondance.resize(correspondance.size()+1);
                                correspondance[cpt-1] = src;
                                cpt++;

			}}

			if(dest < numeric_limits<unsigned int>::max()) {
                                if(corres[dest] == 0) {

                                        corres[dest] = cpt;

                                correspondance.resize(correspondance.size()+1);
                                correspondance[cpt-1] = dest;
                                cpt++;
                        }}

                        else {

                                if(corres_big_ids.find(dest) == corres_big_ids.end()) {

                                corres_big_ids.insert(make_pair(src, cpt));

                                correspondance.resize(correspondance.size()+1);
                                correspondance[cpt-1] = dest;
                                cpt++;

                        }}

			unsigned int pos_src, pos_dest;
                                           if(src < numeric_limits<unsigned int>::max()) {
                                                pos_src = corres[src] - 1;
                                          }
                                           else {
                                                map<unsigned long long int, unsigned int>::iterator it_src;
                                                it_src = corres_big_ids.find(src);
                                                pos_src = it_src->second - 1;
                                        }

                                           if(dest < numeric_limits<unsigned int>::max()) {
                                                pos_dest = corres[dest] - 1;
                                           } else{
                                                map<unsigned long long int, unsigned int>::iterator it_dest;
                                                it_dest = corres_big_ids.find(dest);
                                                pos_dest = it_dest->second - 1;
                                        }


					      nb_links++;
	
		// TODO: count nb_links while parsing for the maximum, and pourcent the progression 			
		if(nb_links % 500 == 0) cout << "50000000 ecrits" << endl;
						foutput << pos_src << " " << pos_dest << endl;

                     	src_prec = src;
                        dest_prec = dest;

                }

	finput.close();
	cout << "Renumerotation ends..." << endl;
	cout << "Building the graph... links out" << endl; 
	foutput.close();

	/* Release memory */
	delete[] corres;

	src_prec = 1;
	dest_prec = 1;

	// Then we build the graph reading the new graph
	// Out links first
	finput.open(tmp, fstream::in);
	while(finput >> src >> dest) {

                        if(type == WEIGHTED)
                                finput >> weight;

		if (links_out.size()<=max(src,dest)+1) {
                                                links_out.resize(max(src,dest)+1);
                                              }
                                              
                                              links_out[src].push_back(make_pair(dest,weight));

                if (links_in.size()<=max(src,dest)+1) {
                                                links_in.resize(max(src,dest)+1);
                                              }
  
                                              links_in[dest].push_back(make_pair(src,weight));
                                              
                                              nb_links++;

                // TODO: count nb_links while parsing for the maximum, and pourcent the progression
                if(nb_links % 50000000 == 0) cout << "50000000 ecrits" << endl;

	}

	cout << "done." << endl;

	finput.close(); 
	
  
	unsigned int s = links_out.size();					
  // outputs weights in a separate file
  if (type==WEIGHTED) {
    ofstream foutput_w;
    foutput_w.open(filename_w,fstream::out | fstream::binary);
    for (unsigned int i=0 ; i<s ; i++) {
      for (unsigned int j=0 ; j<links_out[i].size() ; j++) {
        float weight = links_out[i][j].second;
        foutput_w.write((char *)(&weight),sizeof(float));
      }
    }
    foutput_w.close();
  }
  cout << "done" << endl;

}

}
else {

int src, dest, cpt;
  cpt = 0;

        int src_prec, dest_prec;
        src_prec = -1;
        dest_prec = -1;
        
  int nb_links=0;
  float weight = 1.f;

                
        finput.open(in_filename, fstream::in);
	while(finput >> src >> dest) {

                        if(type == WEIGHTED)
                                finput >> weight;

		if (links_out.size()<=max(src,dest)+1) {
                                                links_out.resize(max(src,dest)+1);
                                              }
                                              
                                              links_out[src].push_back(make_pair(dest,weight));

                if (links_in.size()<=max(src,dest)+1) {
                                                links_in.resize(max(src,dest)+1);
                                              }
  
                                              links_in[dest].push_back(make_pair(src,weight));
                                              
                                              nb_links++;

                // TODO: count nb_links while parsing for the maximum, and pourcent the progression
                if(nb_links % 50000000 == 0) cout << "50000000 ecrits" << endl;

	}

	cout << "done." << endl;

	finput.close(); 
	
}


	ofstream fbin;
  fbin.open(filename, fstream::out | fstream::binary);
// Writing out information
	unsigned int s = links_out.size();
	cout << "number of nodes : " << s << endl;
	cout << "writing in binary file..." << endl;
	// outputs number of nodes
  fbin.write((char *)(&s),sizeof(int));

  // outputs cumulative degree sequence
  /* Contient uniquement les degres sortants en oriente */
  long tot=0;
  for (unsigned int i=0 ; i<s ; i++) {
      tot+=(long)links_out[i].size();
    fbin.write((char *)(&tot),sizeof(long));
  }


	// outputs links_out
  for (unsigned int i=0 ; i<s ; i++) {
    for (unsigned int j=0 ; j<links_out[i].size() ; j++) {
      unsigned long long int dest = links_out[i][j].first;
      fbin.write((char *)(&dest),sizeof(unsigned int));
    }
  }
	cout << "done." << endl;
	cout << "releasing memory..." << endl;
	// Release memory
	for(unsigned int i = 0; i < links_out.size(); i++) {

		links_out[i].clear();
		vector<pair<unsigned int, float> >().swap( links_out[i] );

	}

	links_out.clear();
	vector<vector<pair<unsigned int,float> > >().swap( links_out );
	
	// Writing information
  /* Contient uniquement les degres entrants en oriente */
 long tot_in=0;
  for (unsigned int i=0 ; i<s ; i++) {
        tot_in+=(long)links_in[i].size();
    fbin.write((char *)(&tot_in),sizeof(long));
  }

	// outputs links_in
  for (unsigned int i=0 ; i<s ; i++) {
    for (unsigned int j=0 ; j<links_in[i].size() ; j++) {
      unsigned long long int dest = links_in[i][j].first;
        fbin.write((char *)(&dest),sizeof(unsigned int));
    }
  }
	cout << "releasing memory..." << endl;
	// Release memory
	for(unsigned int i = 0; i < links_in.size(); i++) {

		links_in[i].clear();
		vector<pair<unsigned int, float> >().swap( links_in[i] );

	}

	links_in.clear();
	vector<vector<pair<unsigned int,float> > >().swap( links_in );

	 // outputs correspondance
  if(do_renumber) {
  for(unsigned int i = 0; i <s ; i++) {

                unsigned long long int corr = correspondance[i];
                fbin.write((char *)(&corr),sizeof(unsigned long long int));

  } }

  fbin.close();
  finput.close();

}

/* Permet de gerer les multigraphes : si on a dans le fichier 1 4 3 et 1 4 6, 
 * le graphe original aura deux paires (4,3) et (4,6) dans la liste d'adjacence 
 * de 1. Du coup, cette procedure remplacera ca par (4,9). 
 */
void
Graph::clean(int type) {
  for (unsigned int i=0 ; i<links_out.size() ; i++) {
    map<unsigned long long int, float> m;
    map<unsigned long long int, float>::iterator it;

    for (unsigned int j=0 ; j<links_out[i].size() ; j++) {
      it = m.find(links_out[i][j].first);
      if (it==m.end())
	m.insert(make_pair(links_out[i][j].first, links_out[i][j].second));
      else if (type==WEIGHTED)
      	it->second+=links_out[i][j].second;
    }
    
    vector<pair<unsigned int,float> > v;
    for (it = m.begin() ; it!=m.end() ; it++)
      v.push_back(*it);
    links_out[i].clear();
    links_out[i]=v;
  }
}

void
Graph::display(int type) {
  for (unsigned int i=0 ; i<links_out.size() ; i++) {
    for (unsigned int j=0 ; j<links_out[i].size() ; j++) {
      int dest   = links_out[i][j].first;
      float weight = links_out[i][j].second;
      if (type==WEIGHTED)
	cout << i << " " << dest << " " << weight << endl;
      else
	cout << i << " " << dest << endl;
    }
  }
}

void
Graph::display_binary(char *filename, char *filename_w, int type, bool do_renumber) {
  ofstream foutput;
  foutput.open(filename, fstream::out | fstream::binary);

  unsigned int s = links_out.size();

  // outputs number of nodes
  foutput.write((char *)(&s),sizeof(int));

  // outputs cumulative degree sequence
  /* Contient uniquement les degres sortants en oriente */
  long tot=0;
  for (unsigned int i=0 ; i<s ; i++) {
      tot+=(long)links_out[i].size();
    foutput.write((char *)(&tot),sizeof(long));
  }
  cout << tot << endl;

  // outputs cumulative degree sequence
  /* Contient uniquement les degres entrants en oriente */
  long tot_in=0;
  for (unsigned int i=0 ; i<s ; i++) {
    	tot_in+=(long)links_in[i].size();
    foutput.write((char *)(&tot_in),sizeof(long));
  }
cout << tot_in << endl;
  
  // outputs correspondance 
  if(do_renumber) {
  for(unsigned int i = 0; i <s ; i++) {
  	
  		unsigned long long int corr = correspondance[i];
		foutput.write((char *)(&corr),sizeof(unsigned long long int));
  	
  } 
}

  // outputs links_out
  for (unsigned int i=0 ; i<s ; i++) {
    for (unsigned int j=0 ; j<links_out[i].size() ; j++) {
      unsigned long long int dest = links_out[i][j].first;
      foutput.write((char *)(&dest),sizeof(unsigned int));
    }
  }

  // outputs links_in
  for (unsigned int i=0 ; i<s ; i++) {
    for (unsigned int j=0 ; j<links_in[i].size() ; j++) {
      unsigned long long int dest = links_in[i][j].first;
	foutput.write((char *)(&dest),sizeof(unsigned int));
    }
  }
  foutput.close();
  
  // outputs weights in a separate file
  if (type==WEIGHTED) {
    ofstream foutput_w;
    foutput_w.open(filename_w,fstream::out | fstream::binary);
    for (unsigned int i=0 ; i<s ; i++) {
      for (unsigned int j=0 ; j<links_out[i].size() ; j++) {
	float weight = links_out[i][j].second;
	foutput_w.write((char *)(&weight),sizeof(float));
      }
    }
    foutput_w.close();
  }
  cout << "done" << endl;
}

