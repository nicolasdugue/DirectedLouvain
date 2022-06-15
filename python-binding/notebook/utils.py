import urllib
import networkx as nx
from os import system
import subprocess
import networkit as nk
import numpy as np

from sklearn.metrics import normalized_mutual_info_score, adjusted_mutual_info_score
import networkit as nk

import contextlib

#C++ binding of our code
import directedlouvain as dl



#### IO - loading a graph from an external URL
def load_graph(url, nom_graphe):
    sock = urllib.request.urlopen(url)  # open URL
    s = sock.read()  # read into BytesIO "file"
    sock.close()
    gml = s.decode()  # read gml data
    G = nx.parse_gml(gml)  # parse gml data
    G=nx.relabel_nodes(G, {val:idx for idx,val in enumerate(G.nodes())})
    G.remove_edges_from(nx.selfloop_edges(G))
    #print("directed ?",G.is_directed())
    nx.write_edgelist(G,nom_graphe,data=False)
    return G

#### IO - loading a graph from a local path
def load_graph_local(url, nom_graphe):
    G=nx.read_gml(url)
    #For citeseer, some nodes don't have any membership, I remove them
    nodes=list(G.nodes()).copy()
    for u in nodes:
        try:
            G.nodes[u]["gt"]
        except:
            G.remove_node(u)
    G=nx.relabel_nodes(G, {val:idx for idx,val in enumerate(G.nodes())})
    G.remove_edges_from(nx.selfloop_edges(G))
    #print("directed ?",G.is_directed())
    nx.write_edgelist(G,nom_graphe,data=False)
    return G
    

def get_membership(G, kw_membership="gt"):
    return [G.nodes[u][kw_membership] for u in G.nodes]


# Running directed Louvain algorithm using our Python binding
def run_directed_louvain(nom_graphe, nom_hierarchy, display_level=-2, weighted=False):
    dl_obj=dl.Community(filename=nom_graphe, weighted=weighted)
    output=open(nom_hierarchy, "w")
    with contextlib.redirect_stdout(output):
        if display_level != -2:
            maximum=dl_obj.run(verbose=False, display_level=-1)
        else:
            maximum=dl_obj.run(verbose=False, display_level=display_level)
    output.close()
    fichier=open(nom_hierarchy)
    lignes = fichier.readlines()
    membership = np.zeros(len(lignes))
    if display_level == -2:
        for ligne in lignes:
                tab = ligne.split(" ")
                membership[int(tab[0])]=int(tab[1].strip())
        return membership
    else:
        return get_partition(dl_obj,nom_hierarchy, maximum, level=display_level)
            
    
# Used when one do not want to extract the highest level of the Louvain algorithm
# It is for example useful when running the ensemble clustering for graphs algorithm
def get_partition(dl_obj, nom_hierarchy, maximum, level=-1 ):
    nom_level=nom_hierarchy+"_lvl"+str(level)
    output=open(nom_level, "w")
    with contextlib.redirect_stdout(output):
        if level == -1:
            dl_obj.print_level(maximum -1)
        else:
            dl_obj.print_level(level)
    output.close()
    fichier=open(nom_level)
    lignes = fichier.readlines()
    membership = np.zeros(len(lignes))
    for ligne in lignes:
            tab = ligne.split(" ")
            membership[int(tab[0])]=int(tab[1].strip())
    return membership

def experiment_louvain(url, load_graph=load_graph, weighted=False, nb_runs=20):
    nom_graphe, nom_hierarchy=get_names(url)
    G=load_graph(url, nom_graphe)
    gt=get_membership(G)
    
    nmi=[]
    for i in range(nb_runs):
        membership=run_directed_louvain(nom_graphe, nom_hierarchy)
        nmi.append(normalized_mutual_info_score(gt, membership))
    print(np.mean(nmi), np.std(nmi))

    
def count_com(labels):
    return len(set(labels))

def get_dico(liste):
    dico=dict()
    for (name, measure) in liste:
        dico[name]=[]
    return dico

def get_metrics():
    return [("nmi", normalized_mutual_info_score), ("ami", adjusted_mutual_info_score), ("count_com", count_com)]
    
##
def compare_louvain(url, load_graph=load_graph, weighted=False, nb_runs=50):
    nom_graphe, nom_hierarchy=get_names(url)
    G=load_graph(url, nom_graphe)
    gt=get_membership(G)  
    metrics = get_metrics()
    dir_ = get_dico(metrics)
    undir_=get_dico(metrics)
    
    undirected = G.to_undirected()
    G_nk= nk.nxadapter.nx2nk(undirected, weightAttr=None)
    for i in range(nb_runs):
        communities = nk.community.detectCommunities(G_nk, algo=nk.community.PLM(G_nk, True), inspect=False)
        membership=run_directed_louvain(nom_graphe, nom_hierarchy)
        #map_eq=nk.community.detectCommunities(G_nk, algo=nk.community.LouvainMapEquation(G_nk, True), inspect=False)
        for (name, m) in metrics:
            if "count" in name:
                dir_[name].append(m(membership))
            else:
                dir_[name].append(m(membership, gt))
        
        for (name, m) in metrics:
            if "count" in name:
                undir_[name].append(m(communities.getVector()))
            else:
                undir_[name].append(m(communities.getVector(), gt))
                
    for name in dir_:
        print("Louvain, ", name, " : ", np.mean(dir_[name]), np.std(dir_[name]), "directed : True")
    for name in undir_:
        print("Louvain, ", name, " : ", np.mean(undir_[name]), np.std(undir_[name]), "directed : False")
    #print("Map Equation",np.mean(nmi_map), np.std(nmi_map))

def get_names(url):
    tab=url.split("/")
    nom_xp=tab[-1].split(".")[0]
    nom_graphe=nom_xp+".edge"
    nom_hierarchy=nom_xp+".tree"
    return nom_graphe,nom_hierarchy



##### ECG method comes from
#V. Poulin and F. Théberge, Ensemble clustering for graphs: comparisons and applications, Network Science (2019) 4:51 https://doi.org/10.1007/s41109-019-0162-z or #https://rdcu.be/bLn9i
# It is adapted to directed graph with our directed Louvain algorithm, but please, take a look at https://github.com/ftheberge/Ensemble-Clustering-for-Graphs for the original code

from networkx.algorithms.core import core_number
from collections import namedtuple

#Ensemble Clustering for Graphs
def ecg(G, nom_graphe, level=-1, ens_size = 16, min_weight = 0.05, directed=True):
    W = {k:0 for k in G.edges()}
    ## Ensemble of level-1 Louvain 
    nom_hierarchy = "tmp_hierarchy"
    if not directed:
        G_nk= nk.readGraph(nom_graphe, nk.Format.EdgeListSpaceZero)
        
    for i in range(ens_size):
        if directed:
            l = run_directed_louvain(nom_graphe, nom_hierarchy, display_level=level)
        else:
            l = nk.community.detectCommunities(G_nk, algo=nk.community.PLM(G_nk), inspect=False).getVector()
        #l = get_partition(nom_hierarchy, level=1)
        for e in G.edges():
            if ((e[1], e[0]) not in G.edges()):
                    W[e] += int(l[e[0]] == l[e[1]]) *(1 + min_weight)
            W[e] += int(l[e[0]] == l[e[1]])
    ## vertex core numbers
    core = core_number(G)
    ## set edge weights
    for e in G.edges():
        m = min(core[e[0]],core[e[1]])
        if m > 1:
            W[e] = min_weight + (1-min_weight)*W[e]/ens_size
        else:
            W[e] = min_weight

    fichier=open(nom_graphe+"_ecg", "w")
    for u,v in W:
        fichier.write(str(u)+" "+str(v)+" "+ str(W[(u,v)])+"\n")
    fichier.close()
    
    G_nk= nk.readGraph(nom_graphe+"_ecg", nk.Format.EdgeListSpaceZero, weighted=True)
    communities = nk.community.detectCommunities(G_nk, algo=nk.community.PLM(G_nk, True, recurse=True), inspect=False).getVector()
    return communities




def experiment_ecg(url, load_graph=load_graph,level=-1,ens_size=64, nb_run=50, directed=True):
    nom_graphe, nom_hierarchy=get_names(url)
    G=load_graph(url, nom_graphe)
    gt=get_membership(G)
    
    metrics= get_metrics()
    
    values=get_dico(metrics)
    for i in range(nb_run):
        G=load_graph(url, nom_graphe)
        labels = ecg(G, nom_graphe, ens_size=ens_size, level=level, directed=directed)
        for (name, m) in metrics:
            if "count" in name:
                values[name].append(m(labels))
            else:
                values[name].append(m(labels, gt))
    for name in values:
        print("ecg, ", name, " : ", np.mean(values[name]), np.std(values[name]), "directed : ", directed)


    