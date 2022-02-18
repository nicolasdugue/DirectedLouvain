# [Python binding](#pythonbinding)                       

---
## Add pybind11 submodule

`pybind11` is included as a submodule of this project. Thus, in order to bind our `C++` implementation with `python3`, first do the following:

    git submodule init
    git submodule update

### Computing communities

The graph **must** be in edgelist format, that is one edge per line as follows (the `weight` being optional):  

    src dest [weight]

Moreover, it is **mandatory** that vertices of the input graph are numbered from `0` to `n-1`. 
To ensure a proper computation of the communities, the default computation encompasses a renumbering of the input graph. 
The option `-n` indicates that the graph is already numbered from `0` to `n-1` and hence avoids renumbering. 
**Important**: communities are written using the **original label nodes**.

The standard command is:

    ./bin/community -f graph/graph.txt -l -1 -v > graph.tree

Another possibility is to pass a binary file containing all information regarding the graph. 
This file **must** be generated using the our program in a first place, or follows the CSR format (see below). 
In this case, the only mandatory option is `-w` to indicate whether the graph is weighted. 

Several options are available, among which:
+ `-f` path to the input graph (edgelist or binary format (`.bin`))
+ `-w` to indicate that the input graph is weighted
+ `-n` to indicate that the input graph is correctly numbered (from `0` to `n-1`)
+ `-r` for reproducibility purposes: the renumbered graph is stored on hard drive. 

More options and information are provided using `./bin/community`

### [Graph representation: CSR format](#CSR)

Graphs are stored under the [Compressed Sparse Row (CSR)](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)) format.  
Several structures are containing the whole graph information: 
+ two arrays of **cumulative degrees** (out- and in-degrees): for out-degrees, each index _i_ contains the sum of all out-degrees for nodes from _0_ to _i_. 
+ two arrays of outcoming and incoming **arcs**: the _d(0)_ (out- or in-degree of node _0_) first values contain neighbors of node _0_ and so on.   
To find the first neighbor of a given node _i_ one simply needs to consider the difference between cumulative degrees of _i_ and _i-1_.
+ two array of outcoming and incoming **weights**: similar to the previous ones but store weights instead of node identifiers. 

Example of CSR format for a directed graph. The displayed arrays contain information regarding out-neighbors and weighted out--degrees only.
![CSR example](docs/CSR.png "Example of CSR format for a directed graph. The displayed arrays contain information regarding out-neighbors and weighted out--degrees only.")

## Examples 
Using `graph/graph.txt` one obtains: 

    ./bin/community -f graph/graph.txt -l -1 -v -w -r > graph/graph.tree

to compute hierarchical community structure (using **original** label nodes) 
by first renumbering the graph, and 
then writing files for reproducibility. The next runs would thus be: 

    ./bin/community -f graph/graph.bin -w -l -1 > graph/graph.tree

Finally, using an already renumbered graph one gets: 

    ./bin/community -f graph/graph_renum.txt -l -1 -v -w -n > graph/graph.tree

The program can also start with any given partition using -p option

    ./community graph.bin -p graph.part -v

### Improvements

To ensure a faster computation (with a loss of quality), one can use
the -q option to specify that the program must stop if the increase of
modularity is below epsilon for a given iteration or pass:

    ./bin/community graph/graph.bin -w -l -1 -q 0.0001 > graph/graph.tree

-----------------------------------------------------------------------------
**Display communities information**

Displays information on the tree structure (number of hierarchical
levels and nodes per level):

    ./hierarchy graph.tree

Displays the belonging of nodes to communities for a given level of
the tree:

    ./hierarchy graph.tree -l 2 > graph_node2comm_level2

-----------------------------------------------------------------------------
## References
* **[1]** Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre. [Fast unfolding of communities in large networks](https://arxiv.org/pdf/0803.0476.pdf). Journal of Statistical Mechanics: Theory and Experiment, 2008, vol. 2008, no 10, p. P10008.
* **[2]** Alexandre Arenas, Jordi Duch, Alberto Fern\'andez, Sergio G\'omez, 2007. [Size reduction of complex networks preserving modularity](https://iopscience.iop.org/article/10.1088/1367-2630/9/6/176/pdf). New Journal of Physics 9, 176.
* **[3]** Nicolas Dugué, Anthony Perez. [Directed Louvain: maximizing modularity in directed networks](https://hal.archives-ouvertes.fr/hal-01231784/document). [Research Report] Université d'Orléans. 2015.