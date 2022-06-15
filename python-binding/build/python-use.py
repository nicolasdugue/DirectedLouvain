# Module must first be binded using pybind11
# Usage: python3.8 python-use.py path/to/graph > path/to/hierarchy_file
# See README: 

import directedlouvain
import sys

if __name__=="__main__":
    '''this method uses five arguments:
            + [filename] string value indicating the path to the graph (edgelist or binary format)
            + [weighted] boolean value indicating whether the graph is weighted (default: false)
            + [precision] double value indicating the required precision for modularity improvement (default: 0.0001)
            + [reproducibility] boolean value indicating whether to write renumbered graph on hard drive (default: false)
            + [renumbering] boolean value indicating whether to renumber the input grah (default: true)
    '''
    dl = directedlouvain.Community(sys.argv[1], weighted=True)
    # or for instance for an unweighted graph with nodes not ranging from 0 to N-1
    # dl = directedlouvain.create(sys.argv[1], renumbering=False)
    '''this method uses three arguments:
            + [verbose] boolean value indicating whether to print information on stderr (default: false)
            + [display_level] int value indicating which hierarchical level to display on stdout (default: -1, i.e. hierarchical community structure)
            + [filename_part] string value indicating the path to a partition file for initialization (default: "")
        and returns the number of levels of the hierarchical communit structure
    '''
    levels = dl.run(verbose=True)
    print("modularity of the last computed hierarchical level: {:.6f}".format(dl.modularity()),file=sys.stderr)
    print("last level (between 0 and levels-1)",file=sys.stderr)
    dl.print_level(levels-1)
