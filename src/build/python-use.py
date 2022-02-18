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
    dl = directedlouvain.init(sys.argv[1])
    # or for instance
    # dl = directedlouvain.create(sys.argv[1], weighted=True, renumbering=False)
    '''this method uses three arguments:
            + [verbose] boolean value indicating whether to print information on stderr (default: false)
            + [display_level] int value indicating which hierarchical level to display on stdout (default: -1, i.e. all levels)
            + [filename_part] string value indicating the path to a partition file for initialization (default: "")
    '''
    # indicating an output file for all stdout output
    dl.run(verbose=True)
    print("modularity of the last computed hierarchical level: {:.6f}".format(dl.modularity()),file=sys.stderr)