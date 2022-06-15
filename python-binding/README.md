# [Python binding](#pythonbinding)                       

---
## Add [pybind11](https://github.com/pybind/pybind11) submodule

`pybind11` is included as a submodule of this project. Thus, in order to bind our `C++` implementation with `python3`, first do the following:

    git submodule init
    git submodule update

This is the only step needed to be able to use [`pybind11`](https://github.com/pybind/pybind11).

### Create `python` library 

The next step requires [`cmake`](https://cmake.org/) and can be easily reproduced using the following steps:

    cd build
    cmake ..
    make

That's all folks! After these steps you should see a `.so` file added to the `build` directory. Now any `python` program 
**started from this directory** can import `directedlouvain`. To work from another directory you can simply copy the `.so` file into it. The library is minimalist on purpose, with the following methods: 

+ `init`: creates a `Community` object from a graph (edgelist or binary)
+ `run`: computes communities using Directed Louvain algorithm and returns the number of levels computed
+ `modularity`: computes he modularity of the last hierarchical level
+ `print_level`: displays the community structure at some given level (between 0 and the value returned by run()-1, which is the last level of the structure)

More information about these methods (and their arguments) can be found inside the `python-use.py` file. 
