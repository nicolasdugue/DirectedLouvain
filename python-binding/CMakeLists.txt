cmake_minimum_required(VERSION 3.1.0)
project (directedlouvain)
add_subdirectory(../extern/pybind11 REQUIRED)
pybind11_add_module(directedlouvain pybind_main.cpp community.cpp graph.cpp)
