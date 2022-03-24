CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)
PROJECT(polyclipping)


ADD_LIBRARY(polyclipping cpp/clipper.cpp)
target_include_directories(polyclipping PUBLIC cpp)
