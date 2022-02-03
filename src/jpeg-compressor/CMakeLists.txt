# PrusaSlicer specific CMake

cmake_minimum_required(VERSION 2.8.12)
project(jpeg-compressor)

add_library(jpeg-compressor STATIC
    jpge.h
    jpge.cpp
)
