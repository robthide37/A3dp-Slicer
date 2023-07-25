cmake_minimum_required(VERSION 2.8.12)
project(heatshrink)

add_library(heatshrink STATIC
    heatshrink_common.h
    heatshrink_config.h
    heatshrink_decoder.c
    heatshrink_decoder.h
    heatshrink_encoder.c
    heatshrink_encoder.h
)
