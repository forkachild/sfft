cmake_minimum_required(VERSION 3.25)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Werror")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(sfft C)

add_library(sfft
    fft.c)

target_include_directories(sfft PUBLIC .)