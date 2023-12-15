#!/bin/bash

set -uex
trap 'echo "Compilation failed with an error" >> /dev/stderr' ERR

CONDA_PREFIX="${CONDA_PREFIX:?No CONDA_PREFIX -- no active Conda environment}"

install_strtk() {
    wget -c https://github.com/ArashPartow/strtk/raw/master/strtk.hpp -O ../include/strtk.hpp
}

install_strtk

CPP=x86_64-conda_cos6-linux-gnu-g++
INCLUDES="-I../include -I$CONDA_PREFIX/include"

CPP_OPTS="-L$CONDA_PREFIX/lib -std=c++17 $INCLUDES -O3 -Wall -Wextra -flto -c -fmessage-length=0 -Wno-attributes"

LD_OPTS=""
if [[ "${STATIC:-false}" == "true" ]]; then
    LD_OPTS="-static -static-libgcc -static-libstdc++"
fi

$CPP $CPP_OPTS -o "Alignment.o" "../src/Alignment.cpp"
$CPP $CPP_OPTS -o "Breakpoint.o" "../src/Breakpoint.cpp"
$CPP $CPP_OPTS -o "ChosenBp.o" "../src/ChosenBp.cpp"
$CPP $CPP_OPTS -o "GlobalAppConfig.o" "../src/GlobalAppConfig.cpp"
$CPP $CPP_OPTS -o "ChrConverter.o" "../src/ChrConverter.cpp"
$CPP $CPP_OPTS -o "Hg37ChrConverter.o" "../src/Hg37ChrConverter.cpp"
$CPP $CPP_OPTS -o "Hg38ChrConverter.o" "../src/Hg38ChrConverter.cpp"
$CPP $CPP_OPTS -o "SamSegmentMapper.o" "../src/SamSegmentMapper.cpp"
$CPP $CPP_OPTS -o "Sdust.o" "../src/Sdust.cpp"
$CPP $CPP_OPTS -o "SuppAlignment.o" "../src/SuppAlignment.cpp"
$CPP $CPP_OPTS -o "HelperFunctions.o" "../src/HelperFunctions.cpp"
$CPP $CPP_OPTS -o "sophia.o" "../src/sophia.cpp"

$CPP -L$CONDA_PREFIX/lib -o "sophia"  \
    Alignment.o \
    Breakpoint.o \
    ChosenBp.o \
    ChrConverter.o \
    Hg37ChrConverter.o \
    Hg38ChrConverter.o \
    SamSegmentMapper.o \
    Sdust.o \
    SuppAlignment.o \
    HelperFunctions.o \
    GlobalAppConfig.o \
    sophia.o \
    $LD_OPTS -lboost_program_options
