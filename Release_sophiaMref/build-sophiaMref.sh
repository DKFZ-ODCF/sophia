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

$CPP $CPP_OPTS -o "GlobalAppConfig.o" "../src/GlobalAppConfig.cpp"
$CPP $CPP_OPTS -o "ChrConverter.o" "../src/ChrConverter.cpp"
$CPP $CPP_OPTS -o "Hg37ChrConverter.o" "../src/Hg37ChrConverter.cpp"
$CPP $CPP_OPTS -o "Hg38ChrConverter.o" "../src/Hg38ChrConverter.cpp"
$CPP $CPP_OPTS -o "HelperFunctions.o" "../src/HelperFunctions.cpp"
$CPP $CPP_OPTS -o "sophiaMref.o" "../src/sophiaMref.cpp"
$CPP $CPP_OPTS -o "SuppAlignment.o" "../src/SuppAlignment.cpp"
$CPP $CPP_OPTS -o "SuppAlignmentAnno.o" "../src/SuppAlignmentAnno.cpp"
$CPP $CPP_OPTS -o "MrefEntry.o" "../src/MrefEntry.cpp"
$CPP $CPP_OPTS -o "MrefEntryAnno.o" "../src/MrefEntryAnno.cpp"
$CPP $CPP_OPTS -o "MrefMatch.o" "../src/MrefMatch.cpp"
$CPP $CPP_OPTS -o "MasterRefProcessor.o" "../src/MasterRefProcessor.cpp"
$CPP $CPP_OPTS -o "Breakpoint.o" "../src/Breakpoint.cpp"
$CPP $CPP_OPTS -o "BreakpointReduced.o" "../src/BreakpointReduced.cpp"
$CPP $CPP_OPTS -o "GermlineMatch.o" "../src/GermlineMatch.cpp"
$CPP $CPP_OPTS -o "DeFuzzier.o" "../src/DeFuzzier.cpp"

$CPP -L$CONDA_PREFIX/lib -o "sophiaMref" \
    GlobalAppConfig.o \
    ChrConverter.o \
    Hg37ChrConverter.o \
    Hg38ChrConverter.o \
    HelperFunctions.o \
    SuppAlignment.o \
    SuppAlignmentAnno.o \
    MrefEntry.o \
    MrefEntryAnno.o \
    MrefMatch.o \
    MasterRefProcessor.o \
    Breakpoint.o \
    BreakpointReduced.o \
    GermlineMatch.o \
    DeFuzzier.o \
    sophiaMref.o \
    $LD_OPTS -lz -lboost_system -lboost_iostreams -lboost_program_options
