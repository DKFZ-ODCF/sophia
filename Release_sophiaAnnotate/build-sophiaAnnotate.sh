#!/bin/bash

set -uex
trap 'echo "Compilation failed with an error" >> /dev/stderr' ERR

install_strtk() {
    wget -c https://github.com/ArashPartow/strtk/raw/master/strtk.hpp -O ../include/strtk.hpp
    # sed -ir 's/#include <boost\/lexical_cast.hpp>/#include <boost\/convert\/lexical_cast.hpp>/' ../include/strtk.hpp
}

install_strtk

CPP=x86_64-conda_cos6-linux-gnu-g++
CONDA_ENV_ROOT=/data/kensche/work/share/miniconda3/envs/tools
INCLUDES="-I../include -I$CONDA_ENV_ROOT/include"
CPP_OPTS="-L$CONDA_ENV_ROOT/lib -std=c++1z $INCLUDES -O3 -Wall -Wextra -static -static-libgcc -static-libstdc++ -flto -c -fmessage-length=0 -Wno-attributes"

if [[ "${STATIC:-false}" == "true" ]]; then
    CPP_OPTS="-static -static-libgcc -static-libstdc++ $CPP_OPTS"
fi

$CPP $CPP_OPTS -o "AnnotationProcessor.o" "../src/AnnotationProcessor.cpp"
$CPP $CPP_OPTS -o "Breakpoint.o" "../src/Breakpoint.cpp"
$CPP $CPP_OPTS -o "BreakpointReduced.o" "../src/BreakpointReduced.cpp"
$CPP $CPP_OPTS -o "ChrConverter.o" "../src/ChrConverter.cpp"
$CPP $CPP_OPTS -o "DeFuzzier.o" "../src/DeFuzzier.cpp"
$CPP $CPP_OPTS -o "GermlineMatch.o" "../src/GermlineMatch.cpp"
$CPP $CPP_OPTS -o "MrefEntry.o" "../src/MrefEntry.cpp"
$CPP $CPP_OPTS -o "MrefEntryAnno.o" "../src/MrefEntryAnno.cpp"
$CPP $CPP_OPTS -o "MrefMatch.o" "../src/MrefMatch.cpp"
$CPP $CPP_OPTS -o "SuppAlignment.o" "../src/SuppAlignment.cpp"
$CPP $CPP_OPTS -o "SuppAlignmentAnno.o" "../src/SuppAlignmentAnno.cpp"
$CPP $CPP_OPTS -o "SvEvent.o" "../src/SvEvent.cpp"
$CPP $CPP_OPTS -o "sophiaAnnotate.o" "../sophiaAnnotate.cpp"

$CPP -L$CONDA_ENV_ROOT/lib -flto -o "sophiaAnnotate"  AnnotationProcessor.o Breakpoint.o BreakpointReduced.o ChrConverter.o DeFuzzier.o GermlineMatch.o MrefEntry.o MrefEntryAnno.o MrefMatch.o SuppAlignment.o SuppAlignmentAnno.o SvEvent.o  sophiaAnnotate.o   -lz -lboost_system -lboost_iostreams
