#!/bin/bash
# test cases for example-0404 (1D problem with Hodgkin-Huxley)

echo "compiling and running example $(pwd)"

folder=$1

mkdir -p $folder

echo "  compiling $folder"
cd $folder
cmake -DCMAKE_BUILD_TYPE=$folder -DOPENCMISS_BUILD_TYPE=$folder ..
make
cd ..
echo "  running $folder"
# <number elements X> <interpolation type> <solver type> <PDE step size> <stop time> <output frequency> <CellML Model URL> <slow-twitch> <ODE time-step>

mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt005 && ./$folder/src/example 2047 1 0 0.005 10 10 slow_TK_2014_12_08.xml T 0.001 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt005
mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt004 && ./$folder/src/example 2047 1 0 0.004 10 10 slow_TK_2014_12_08.xml T 0.0008 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt004
mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt002 && ./$folder/src/example 2047 1 0 0.002 10 25 slow_TK_2014_12_08.xml T 0.0004 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt002
mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt001 && ./$folder/src/example 2047 1 0 0.001 10 50 slow_TK_2014_12_08.xml T 0.0002 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt001
mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt0005 && ./$folder/src/example 2047 1 0 0.0005 10 100 slow_TK_2014_12_08.xml T 0.0001 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt0005
mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt0004 && ./$folder/src/example 2047 1 0 0.0004 10 100 slow_TK_2014_12_08.xml T 0.00008 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt0004
mkdir -p results/current_run/l1x1_n2047_i1_s0_O2_dt00025 && ./$folder/src/example 2047 1 0 0.00025 10 200 slow_TK_2014_12_08.xml T 0.00005 O2 && mv *.ex* results/current_run/l1x1_n2047_i1_s0_O2_dt00025
