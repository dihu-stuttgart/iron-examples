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

mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.005 10 25 slow_TK_2014_12_08.xml T 0.005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt0025 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0025 10 50 slow_TK_2014_12_08.xml T 0.0025 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt0025
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.001 10 125 slow_TK_2014_12_08.xml T 0.001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt0005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0005 10 250 slow_TK_2014_12_08.xml T 0.0005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt0005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt00025 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00025 10 500 slow_TK_2014_12_08.xml T 0.00025 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt00025
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt0001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0001 10 1250 slow_TK_2014_12_08.xml T 0.0001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt0001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt00005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00005 10 2500 slow_TK_2014_12_08.xml T 0.00005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt00005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt00001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00001 10 12500 slow_TK_2014_12_08.xml T 0.00001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt00001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt000005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.000005 10 25000 slow_TK_2014_12_08.xml T 0.000005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt000005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M1_dt000001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.000001 10 125000 slow_TK_2014_12_08.xml T 0.000001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M1_dt000001

mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.005 10 25 slow_TK_2014_12_08.xml T 0.001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt0025 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0025 10 50 slow_TK_2014_12_08.xml T 0.0005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt0025
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.001 10 125 slow_TK_2014_12_08.xml T 0.0002 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt0005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0005 10 250 slow_TK_2014_12_08.xml T 0.0001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt0005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt00025 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00025 10 500 slow_TK_2014_12_08.xml T 0.00005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt00025
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt0001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0001 10 1250 slow_TK_2014_12_08.xml T 0.00002 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt0001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt00005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00005 10 2500 slow_TK_2014_12_08.xml T 0.00001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt00005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt00001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00001 10 12500 slow_TK_2014_12_08.xml T 0.000002 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt00001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt000005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.000005 10 25000 slow_TK_2014_12_08.xml T 0.000001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt000005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M5_dt000001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.000001 10 125000 slow_TK_2014_12_08.xml T 0.0000002 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M5_dt000001

mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.005 10 25 slow_TK_2014_12_08.xml T 0.0005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt0025 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0025 10 50 slow_TK_2014_12_08.xml T 0.00025 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt0025
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.001 10 125 slow_TK_2014_12_08.xml T 0.0001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt0005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0005 10 250 slow_TK_2014_12_08.xml T 0.00005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt0005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt00025 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00025 10 500 slow_TK_2014_12_08.xml T 0.000025 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt00025
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt0001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.0001 10 1250 slow_TK_2014_12_08.xml T 0.00001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt0001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt00005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00005 10 2500 slow_TK_2014_12_08.xml T 0.000005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt00005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt00001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.00001 10 12500 slow_TK_2014_12_08.xml T 0.000001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt00001
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt000005 && mpirun -n 4 ./$folder/src/example 64 1 0 0.000005 10 25000 slow_TK_2014_12_08.xml T 0.0000005 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt000005
mkdir -p results/current_run/l1x1_n64_i1_s0_O2_M10_dt000001 && mpirun -n 4 ./$folder/src/example 64 1 0 0.000001 10 125000 slow_TK_2014_12_08.xml T 0.0000001 O2 && mv *.ex* results/current_run/l1x1_n64_i1_s0_O2_M10_dt000001
