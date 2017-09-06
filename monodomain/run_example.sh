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

#mkdir -p results/current_run/l1x1_n31_i1_s0_O1_dt05 && ./$folder/src/example 31 1 0 0.05 10 1 slow_TK_2014_12_08.xml T 0.05 O1 && mv *.ex* results/current_run/l1x1_n31_i1_s0_O1_dt05
mkdir -p results/current_run/l1x1_n31_i1_s0_O2_dt05 && ./$folder/src/example 31 1 0 0.05 10 1 slow_TK_2014_12_08.xml T 0.05 O2 && mv *.ex* results/current_run/l1x1_n31_i1_s0_O2_dt05

#mkdir -p results/current_run/l1x1_n8_i1_s0_dt05 && ./$folder/src/example 8 1 0 0.05 10 1 slow_TK_2014_12_08.xml T 0.05 && mv *.ex* results/current_run/l1x1_n8_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n16_i1_s0_dt05 && ./$folder/src/example 16 1 0 0.05 10 1 slow_TK_2014_12_08.xml T 0.05 && mv *.ex* results/current_run/l1x1_n16_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n32_i1_s0_dt05 && ./$folder/src/example 32 1 0 0.05 10 1 slow_TK_2014_12_08.xml T 0.05 && mv *.ex* results/current_run/l1x1_n32_i1_s0_dt05

#mkdir -p results/current_run/l1x1_n8_i1_s0_dt05 && ./$folder/src/example 8 1 0 0.05 10 1 slow_TK_2014_12_08.xml F 0.05 && mv *.ex* results/current_run/l1x1_n8_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n16_i1_s0_dt05 && ./$folder/src/example 16 1 0 0.05 10 1 slow_TK_2014_12_08.xml F 0.05 && mv *.ex* results/current_run/l1x1_n16_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n32_i1_s0_dt05 && ./$folder/src/example 32 1 0 0.05 10 1 slow_TK_2014_12_08.xml F 0.05 && mv *.ex* results/current_run/l1x1_n32_i1_s0_dt05

#mkdir -p results/current_run/l1x1_n8_i1_s0_dt05 && ./$folder/src/example 8 1 0 0.05 10 1 hodgkin_huxley_1952.cellml T 0.05 && mv *.ex* results/current_run/l1x1_n8_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n16_i1_s0_dt05 && ./$folder/src/example 16 1 0 0.05 10 1 hodgkin_huxley_1952.cellml T 0.05 && mv *.ex* results/current_run/l1x1_n16_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n32_i1_s0_dt05 && ./$folder/src/example 32 1 0 0.05 10 1 hodgkin_huxley_1952.cellml T 0.05 && mv *.ex* results/current_run/l1x1_n32_i1_s0_dt05

#mkdir -p results/current_run/l1x1_n16_i1_s0_dt05 && ./$folder/src/example 16 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.05 && mv *.ex* results/current_run/l1x1_n16_i1_s0_dt05
#mkdir -p results/current_run/l1x1_n32_i1_s0_dt05 && ./$folder/src/example 32 1 0 0.05 10 1 hodgkin_huxley_1952.cellml F 0.05 && mv *.ex* results/current_run/l1x1_n32_i1_s0_dt05
