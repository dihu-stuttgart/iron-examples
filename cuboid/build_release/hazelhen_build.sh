# build script for hazel hen

export CRAYPE_LINK_TYPE=dynamic
export OPENCMISS_INSTALL_DIR=/lustre/cray/ws8/ws/icbbnmai-iron/install

# load modules
module restore /lustre/cray/ws8/ws/icbbnmai-iron/manage/build_release/gcc49.module_snapshot

ccc && cmake -DOPENCMISS_BUILD_TYPE=RELEASE -DCMAKE_BUILD_TYPE=RELEASE  .. && make clean && make all

export OPENCMISS_REL_DIR=/zhome/academic/HLRS/cbm/icbbnmai/opencmiss/examples/cuboid/build_release/Fortran
