#!/bin/bash -e

# ================== #
# compiling defaults
# ================== #
BUILD_3PL=0
BUILD_SOLVER=0
BUILD_TYPE=0
CLEAN_DIST=0

# ============== #
# print strings
# ============== #
opt_str="[OPTION] "

eC="\x1B[0m"
bC="\x1B[0;34m"
GC="\x1B[1;32m"
yC="\x1B[0;33m"
aC="\x1B[0;96m"
rC="\x1B[0;41m"
gC="\x1B[0;32m"
oC="\x1B[3;93m"
mC="\x1B[0;43m"

help() {
    echo -e " ==================================================== "
    echo -e " ${GC} >>>>>>  Easy Build Option:  ./makescript.sh -go  <<<<<<${eC}"
    echo -e " ==================================================== "
    echo " "
    echo -e "${GC} Usage:${eC} $0 [OPTION]...[COMPILER OPTIONS]...[3PL OPTIONS]"
    echo " "
    echo -e " ${aC}Recommended Options:${eC}"
    echo -e "    Default: ./makescript.sh -go"
    echo "      ./makescript.sh CC=mpicc CXX=mpicxx FC=mpif90"
    echo " "
    echo "    NWSC Cheyenne: ./makescript.sh -nwsc "
    echo "      ./makescript.sh CC=mpicc CXX=mpicxx FC=mpif90 --avx2 -mkl"
    echo " "
    echo "    Intel MPI: ./makescript.sh -impi"
    echo "      ./makescript.sh CC=mpiicc CXX=mpiicpc FC=mpiifort"
    echo " "
    echo -e " ${aC}Options List:${eC}"
    echo "  [OPTION]:"
    echo "    --3pl       -3pl       build the 3rd party libraries: metis, t8code, occa, googletest"
    echo "    --t8wyo     -t8wyo     build dg4est amr library"
    echo " "
    echo "    --help      -h      displays this help message"
    echo "    --clean     -c      removes local build directories"
    echo "    --distclean -dc     removes builds and install directories"
    echo "    --release   -opt    compile the project in optimized mode"
    echo "    --debug     -deb    compile the project in debug mode"
    echo "    --testsON   -ton    turn on unit tests (google tests)"
    echo " "
    echo "  [COMPILER OPTIONS]:"
    echo "     CC=<arg>   cc=<arg>    sets the C compiler"
    echo "    CXX=<arg>  cxx=<arg>    sets the C++ compiler"
    echo "     FC=<arg>   fc=<arg>    sets the Fortran compiler"
    echo " "
    echo "      C_FLAGS=<arg>    c_flags=<arg>    sets the C compiler flags"
    echo "    CXX_FLAGS=<arg>  cxx_flags=<arg>    sets the C++ compiler flags"
    echo "     FC_FLAGS=<arg>   fc_flags=<arg>    sets the Fortran compiler flags"
    echo " "
    echo "  [Intel Flags]:"
    echo "     --avx     -avx         sets Intel AVX Instructions"
    echo "     --avx2    -avx2        sets Intel AVX-2 Instructions"
    echo "     --avx512  -avx512      sets Intel AVX-512 Instructions"
    echo "     --mkl     -mkl         sets Intel MKL Library"
    echo " "
    echo "  [3PL OPTIONS]:"
    echo "    --ALL3P   -all3p  compile all 3rd party libraries"
    echo "    --T8CODE  -t8     compile t8code"
    echo "    --METIS   -metis  compile metis"
    echo " "
}

# ============ #
# parse inputs
# ============ #

if [[ $# -lt 1 ]]; then
  help
  exit 1
fi

for var in "$@"
do
  if [ "$var" == "--help" -o "$var" == "-help" -o "$var" == "-h" ]; then
    help
    exit 0

  elif [ "$var" == "--distclean" -o "$var" == "-distclean" -o "$var" == "-dc" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    echo ${opt_str} "Cleaning the distribution..."
    CLEAN_DIST=1

  elif [ "$var" == "--3pl" -o "$var" == "-3pl" -o "$var" == "-3pl" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    BUILD_3PL=1

  elif [ "$var" == "--T8CODE" -o "$var" == "-t8" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    BUILD_3PL=1

  elif [ "$var" == "--METIS" -o "$var" == "-metis" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    BUILD_3PL=1

  elif [ "$var" == "--t8wyo" -o "$var" == "-t8wyo" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    BUILD_SOLVER=1

  elif [ "$var" == "--release" -o "$var" == "-release" -o "$var" == "-opt" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    BUILD_TYPE=0

  elif [ "$var" == "--debug" -o "$var" == "-debug" -o "$var" == "-deb" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    BUILD_TYPE=1

  elif [ "$var" == "--clean" -o "$var" == "-clean" -o "$var" == "-c" -o \
         "$var" == "--testsON" -o "$var" == "-testsON" -o "$var" == "-ton" -o \
         "${var:0:3}" == "CC=" -o "${var:0:3}" == "cc=" -o \
         "${var:0:4}" == "CXX=" -o "${var:0:4}" == "cxx=" -o \
         "${var:0:3}" == "FC=" -o "${var:0:3}" == "fc=" -o \
         "${var:0:8}" == "C_FLAGS=" -o "${var:0:8}" == "c_flags=" -o \
         "${var:0:10}" == "CXX_FLAGS=" -o "${var:0:10}" == "cxx_flags=" -o \
         "${var:0:9}" == "FC_FLAGS=" -o "${var:0:9}" == "fc_flags=" -o \
         "${var}" == "--avx" -o  "${var}" == "-avx" -o \
         "${var}" == "--avx2" -o  "${var}" == "-avx2" -o \
         "${var}" == "--avx512" -o  "${var}" == "-avx512" -o \
         "${var:0:13}" == "solver_so_2d=" -o \
         "${var:0:13}" == "solver_so_3d=" -o \
         "${var:0:11}" == "solver_dir=" -o \
         "$var" == "--ALL3PL" -o "$var" == "--all3pl" -o "$var" == "-all3pl" -o \
         "$var" == "--T8CODE" -o "$var" == "-t8code" -o \
         "$var" == "--METIS" -o "$var" == "-metis" -o \
         "${var}" == "--mkl" -o  "${var}" == "-mkl" -o \
         "${var}" == "-go" -o \
         "${var}" == "-nwsc" -o \
         "${var}" == "-impi" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"

  else
    echo -e "${oC}Unknown option:${eC}  ${rC}$var${eC}"
    echo "See available options: ./makescript.sh -help"
    echo "Using Default Options..."

  fi
done

# ========================= #
# display command line args
# ========================= #
echo "$0 $@"
cmd_args="${@:1}"

# =================================================================== #
if [ $CLEAN_DIST == 1 ]; then
  rm -rf install
  rm -rf builds
  exit 0
fi
# =================================================================== #

# =================================================================== #
if [ $BUILD_SOLVER == 0 -a $BUILD_3PL == 0 ]; then
  echo "=================================================================="
  echo "Building the Metis, t8code, and t8wyo..."
  echo "=================================================================="
  echo " "

  cd scripts

  # build 3PL libraries
  ./build_3PL.sh $cmd_args

  # build solver
  ./build_solver.sh $cmd_args

  cd ..

  echo
  echo "================================================"
  echo -e "${gC} Finished Successfully...${eC}"
  echo -e " Executable: ${GC}builds/t8wyo_release/bin/t8wyo?d.mpi${eC}"
  echo "================================================"
  exit 0
fi

# ================================ #
# BUILD INDIVIDUAL COMPONENTS ONLY #
# ================================ #
if [ $BUILD_3PL == 1 ]; then
  echo "========================================================"
  echo "Building the 3rd Party Libraries..."
  echo "========================================================"

  cd scripts
  ./build_3PL.sh $cmd_args
  cd ..
fi

if [ $BUILD_SOLVER == 1 ]; then
  echo "================="
  echo "Building t8wyo..."
  echo "================="

  cd scripts
  ./build_solver.sh $cmd_args
  cd ..

  echo
  echo "================================================"
  echo -e "${gC} Finished Successfully...${eC}"
  echo -e " Executable: ${GC}builds/t8wyo_release/bin/t8wyo?d.mpi${eC}"
  echo "================================================"
fi
# =================================================================== #