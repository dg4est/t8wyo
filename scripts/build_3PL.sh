#!/bin/bash -e
# Folder structure:
# ================
# project_root
#   install/        # install directory
#     3PL/          # install to this folder
#   builds/         # build directory
#   3PL/.           # location of 3rd party library source code
#   src/            # source codes
#   scripts/.       # location of this script


# ======================= #
# project directory paths
# ======================= #
CURRENT_PATH="$(pwd)"
MY_PATH="$( cd "$( dirname "$0" )" && pwd )"
PROJECT_ROOT=${MY_PATH}/..

# ====================== #
# folder directory paths
# ====================== #
INSTALL_DIRECTORY=${PROJECT_ROOT}/install
INSTALL_3PL_DIRECTORY=${INSTALL_DIRECTORY}/3PL

BUILD_DIRECTORY=${PROJECT_ROOT}/builds

# =========== #
# 3PL sources #
# =========== #
SOURCES_3PL_DIRECTORY=${PROJECT_ROOT}/3PL

# ================== #
# compiling defaults #
# ================== #
BUILD_T8CODE=0
BUILD_METIS=0
BUILD_CLEAN=0

# ================= #
# compiler defaults #
# ================= #
CC=mpicc
CXX=mpicxx
FC=mpif90
F77=mpif90

CFLAGS="-fPIC -O2 -std=gnu99"

# ============ #
# BLAS library #
# ============ #
INSTALL_BLAS_DIRECTORY="${MKLROOT}"
BLAS_PATH="-Wl,--start-group \
${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
${MKLROOT}/lib/intel64/libmkl_core.a \
${MKLROOT}/lib/intel64/libmkl_sequential.a \
-Wl,--end-group -lpthread -lm"


# ======================== #
# compiler option defaults #
# ======================== #
BUILD_SUFFIX="_release"
BUILD_TYPE="Release"

# ======================== #
# make and install command #
# ======================== #
MAKE_CMD="make -j4 install"

# ============= #
# print strings #
# ============= #
opt_str="[OPTION] "

eC="\x1B[0m"
rC="\x1B[0;41m"
gC="\x1B[0;42m"
yC="\x1B[0;33m"
mC="\x1B[0;43m"

help() {
    echo "Usage: $0 [OPTION]...[COMPILER OPTIONS]...[3PL OPTIONS]"
    echo " "
    echo "  This script builds the 3rd party libraries: "
    echo "       t8code, metis"
    echo " "
    echo "  [OPTION]:"
    echo "    --help    -h      displays this help message"
    echo "    --clean   -c      removes build directory: dg4est/builds/dg4est_{}"
    echo "    --release -opt    compile the project in optimized mode"
    echo "    --debug   -deb    compile the project in debug mode"
    echo " "
    echo "  [COMPILER OPTIONS]:"
    echo "     CC=<arg>     cc=<arg>    sets the C compiler"
    echo "    CXX=<arg>    cxx=<arg>    sets the C++ compiler"
    echo "     FC=<arg>     fc=<arg>    sets the Fortran compiler"
    echo " "
    echo "      C_FLAGS=<arg>    c_flags=<arg>    sets the C compiler flags"
    echo "    CXX_FLAGS=<arg>  cxx_flags=<arg>    sets the C++ compiler flags"
    echo "     FC_FLAGS=<arg>   fc_flags=<arg>    sets the Fortran compiler flags"
    echo " "
    echo "  [Intel Flags]:"
    echo "     --avx     -avx         sets Intel AVX Instructions"
    echo "     --avx2    -avx2        sets Intel AVX-2 Instructions"
    echo "     --avx512  -avx512      sets Intel AVX-512 Instructions"
    echo " "
    echo "  [3PL OPTIONS]:"
    echo "    --ALL3PL  -all3pl comile all 3rd party libraries"
    echo "    --T8CODE  -t8     compile t8code"
    echo "    --METIS   -metis  compile metis"
    echo " "
    echo -e "  ${aC}Recommended Options:${eC}"
    echo -e "    Default (-go): ${yC}./build_3PL.sh -go${eC}"
    echo "        CC=mpicc CXX=mpicxx FC=mpif90 -all3pl"
    echo "  "
    echo -e "    Intel MPI (-impi): ${yC}./build_3PL.sh -impi${eC}"
    echo "        CC=mpiicc CXX=mpiicpc FC=mpiifort -all3pl"
}

# ----------------------------- #
# Start the compilation process #
# ----------------------------- #
cd $PROJECT_ROOT

# ============ #
# parse inputs #
# ============ #
for var in "$@"
do
  if [ "$var" == "--help" -o "$var" == "-help" -o "$var" == "-h" ]; then
    help
    exit 0
  elif [ "$var" == "--clean" -o "$var" == "-clean" -o "$var" == "-c" ]; then
    echo ${opt_str} "Clean and rebuild"
    BUILD_CLEAN=1

  elif [ "$var" == "--release" -o "$var" == "-release" -o "$var" == "-opt" ]; then
    echo ${opt_str} "Compiling in optimized mode"
    BUILD_SUFFIX="_release"
    BUILD_TYPE="Release"

  elif [ "$var" == "--debug" -o "$var" == "-debug" -o "$var" == "-deb" ]; then
    echo ${opt_str} "Compiling in debug mode"
    BUILD_SUFFIX="_debug"
    BUILD_TYPE="Debug"

  elif [ "${var:0:3}" == "CC=" -o "${var:0:3}" == "cc=" ]; then
    CC=${var:3}
    echo -e "[OPTION]       C Compiler: $yC$CC$eC"

  elif [ "${var:0:4}" == "CXX=" -o "${var:0:4}" == "cxx=" ]; then
    CXX=${var:4}
    echo -e "[OPTION]     CXX Compiler: $yC$CXX$eC"

  elif [ "${var:0:3}" == "FC=" -o "${var:0:3}" == "fc=" ]; then
    FC=${var:3}
    F77=${FC}
    echo -e "[OPTION] Fortran Compiler: $yC$FC$eC"

  elif [ "$var" == "--T8CODE" -o "$var" == "-t8" ]; then
    BUILD_T8CODE=1
  elif [ "$var" == "--METIS" -o "$var" == "-metis" ]; then
    BUILD_METIS=1
  elif [ "$var" == "-go" ]; then
    CC=mpicc
    CXX=mpicxx
    FC=mpif90
  elif [ "$var" == "-impi" ]; then
    CC=mpiicc
    CXX=mpiicpc
    FC=mpiifort
  elif [ "$var" == "--ALL3PL" -o "$var" == "--all3pl" -o "$var" == "-all3pl" ]; then
    BUILD_METIS=1
    BUILD_T8CODE=1
  fi
done

# if no 3PL are selected, compile all of them
if [ $BUILD_METIS == 0 -a $BUILD_T8CODE == 0 ]; then
  BUILD_METIS=1
  BUILD_T8CODE=1
fi

# ========================= #
# display command line args #
# ========================= #
echo " "
echo "$0 $@"

# ----------------------------------------------------- #
# After reading in cmd arg options, set remaining paths #
# ----------------------------------------------------- #

# ====================================== #
# install/build location compiled source #
# ====================================== #
COMPILE_INSTALL_3PL_DIRECTORY="${INSTALL_3PL_DIRECTORY}${BUILD_SUFFIX}"
COMPILE_BUILD_3PL_DIRECTORY="${BUILD_3PL_DIRECTORY}${BUILD_SUFFIX}"

# ============== #
# compiler paths #
# ============== #
FC_PATH="`which $FC`"
CC_PATH="`which $CC`"
CXX_PATH="`which $CXX`"
LD_PATH="`which ld`"

# ====================== #
# check source directory
# ====================== #
if [ ! -d "${SOURCES_3PL_DIRECTORY}" ]; then
  echo "${rC}ERROR: {SOURCES_3PL_DIRECTORY} does not exist.${eC}"
  exit 1
fi

# ======================= #
# check install directory
# ======================= #
if [ ! -d "${INSTALL_DIRECTORY}" ]; then
  echo  "${INSTALL_DIRECTORY} does not exist. Making it..."
  mkdir "${INSTALL_DIRECTORY}"
fi

# ====================== #
# check builds directory
# ====================== #
if [ ! -d "${BUILD_DIRECTORY}" ]; then
  echo  "${BUILD_DIRECTORY} does not exist. Making it..."
  mkdir "${BUILD_DIRECTORY}"
fi


# =================================================================== #
# METIS BUILD
# =================================================================== #
COMPILE_FAIL=0
INSTALL_METIS_DIRECTORY=${INSTALL_3PL_DIRECTORY}/metis
if [ ${BUILD_METIS} == 1 ]; then
  echo " "
  echo -e "${mC} ==== Building Metis ==== ${eC}"
  echo " Compiling Options:"
  echo "        Build Type: ${BUILD_TYPE}"
  echo "  Install Location: ${INSTALL_METIS_DIRECTORY}"
  echo " "
  echo "                CC: ${CC}"
  echo "               CXX: ${CXX}"
  echo "                FC: ${FC}"
  echo -e "${mC} ========================= ${eC}"
  echo " "

  cd ${SOURCES_3PL_DIRECTORY}/metis-5.0.3
  rm -rf buildW
  rm -rf ${INSTALL_METIS_DIRECTORY}
  mkdir buildW

  cd buildW
  cmake -D CMAKE_INSTALL_PREFIX:PATH=${INSTALL_METIS_DIRECTORY}       \
        -D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}                      \
        -D GKLIB_PATH:PATH=${SOURCES_3PL_DIRECTORY}/metis-5.0.3/GKlib \
        -D SHARED:BOOL=TRUE                                           \
        -G "Unix Makefiles" ../

  # make metis
  ${MAKE_CMD}

  # remove build directory
  cd ..
  rm -rf buildW
  cd ${CURRENT_PATH}

  if [ ! -d "${INSTALL_METIS_DIRECTORY}" ]; then
    echo "ERROR:"
    echo "${INSTALL_METIS_DIRECTORY} does not exist."
    COMPILE_FAIL=1
  fi

  if [ ${COMPILE_FAIL} == 0 ]; then
    echo " "
    echo         "========================="
    echo -e "${gC} Metis build successful! ${eC}"
    echo         "========================="
    echo " "
  else
    echo " "
    echo         "==========================="
    echo -e "${rC} Metis build FAILED! ${eC}"
    echo         "==========================="
    echo " "
    exit 1
  fi
fi
# =================================================================== #

# =================================================================== #
# T8CODE BUILD
# =================================================================== #
COMPILE_FAIL=0
INSTALL_T8CODE_DIRECTORY=${INSTALL_3PL_DIRECTORY}/t8code
if [ ${BUILD_T8CODE} -eq 1 ]; then
  echo " "
  echo -e "${mC} ==== Building t8code ==== ${eC}"
  echo " Compiling Options:"
  echo "        Build Type: ${BUILD_TYPE}"
  echo "  Install Location: ${INSTALL_T8CODE_DIRECTORY}"
  echo " "
  echo "                CC: ${CC}"
  echo "               CXX: ${CXX}"
  echo "                FC: ${FC}"
  echo -e "${mC} ========================= ${eC}"
  echo " "

  if [ ! -d "${INSTALL_METIS_DIRECTORY}" ]; then
    echo "Error:"
    echo "${INSTALL_METIS_DIRECTORY} does not exist."
    echo " Please build Metis first: ./build_3PL.sh metis"
    exit
  fi

  cd ${PROJECT_ROOT}
  git submodule init    # init t8wyo submodules
  git submodule update  # update t8wyo submodules

  cd ${SOURCES_3PL_DIRECTORY}/t8code
  rm -rf ${INSTALL_T8CODE_DIRECTORY}

  git submodule init    # init t8code submodules
  git submodule update  # update t8code submodules

  # build t8code #
  cd ${SOURCES_3PL_DIRECTORY}/t8code
  ./bootstrap
  ./configure --prefix=${INSTALL_T8CODE_DIRECTORY}                  \
              --exec-prefix=${INSTALL_T8CODE_DIRECTORY}             \
              --enable-shared --disable-static                      \
              --disable-debug --enable-mpi                          \
              --without-blas                                        \
              --with-metis=${INSTALL_METIS_DIRECTORY}               \
              "CC=$CC" "CXX=$CXX" "FC=$FC" "F77=$FC"                \
              "CFLAGS=$CFLAGS"                                      \
              | tee config.out

  #           --with-vtk=<VTK_LIBS>
  #           --with-occ=<OCC_LIBS> # OpenCASCADE High Order Geometry

  ${MAKE_CMD}
  cd ${CURRENT_PATH}

  if [ ! -d "${INSTALL_T8CODE_DIRECTORY}" ]; then
    echo "ERROR:"
    echo "${INSTALL_T8CODE_DIRECTORY} does not exist."
    COMPILE_FAIL=1
  fi

  if [ ${COMPILE_FAIL} == 0 ]; then
    echo " "
    echo         "=========================="
    echo -e "${gC} t8code build successful! ${eC}"
    echo         "=========================="
    echo " "
  else
    echo " "
    echo         "============================"
    echo -e "${rC} t8code build FAILED! ${eC}"
    echo         "============================"
    echo " "
    exit 1
  fi
fi
# =================================================================== #
echo -e "${gC}Build 3PL Script Completed Successfully!${eC}"