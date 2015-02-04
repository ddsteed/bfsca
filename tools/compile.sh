#! /bin/bash
#
# ==================================================
# compile bfsca software suite
#
# Hao Feng (RDS), XHU
#
# ---
# Usage:
#       ./compile.sh
#  or   ./compile.sh   -c mac|intel|gnu  -d obj_dir

# ==================================================
# Options
while test "X$1" != "X" ; do
    case $1 in
	-c) shift; comp=$1     ;;        # compiler vendors
	-d) shift; dest_dir=$1 ;;        # destination of bin|lib files
	-r)        clear=true  ;;        # if clear the OLD temporaries
	*)                     ;;
    esac
    shift
done

# =======
if [ "X$dest_dir" == "X" ]; then  dest_dir=`pwd`/..;   fi
dest_bin=$dest_dir/bin
dest_lib=$dest_dir/lib

# =======
if [ "X$clear" == "X" ]; then clear=false; fi
if [ $clear == true ]; then rm -rf ../build; fi

# compile from out-of-source-tree
if [ ! -d ../build ]; then
    mkdir ../build
fi

cd ../build

# ======
os=`uname`
if [ "X$comp" == "X" ]; then
    case $os in
	Darwin) CC=clang; CXX=clang++; F77=ifort ;;
	Linux)  CC=icc;   CXX=icpc;    F77=ifort ;;
    esac
else
    case $comp in
	mac)   CC=clang; CXX=clang++; F77=gfortran ;;
	gnu)   CC=gcc;   CXX=g++;     F77=gfortran; CCFLAGS="-msse2" ;;
	intel) CC=icc;   CXX=icpc;    F77=ifort ;;
    esac
fi
export CC
export CXX
export F77
export CCFLAGS

# ======
# Let's do it!

echo " --------------------------------------------------------"
echo
echo "       CC = $CC   CXX = $CXX   F77 = $F77                "
echo
echo " --------------------------------------------------------"
echo

cmake ..                            \
    -DCMAKE_Fortran_COMPILER=${F77} \
    -DCMAKE_C_COMPILER=${CC}        \
    -DCMAKE_CXX_COMPILER=${CXX}     \
    -DINSTA_BIN=${dest_bin}         \
    -DINSTA_LIB=${dest_lib}

# (Options
#    -DBLA_VENDOR=${blas_vendor} \
# )

make
make install
