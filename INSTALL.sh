#!/bin/bash
#This code follows the Cluster_install.odt instructions. 
#Currently it is set up to install all packages locally using the Code_Cluster folder.
#To install just run this script using: $bash INSTALL.sh 
#(do not copy the $ or this parethesis, also make sure you are in the Code_Cluster folder)

GCC=false
CMAKE=false
GMP=false
MPFR=false
BOOST=false
CGAL=false
ZLIB=false
HDF5=false
NETCDFC=false
NETCDF42=false
NETCDF43=false
DOALL=false
DOWNLOAD=false

helpFunction()
{
   echo ""
   echo "Usage: bash $0 -OPTION1 STATEMENT1 -OPTION2 STATEMENT2 ..."
   echo "The option command controls whether a specific package is installed (true) or not (false)."
   echo "The STATEMENT is either true/false for the packages. It is considered false as default."
   echo -e "\t-a for gcc-4.8.5"
   echo -e "\t-b for cmake-3.13.1"
   echo -e "\t-c for gmp-6.1.2"
   echo -e "\t-d for mpfr-4.0.1"
   echo -e "\t-e for boost-1.62.0"
   echo -e "\t-f for CGAL-4.9"
   echo -e "\t-g for zlib-1.2.11"
   echo -e "\t-h for hdf5-1.10.4"
   echo -e "\t-i for netcdf-c-4.6.2"
   echo -e "\t-j for netcdf-cxx-4.2"
   echo -e "\t-k for netcdf-cxx4-4.3.0"
   echo -e "\t-u if the download of the tar files is also required"
   echo -e "\t-s if netcdf-c-4.6.2 is false, then give its install directory for netcdf-cxx-4.2 or netcdf-cxx4-4.3.0"
   echo -e "\t-t if hdf5-1.10.4 is false, then give its install directory for netcdf-c-4.6.2"
   echo -e "\t-v if zlib-1.2.11 is false, then give its install directory for hdf5-1.10.4 or netcdf-c-4.6.2"
   echo -e "\t-w if boost-1.62.0 is false, then give its install directory for CGAL-4.9"
   echo -e "\t-x if gmp-6.1.2 is false, then give its install directory for mpfr-4.0.1"
   echo -e "\t-y install all packages"
   echo -e "\t-z for help"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:d:e:f:g:h:i:j:k:s:t:v:w:x:y:z:" opt
do
   case "$opt" in
      a ) GCC="$OPTARG" ;;
      b ) CMAKE="$OPTARG" ;;
      c ) GMP="$OPTARG" ;;
      d ) MPFR=$OPTARG ;;
      e ) BOOST="$OPTARG" ;;
      f ) CGAL="$OPTARG" ;;
      g ) ZLIB=$OPTARG ;;
      h ) HDF5="$OPTARG" ;;
      i ) NETCDFC="$OPTARG" ;;
      j ) NETCDF42=$OPTARG ;;
      k ) NETCDF43="$OPTARG" ;;
      s ) NETCDFDIR="$OPTARG" ;;
      t ) HDF5DIR="$OPTARG" ;;
      v ) ZLIBDIR="$OPTARG" ;;
      w ) BOOSTDIR="$OPTARG" ;;
      x ) GMPDIR="$OPTARG" ;;
      y ) DOALL="$OPTARG" ;;
      z ) helpFunction ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ "$DOALL" = true ] ;then
	GCC=true
	CMAKE=true
	GMP=true
	MPFR=true
	BOOST=true
	CGAL=true
	ZLIB=true
	HDF5=true
	NETCDFC=true
	NETCDF42=true
	NETCDF43=true
fi

# Print helpFunction in case parameters are empty
#if [ -z "$parameterA" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ]
#then
#   echo "Some or all of the parameters are empty";
#   helpFunction
#fi

if [ "$GCC" = true ] ;then
	echo -e "WARNING: THIS INSTALLATION MAY TAKE SEVERAL HOURS. \n"
fi

echo -e "\nHere are the install parameters:\n"

echo "gcc-4.8.5: $GCC"
echo "cmake-3.13.1: $CMAKE"
echo "gmp-6.1.2: $GMP"
echo "mpfr-4.0.1: $MPFR"
echo "boost-1.62.0: $BOOST"
echo "CGAL-4.9: $CGAL"
echo "zlib-1.2.11: $ZLIB"
echo "hdf5-1.10.4: $HDF5"
echo "netcdf-c-4.2: $NETCDFC"
echo "netcdf-cxx-4.2: $NETCDF42"
echo "netcdf-cxx5-4.3.0: $NETCDF43"
echo "download packages: $DOWNLOAD"

echo -e "\nThe install will start in 30 seconds, check if all above variables are correct. If not then type CTRL-C, otherwise, see you on the other side!"
(sleep 30)

BASEDIR=$PWD
mkdir CellGPUPackages #where the install will be
mkdir CellGPULib #where all the downloads will be

if [ "$DOWNLOAD" = true ] ;then
	cd CellGPULib
	wget "https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz"
	tar -xzvf eigen-3.3.7.tar.gz
	cd ..
fi

if [ "$GCC" = true ] ;then

	echo -e "\nStart gcc-4.8.5 install"
	cd "$BASEDIR"/CellGPUPackages/
	rm -r gcc-4.8.5
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://ftp.gnu.org/gnu/gcc/gcc-4.8.5/gcc-4.8.5.tar.gz"
	fi
	tar -xzvf gcc-4.8.5.tar.gz
	cd gcc-4.8.5
	mkdir build && cd build
	../configure --prefix="$BASEDIR"/CellGPUPackages/gcc-4.8.5 --enable-languages=c,c++,fortran,go --disable-multilib
	make
	make install

	echo "" >> ~/.bashrc
	echo "#New gcc version" >> ~/.bashrc
	echo "export PATH=$BASEDIR/CellGPUPackages/gcc-4.8.5/bin:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=$BASEDIR/CellGPUPackages/gcc-4.8.5/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=$BASEDIR/CellGPUPackages/gcc-4.8.5/lib64:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$CMAKE" = true ] ;then

	echo -e "\nStart cmake-3.13.1 install"
	cd "$BASEDIR"/CellGPUPackages/
	rm -r cmake-3.13.1
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://cmake.org/files/v3.13/cmake-3.13.1.tar.gz"
	fi
	tar -xzvf cmake-3.13.1.tar.gz
	cd cmake-3.13.1
	mkdir build && cd build
	../bootstrap --prefix="$BASEDIR"/CellGPUPackages/cmake-3.13.1
	make
	make install

	echo "" >> ~/.bashrc
	echo "#New cmake version" >> ~/.bashrc
	echo "export PATH=$BASEDIR/CellGPUPackages/cmake-3.13.1/bin:\$PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$GMP" = true ] ;then

	echo -e "\nStart gmp-6.1.2 install"
	cd "$BASEDIR"/CellGPUPackages/
	rm -r gmp-6.1.2
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz"
	fi
	tar -xf gmp-6.1.2.tar.xz
	cd gmp-6.1.2
	mkdir build
	cd build
	../configure --prefix="$BASEDIR"/CellGPUPackages/gmp-6.1.2
	make
	make install

	GMPDIR="$BASEDIR"/CellGPUPackages/gmp-6.1.2
	echo "" >> ~/.bashrc
	echo "#New gmp version" >> ~/.bashrc
	echo "export PATH=$BASEDIR/CellGPUPackages/gmp-6.1.2:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=$BASEDIR/CellGPUPackages/gmp-6.1.2/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$MPFR" = true ] ;then

	echo -e "\nStart mpfr-4.0.1 install"

	if [ -z "$GMPDIR" ]; then
		echo "Need gmp-6.1.2 install directory"
		helpFunction
	fi

	cd "$BASEDIR"/CellGPUPackages/
	rm -r mpfr-4.0.1
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://www.mpfr.org/mpfr-4.0.1/mpfr-4.0.1.tar.xz"
	fi
	tar -xf mpfr-4.0.1.tar.xz
	cd mpfr-4.0.1
	mkdir build && cd build
	../configure --prefix="$BASEDIR"/CellGPUPackages/mpfr-4.0.1 --with-gmp="$GMPDIR"
	make
	make install

	echo "" >> ~/.bashrc
	echo "#New mpfr version" >> ~/.bashrc
	echo "export PATH=~/Code_Cluster/CellGPUPackages/mpfr-4.0.1:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/mpfr-4.0.1/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$BOOST" = true ] ;then

	echo -e "\nStart boost-1.62.0 install"
	cd "$BASEDIR"/CellGPUPackages/
	rm -r boost-1.62.0
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz/download"
	fi
	tar -xzf boost_1_62_0.tar.gz
	cd boost_1_62_0
	./bootstrap.sh --prefix="$BASEDIR"/CellGPUPackages/boost-1.62.0 
	./b2
	./b2 install

	BOOSTDIR="$BASEDIR"/CellGPUPackages/boost-1.62.0 
	echo "" >> ~/.bashrc
	echo "#New boost version" >> ~/.bashrc
	echo "export PATH=~/Code_Cluster/CellGPUPackages/boost-1.62.0:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/boost-1.62.0/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$CGAL" = true ] ;then

	echo -e "\nStart CGAL-4.9 install"

	if [ -z "$BOOSTDIR" ]; then
		echo "Need boost-1.62.0 install directory"
		helpFunction
	fi

	cd "$BASEDIR"/CellGPUPackages/
	rm -r CGAL-4.9
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.9/CGAL-4.9.tar.xz"
	fi
	tar -xf CGAL-4.9.tar.xz
	cd  CGAL-4.9
	mkdir build && cd build
	cmake -DCMAKE_INSTALL_PREFIX="$BASEDIR"/CellGPUPackages/CGAL-4.9 -DBOOST_INCLUDEDIR="$BOOSTDIR" -DBOOST_LIBRARYDIR="$BOOSTDIR"/lib/  -DCGAL_HEADER_ONLY=OFF -DCMAKE_BUILD_TYPE=Release .. 
	make
	make install

	echo "" >> ~/.bashrc
	echo "#New CGAL version" >> ~/.bashrc
	echo "export PATH=~/Code_Cluster/CellGPUPackages/CGAL-4.9:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/CGAL-4.9/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/CGAL-4.9/lib64:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$ZLIB" = true ] ;then

	echo -e "\nStart zlib-1.2.11 install"
	cd "$BASEDIR"/CellGPUPackages/
	rm -r zlib-1.2.11
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://www.zlib.net/zlib-1.2.11.tar.gz"
	fi
	tar -xzf zlib-1.2.11.tar.gz
	cd zlib-1.2.11
	mkdir build && cd build
	../configure --prefix="$BASEDIR"/CellGPUPackages/zlib-1.2.11
	make
	make install
	ZLIBDIR="$BASEDIR"/CellGPUPackages/zlib-1.2.11
fi

if [ "$HDF5" = true ] ;then

	echo -e "\nStart hdf5-1.10.4 install"

	if [ -z "$ZLIBDIR" ]; then
		echo "Need zlib-1.2.11 install directory"
		helpFunction
	fi

	cd "$BASEDIR"/CellGPUPackages/
	rm -r hdf5-1.10.4
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz"
	fi
	tar -xzf hdf5-1.10.4.tar.gz
	cd hdf5-1.10.4
	mkdir build && cd build
	../configure --prefix="$BASEDIR"/CellGPUPackages/hdf5-1.10.4 --with-zlib="$ZLIBDIR" --enable-hl
	make
	make install
	HDF5DIR="$BASEDIR"/CellGPUPackages/hdf5-1.10.4
fi

if [ "$NETCDFC" = true ] ;then

	echo -e "\nStart netcdf-c-4.6.2 install"

	if [ -z "$ZLIBDIR" ] ; then
		echo "Need zlib-1.2.11 install directory"
		helpFunction
	fi

	if [ -z "$HDF5DIR" ] ; then
		echo "Need hdf5-1.10.4 install directory"
		helpFunction
	fi

	cd "$BASEDIR"/CellGPUPackages/
	rm -r "netcdf-c-4.6.2"
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://github.com/Unidata/netcdf-c/archive/v4.6.2.tar.gz"
	fi
	tar -xzf "netcdf-c-4.6.2.tar.gz"
	cd "netcdf-c-4.6.2"
	mkdir build && cd build
	CPPFLAGS='-I"$HDF5DIR"/include -I"$ZLIBDIR"/include' LDFLAGS='-L"$HDF5DIR"/lib -L"$ZLIBDIR"/lib' ../configure --prefix="$BASEDIR/CellGPUPackages/netcdf-c-4.6.2"
	make
	make install
	NETCDFDIR="$BASEDIR/CellGPUPackages/netcdf-c-4.6.2"

	echo "" >> ~/.bashrc
	echo "#New NetCDF version" >> ~/.bashrc
	echo "export PATH=~/Code_Cluster/CellGPUPackages/netcdf-c-4.6.2:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/netcdf-c-4.6.2/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$NETCDF42" = true ] ;then

	echo -e "\nStart netcdf-cxx-4.2 install"

	if [ -z "$NETCDFDIR" ] ; then
		echo "Need netcdf-c-4.6.2 install directory"
		helpFunction
	fi

	cd "$BASEDIR"/CellGPUPackages/
	rm -r netcdf-cxx-4.2
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx-4.2.tar.gz"
	fi
	tar -xzf netcdf-cxx-4.2.tar.gz
	cd netcdf-cxx-4.2
	mkdir build && cd build
	CPPFLAGS='-I"$NETCDFDIR"/include' LDFLAGS='-L"$NETCDFDIR"/lib' ../configure --prefix="$BASEDIR"/CellGPUPackages/netcdf-cxx-4.2
	make
	make install

	echo "" >> ~/.bashrc
	echo "#New NetCDF C++ version" >> ~/.bashrc
	echo "export PATH=~/Code_Cluster/CellGPUPackages/netcdf-cxx-4.2:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/netcdf-cxx-4.2/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

if [ "$NETCDF43" = true ] ;then

	echo -e "\nStart netcdf-cxx4-4.3.0 install"

	if [ -z "$NETCDFDIR" ] ; then
		echo "Need netcdf-c-4.6.2 install directory"
		helpFunction
	fi

	cd "$BASEDIR"/CellGPUPackages/
	rm -r netcdf-cxx4-4.3.0
	cd "$BASEDIR"/CellGPULib/
	if [ "$DOWNLOAD" = true ] ;then
		wget "https://github.com/Unidata/netcdf-cxx4/archive/v4.3.0.tar.gz"
	fi
	tar -xzf netcdf-cxx4-4.3.0.tar.gz
	cd netcdf-cxx4-4.3.0
	mkdir build && cd build
	CPPFLAGS='-I"$NETCDFDIR"/include' LDFLAGS='-L"$NETCDFDIR"/lib' ../configure --prefix="$BASEDIR"/CellGPUPackages/netcdf-cxx4-4.3.0
	make
	make install

	echo "" >> ~/.bashrc
	echo "#New New NetCDF C++ version" >> ~/.bashrc
	echo "export PATH=~/Code_Cluster/CellGPUPackages/netcdf-cxx4-4.3.0:\$PATH" >> ~/.bashrc
	echo "export LD_LIBRARY_PATH=~/Code_Cluster/CellGPUPackages/netcdf-cxx4-4.3.0/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
	. ~/.bashrc
fi

echo -e "\nInstalation Complete! Congratulations!"
