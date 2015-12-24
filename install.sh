#!/bin/bash
LOG_FILE=./LOG_INSTALL.log


# LIBFORBES INSTALLATION SCRIPT
# Licence note:
#
# ForBES is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#   
# ForBES is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with ForBES. If not, see <http://www.gnu.org/licenses/>.


#						 #
# Do not modify this file after this line	 #
#						 #


# Function: install_libforbes
function install_libforbes {

  echo "Checking for dependencies:"
  if hash wget 2>/dev/null; then
       echo "wget .............. OK"
  else
       echo "wget .............. NOT FOUND"
  fi
  if hash g++ 2>/dev/null; then
       echo "g++  .............. OK"
  else
       echo "g++  .............. NOT FOUND"
  fi
  if hash make 2>/dev/null; then
       echo "make .............. OK"
  else
       echo "make .............. NOT FOUND"
  fi
  if  [ -e /usr/lib/libopenblas.a ]; then
       echo "openblas .......... OK";
  else
       echo "openblas .......... NOT FOUND";
  fi
  if [ -d libs/SuiteSparse ]; then
        echo "SuiteSparse ....... OK";
  else
       echo "SuiteSparse ....... NOT FOUND"
  fi
  if [ -e /usr/lib/liblapack.a ]; then
       echo "lapack ............ OK";
  else
       echo "lapack ............ NOT FOUND";
  fi
  if [ -e /usr/lib/liblapacke.a ]; then
       echo "Lapacke ........... OK";
  else
       echo "Lapacke ........... NOT FOUND";
  fi


  mkdir -p libs

  # Check if the SuiteSparse directory exists; if not, download SuiteSparse,
  # untar it and delete the tar file
  echo "Downloading openblas..."
  sudo apt-get install libopenblas-base
  if [ ! -d "./libs/SuiteSparse" ]; then
    echo "SuiteSparse not currently installed - downloading..."
    wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.6.tar.gz -O ./libs/SuiteSparse.tar.gz
    tar xvf ./libs/SuiteSparse.tar.gz -C ./libs
    rm -rf ./libs/SuiteSparse.tar.gz
  fi

  # Check if libopenblas.a exists in /usr/lib, if not check if it is
  # in /usr/lib/openblas-base and create a link.
  if [ ! -e "/usr/lib/libopenblas.a" ]; then
    if [ -r "/usr/lib/openblas-base/libopenblas.a" ]; then
      sudo ln -s /usr/lib/openblas-base/libopenblas.a /usr/lib/;
    fi
  fi

  # Make SuiteSparse
  if [ ! -e "libs/SuiteSparse/CHOLMOD/Lib/libcholmod.a" ]; then
    echo "Making SuiteSparse..."
    cd ./libs/SuiteSparse
    make
    cd ../../
  fi

  # Install lapack
  sudo apt-get -y install libblas-dev checkinstall
  sudo apt-get -y install liblapack-dev checkinstall

  # Check whether liblapacke exists
  if [ ! -e /usr/lib/liblapacke.a ]; then
    echo "LAPACKE not found - donwloading..."
    wget http://www.netlib.org/lapack/lapack-3.6.0.tgz -O ./libs/lapack-3.6.0.tgz

    echo "LAPACKE downloaded to ./libs - now unpacking..."
    tar zxf ./libs/lapack-3.6.0.tgz -C ./libs

    echo "LAPACKE unpacked in ./libs - now removing the tgz file..."
    rm ./libs/lapack-3.6.0.tgz

    echo "Entering libs/lapack..."
    cd ./libs/lapack-3.6.0/
    cp make.inc.example make.inc
    cd LAPACKE

    echo "Making LAPACKE..."
    make

    echo "Copying LAPACKE files into /usr/lib and /usr/include..."
    sudo cp ./include/lapacke.h /usr/include/
    sudo cp ./include/lapacke_mangling.h /usr/include/
    sudo cp ./include/lapacke_utils.h /usr/include/
    cd ..
    sudo cp ./liblapacke.a /usr/lib
    cd ../../
  fi

  make
}

exec 3>&1 1>>${LOG_FILE} 2>&1
install_libforbes | tee /dev/fd/3;

CONFIG_FILE=config.mk
echo "# Auto-generated config file." > $CONFIG_FILE
echo -e "# Path to your SuiteSparse directory" >> $CONFIG_FILE
echo -e "SS_DIR = libs/SuiteSparse" >> $CONFIG_FILE
echo -e "# Extra headers and libraries paths" >> $CONFIG_FILE
echo -e "IEXTRA = /usr/include" >> $CONFIG_FILE
echo -e "LEXTRA = /usr/lib" >> $CONFIG_FILE

