#!/bin/bash
#make clean
#make
#make build-tests
make -j 8 test
lcov --capture --directory ./build/Debug --output-file coverage.info 
lcov --remove coverage.info '/usr/include/*' 'main.cpp' --output-file libforbes-coverage.info
genhtml libforbes-coverage.info --output-directory coverage
