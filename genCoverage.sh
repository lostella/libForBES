#!/bin/bash
make clean
#make
#make build-tests
make -j 8 test
lcov --directory ./build/Debug --capture --output-file coverage.info 
lcov --remove coverage.info '/usr/include/*' 'main.cpp' --output-file libforbes-coverage.info
genhtml -s --legend --title 'LibForBES Unit Tests' libforbes-coverage.info --output-directory coverage
