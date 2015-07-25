#!/bin/bash
make clean
make
make build-tests
make test
lcov --capture --directory ./build/Debug/GNU-Linux-x86 --output-file coverage.info 
lcov --remove coverage.info '/usr/include/*' 'main.cpp' 'Function.*' 'LinearOperator.cpp' '*.h' --output-file libforbes-coverage.info
genhtml libforbes-coverage.info --output-directory coverage
