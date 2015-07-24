#!/bin/bash
make clean
make
make build-tests
make test
lcov --capture --directory ./build/Debug/GNU-Linux-x86 --output-file coverage.info
genhtml coverage.info --output-directory coverage
