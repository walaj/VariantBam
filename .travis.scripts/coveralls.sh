#!/bin/bash

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
if [ "$CXX" == "g++" ];
then
  sudo pip install cpp-coveralls
  cd test
  ./variant_test
  cpp-coveralls --verbose
fi
