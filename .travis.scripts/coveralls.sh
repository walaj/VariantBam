#!/bin/bash

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
if [ "$CXX" == "g++" ];
then
  sudo pip install cpp-coveralls
  cd build
  coveralls -r ../ -e CMakeFiles -e contrib -e test -t ${COVERALLS_TOKEN}
fi
