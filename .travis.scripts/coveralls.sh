#!/bin/bash

echo "...trying coveralls"

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
#if [ "$CXX" == "g++" ];
#then
  sudo pip install cpp-coveralls
  make clean
  ./configure --with-boost=${BOOST_ROOT} LDFLAGS=--coverage CXXFLAGS=--coverage
  make
  cd test
  ./variant_test
  cpp-coveralls -r ../ -e examples -e doxy -e R -e rtdocs --verbose -t ${COVERALLS_TOKEN}
#fi
