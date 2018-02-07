#!/bin/bash
# Install Eigen
wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz
mkdir lib_eigen && tar xzvf 3.3.4.tar.gz -C lib_eigen --strip-components 1
rm -f 3.3.4.tar.gz
