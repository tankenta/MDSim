#!/bin/bash
# Install Eigen
cd ./third_party
wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz
mkdir eigen && tar xzvf 3.3.4.tar.gz -C eigen --strip-components 1
rm -f 3.3.4.tar.gz

git clone https://github.com/g-truc/glm.git
cd ../
