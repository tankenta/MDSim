#!/bin/bash

# c.f. glfw: Better ubuntu (apt-get) install instructions #808
# https://github.com/glfw/glfw/issues/808

sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 3EB9326A7BF6DFCD
sudo apt-get install software-properties-common
sudo add-apt-repository 'deb http://ppa.launchpad.net/keithw/glfw3/ubuntu trusty main'
sudo apt-get update
sudo apt-get install libglfw3 libglfw3-dev
