#! /bin/sh
# This is a very simplistic uninstall scipt. Use with care!

./uninstall.sh
python setup.py install
./cleanfiles.sh
