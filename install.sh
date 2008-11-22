#! /bin/sh
# This is a very simplistic install scipt. Use with care!

./uninstall.sh
(cd ext; python setup.py install)
python setup.py install
./clean.sh
