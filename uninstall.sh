#! /bin/sh
# This is a very simplistic uninstall scipt. Use with care!

if [ -z $1 ] && [ "$1" = "--system" ]; then
  rm -vr /usr/share/molmod
  rm -vr /usr/lib/python*/site-packages/molmod
else
  rm -vr $HOME/share/molmod
  rm -vr $HOME/lib/molmod
fi
