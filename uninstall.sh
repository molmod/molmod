#! /bin/bash
# This is a very simplistic uninstall scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  rm -vr /usr/local/share/molmod
  rm -vr /usr/local/lib/python*/site-packages/molmod
else
  rm -vr $HOME/share/molmod
  rm -vr $HOME/lib/python/molmod
fi
