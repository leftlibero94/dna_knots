#!/bin/zsh

cd 99bp_loose_LJconf

for i in seed*
do
       cd $i
       gfortran knot99bp_LJconf*.f -o knot99bp_LJconf$i.out
       ./knot99bp_LJconf$i.out &
       cd ..
done
