#!/bin/zsh

cd 99bp_loose_LJfree

for i in seed*
do
       cd $i
       gfortran knot99bp_LJfree*.f -o knot99bp_LJfree$i.out
       ./knot99bp_LJfree$i.out &
       cd ..
done
