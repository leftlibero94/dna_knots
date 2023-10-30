#!/bin/zsh
mkdir 99bp_loose_LJfree
cd 99bp_loose_LJfree

for i in {1..10}
do
       mkdir seed$i
       cd seed$i
       rand=$(od -N 4 -t dI -An /dev/urandom | tr -d " ")
       echo $rand > seed.inp
# CHANGE THE FOLLOWING LINE
       cp /Users/souradeep/Coding/dnaknots/knot_sim/{trefoil_short_loose.txt,knot99bp_LJfree.f,test1.inp} .
       cd ..
done
