#!/bin/bash
##  start the timer
SECONDS=0

# change working directory
cd ~/CHIKV_DEG

##  Total runtime
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."