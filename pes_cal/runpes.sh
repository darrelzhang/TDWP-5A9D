#!/bin/bash

exefile=xynh2_pes.exe

potnum=19

for opt in asy int
do
   echo ${opt}	
   sed "s/potopt/\"${opt}\"/" input_pot0 > input_pot 

  for inode in $(seq 0 ${potnum})
  do
    echo ${inode}
    sed "s/abcde/${potnum} ${inode} /" input_pot > input 
    ./${exefile}  >  result/out${inode} &
    sleep 2 
  done
  wait
  sleep 2
done

