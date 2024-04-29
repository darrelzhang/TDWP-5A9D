#!/bin/bash


procs=1

run0()
{
for opt in $@
do
   echo ${opt}  
   sed "s/potopt/\"${opt}\"/" input_sca > input
   sleep 3
   mpirun -hostfile machinefile -n ${procs} ./xynh2.exe >  out.${opt} 2>&1
   sleep 3
done
}

case $1 in
'sca0')
  run0 asy int sca ;;
esac

