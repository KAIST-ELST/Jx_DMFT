#!/bin/bash
#
#
#$ -l h_rt=72:00:0
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -o std.out2
#$ -e std.err
#$ -pe orte 56 -l h=node05




## source ~/.profile_hk
echo $job_id
echo $pe
echo $pe_slots
echo $NSLOTS
echo $pe_hostfile
echo "name:$name name2:$name2 nt:$nt  "
echo WIENROOT = $WIENROOT
echo WIEN_DMFT_ROOT = $WIEN_DMFT_ROOT
echo COMSUITE_BIN = $COMSUITE_BIN
env

node_list=$TMPDIR/machines
/home/users1/xorwndekda/intel/oneapi/mpi/2021.1.1/bin/mpirun -np $NSLOTS hostname > $node_list

echo "Will run command: julia_rsh -np $NSLOTS -machinefile $TMPDIR/machines"
cat $node_list



##~/opt/local/mvapich2/2.1/bin/mpirun -np $NSLOTS  /home/users1/bluehope/work_local/DFTcodes/openmx/openmx/source/openmx  NiO.dat |tee NiO.log

#echo "mpirun -np $NSLOTS" > mpi_prefix.dat
#echo "mpirun -np $NSLOTS" > mpi_prefix.dat2

echo $JOB_ID >> running
#mpirun -np $NSLOTS   ~/comsuite/COMSUITE/bin/rspflapw.exe  |tee -a comlog
#~/opt/local/mvapich2/2.3/bin/mpirun -np $NSLOTS   ~/comsuite/COMSUITE/bin/rspflapw.exe  |tee -a comlog
#python /home/users1/xorwndekda/comsuite/comsuite-master/bin/comdmft.py >&screen.out
#julia -p $NSLOTS Jx_modelDMFT.jl -T nio_J_wannier.toml |& tee screen.out
julia --machine-file $node_list Jx_DMFT.jl -T fe_J_wannier.toml |& tee screen.out

date >> done

##echo $sge_o_workdir >>  ~/work/done.mj2 
#echo $sge_o_workdir >>  ~/work/done.mj2 
echo $pwd >>  ~/work/done.mj2 
