#!/bin/bash
#SBATCH -p gpu         # partition name
#SBATCH -n 1                # number of core
#SBATCH --gres=gpu:1        # request 1 gpu
#SBATCH -J crf      # jobname
#SBATCH -o cuda.out         # standard output
#SBATCH -e cuda.err         # standart error

nstlim=500000000
ntw=10000
ntraj=`ls traj_*md*.* | wc -l`
if [ $ntraj -eq 0 ]; then
    symbol=!
else
    symbol=
fi
cat << EOF > md.in
md script
&cntrl
${symbol}irest = 1, ntx = 5,
imin = 0, ntb = 2, ntp = 1,
!igb = 5,
ntpr = $ntw, ntwr = $ntw, ntwx = $ntw,
ntt = 3, gamma_ln = 1.0, tempi = 300.0, temp0 = 300.0,
ig = -1,
ntf = 2, ntc = 2,
iwrap = 1,
nstlim = $nstlim, dt = 0.002,
cut =12.0,
ioutfm = 1,ntwprt=464,
/
EOF
if [ $ntraj -eq 0 ]; then
    srun --gres=gpu:1 pmemd.cuda -O -i md.in -o md1.out -p sys.top -c equi.rst -r md1.rst -x traj_md1.nc -ref equi.rst
else
    let i=ntraj
    let j=ntraj+1
    srun --gres=gpu:1 pmemd.cuda -O -i md.in -o md${j}.out -p sys.top -c md${i}.rst -r md${j}.rst -x traj_md${j}.nc -ref md${i}.rst
fi
