#!/bin/bash
#SBATCH -p gpu         # partition name
#SBATCH -n 1                # number of core
#SBATCH --gres=gpu:1        # request 1 gpu
#SBATCH -J crf      # jobname
#SBATCH -o cuda.out         # standard output
#SBATCH -e cuda.err         # standart error

nstlim=500000
ntw=10000
nres=12

cat << EOF > md.in
md script
&cntrl
!irest = 1, ntx = 5,
ntb = 2,ntp = 1,
ntpr = $ntw, ntwr = $ntw, ntwx = $ntw,
ntt = 3, gamma_ln = 1.0, 
tempi = 300.0, temp0 = 300.0,
ntf = 2, ntc = 2,
iwrap=1,ntwprt=464,
nstlim = $nstlim, dt = 0.002,
cut =12.0,
ntr=1,
itunnel=2,staysteps=$nstlim,
xminpos=14.5,
/
Keep protein fixed with weak restraints
10.0
RES 2 $nres
END
END
EOF

srun --gres=gpu:1 pmemd.cuda -O -i md.in -o md1.out -p sys.top -c equi.rst -r md1.rst -x traj_md1.nc -ref equi.rst

for i in `seq 3 $nres`
do
    cat << EOF > md.in
md script
&cntrl
irest = 1, ntx = 5,
ntb = 2,ntp = 1,
ntpr = $ntw, ntwr = $ntw, ntwx = $ntw,
ntt = 3, gamma_ln = 1.0, 
tempi = 300.0, temp0 = 300.0,
ntf = 2, ntc = 2,
iwrap=1,ntwprt=464,
nstlim = $nstlim, dt = 0.002,
cut =12.0,
ntr=1,
itunnel=2,staysteps=$nstlim,
xminpos=14.5,
/
Keep protein fixed with weak restraints
10.0
RES $i $nres
END
END
EOF
    let j=i-1
    let k=i-2
    srun --gres=gpu:1 pmemd.cuda -O -i md.in -o md${j}.out -p sys.top -c md${k}.rst -r md${j}.rst -x traj_md${j}.nc -ref md${k}.rst
done
