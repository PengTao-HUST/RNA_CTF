# RNA_CTF
A simplified all-atom molecular dynamics framework to simulate RNA co-transcriptional folding (CTF).

## Usage
Step1: Just put the Fortran code in [code](https://github.com/PengTao-HUST/CotranslationalProteinFoldingSimulations/tree/master/code) folder into the pmemd source code folder in Amber18 ($AMBERHOME/src/pmemd/src), then recompile pmemd.

Step2: Perform CTF simulations by runing the command `sbatch ctf.sh` in your computing cluster.
