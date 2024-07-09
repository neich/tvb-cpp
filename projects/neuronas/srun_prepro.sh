#!/bin/bash -x

#SBATCH --nodes=6
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --output=c-out.%j
#SBATCH --error=c-err.%j
#SBATCH --partition=debug
#SBATCH --distribution=cyclic:*:*
#SBATCH --hint=nomultithread

cd /home/cluster/imartin/tvb-cpp/projects/cmake-build-relwithdebinfo

./prepro_cluster --data-path /home/cluster/imartin/neuronas/Neuronas/Datos/Datasets/DataHCP80 --out-path /home/cluster/imartin/neuronas/Neuronas/Datos/Results/Results_cluster/tvb-cpp