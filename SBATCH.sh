#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="startup"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu


module load cmake/3.12.3
module load cuda/9.1
module load extra
module load GCCcore/6.3.0


srun -p gpu --gres=gpu:1 ./build/bend-model reset src/Linux-Debug/data_Cir_0.0249_30_30_30_0.7.xml build/State__domain_29_plt_count_90_filo_count_10_maxForce_29.00_minForce_2.10_pltForceScale_0_dynamicResp_1_nonlinDynamicResp_0/State_13.000000.sta
