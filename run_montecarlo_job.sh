#!/bin/bash

# run like sbatch --array=0-188 -J "Jobname" run_montecarlo.sh

function run_job {
    ./build/numsim montecarlo.param montecarlo.geom | sed -n "s#S_T_E_P: \(.*\)#\1#p" > "mc/$((${SLURM_ARRAY_TASK_ID} * 8 + $1)).dat.tmp"
    mv "mc/$((${SLURM_ARRAY_TASK_ID} * 8 + $1)).dat.tmp" "mc/$((${SLURM_ARRAY_TASK_ID} * 8 + $1)).dat"
}

function run_tr_job {
    mu=1500
    sigma="(1000/6)"
    ix=$((${SLURM_ARRAY_TASK_ID} * 8 + $1))
    cp montecarlo.param "montecarlo.param.$((${SLURM_ARRAY_TASK_ID} * 8 + $1))"
    echo "sigma_re = 0"  >> "montecarlo.param.$((${SLURM_ARRAY_TASK_ID} * 8 + $1))"
    python3 -c "if True:
        import numpy as np

        def hier_ix(ix, shift_n1p1=True):
            level = int(np.log2(ix+1))
            offset = 0.5**(level + 1)
            stride = 0.5**level
            locix = ix - (2**level - 1)

            if shift_n1p1:
                return (offset + locix*stride)*2 - 1
            else:
                return offset + locix*stride

        if $ix == 0:
            print('re = {}'.format($mu - 3 * $sigma))
        elif $ix == 1:
            print('re = {}'.format($mu + 3 * $sigma))
        else:
            print('re = {}'.format($mu + hier_ix($ix - 2, True) * 3 * $sigma))
    " >> "montecarlo.param.$((${SLURM_ARRAY_TASK_ID} * 8 + $1))"

    ./build/numsim "montecarlo.param.$((${SLURM_ARRAY_TASK_ID} * 8 + $1))" montecarlo.geom | sed -n "s#S_T_E_P: \(.*\)#\1#p" > "mc/$((${SLURM_ARRAY_TASK_ID} * 8 + $1)).dat.tmp"
    mv "mc/$((${SLURM_ARRAY_TASK_ID} * 8 + $1)).dat.tmp" "mc/$((${SLURM_ARRAY_TASK_ID} * 8 + $1)).dat"
}

for i in $(seq 0 7); do
    run_job $i &
    sleep 1
done

wait
