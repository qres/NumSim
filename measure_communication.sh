NUM_P=(1 2 4 8)
GRID=$1
VALS=(t_total t_comm t_red t_black t_res blocked_t_res)
OUT=out_sweep

cp default.geom "measure-"$OUT$GRID".geom"
echo "size $GRID $GRID" >> "measure-"$OUT$GRID".geom"

# clear files
mkdir $OUT
for v in ${VALS[@]}; do
    echo "# p "$v"_min "$v"_avg "$v"_max" > $OUT/$GRID"_$v.dat"
done

# communication times
for p in ${NUM_P[@]}; do
    # run measurement with #processes
    salloc --hint=compute_bound -m block:block -n $p mpirun ./build/numsim default.param "measure-"$OUT$GRID".geom" &> "_measure-"$OUT$GRID".out"
    # append result to files
    for v in ${VALS[@]}; do
        dat=$OUT/$GRID"_$v.dat"
        echo -n "$p " >> $dat
        grep "^$v," "_measure-"$OUT$GRID".out" | cut -d, -f 2 | tr -d '\n' >> $dat; echo -n " " >> $dat
        grep "^$v," "_measure-"$OUT$GRID".out" | cut -d, -f 3 | tr -d '\n' >> $dat; echo -n " " >> $dat
        grep "^$v," "_measure-"$OUT$GRID".out" | cut -d, -f 4              >> $dat
    done
done

#squeue
#scancel -u huberfx
