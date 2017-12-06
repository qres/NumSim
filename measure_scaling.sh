NUM_P=(1 2 4 8 16 32 64)
GRID=$1
OUT=out_sweep

cp default.geom "measure-"$OUT$GRID".geom"
echo "size $GRID $GRID" >> "measure-"$OUT$GRID".geom"

# clear files
v=runtime
mkdir $OUT
echo "# p "$v"_min "$v"_avg "$v"_max" > $OUT/$GRID"_$v.dat"

# communication times
for p in ${NUM_P[@]}; do
    # run measurement with #processes
    salloc --hint=compute_bound -m block:block -n $p mpirun ./build/numsim default.param "measure-"$OUT$GRID".geom" &> "_measure-"$OUT$GRID".out"
    # append result to files
    dat=$OUT/$GRID"_$v.dat"
    echo -n "$p " >> $dat
    grep "^$v," "_measure-"$OUT$GRID".out" | cut -d, -f 2 >> $dat;
done

#squeue
#scancel -u huberfx
