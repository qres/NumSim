set style fill transparent solid 0.25 noborder

set xlabel "#processes"
set ylabel "time [ns]"
unset logscale xy


# Runtime vs Communication
set term png
set output "plots/4096_total_vs_comm.png"
file_tot = "out_sweep/4096_t_total.dat"
file_com = "out_sweep/4096_t_comm.dat"
set title "Total vs Communication"
plot file_tot u 1:2:4 w filledcu lc rgb "blue" t "total min-max", file_tot u 1:3 w l lc rgb "blue" t "total avg", \
     file_com u 1:2:4 w filledcu lc rgb "red"  t "comm min-max" , file_com u 1:3 w l lc rgb "red"  t "comm avg"

set output "plots/4096_total_vs_comm_logxy.png"
set logscale xy
replot
unset logscale xy

set title "Total"
set output "plots/4096_total.png"
plot file_tot u 1:2:4 w filledcu lc rgb "blue" t "total min-max", file_tot u 1:3 w l lc rgb "blue" t "total avg"        # filled curve

set output "plots/4096_total_logxy.png"
set logscale xy
f(x) = a*x + b
fit [3:] f(x) file_tot using (log($1)):(log($3)) via a,b                                                                # log fit to 8... processes
plot file_tot u 1:2:4 w filledcu lc rgb "blue" t "total min-max", file_tot u 1:3 w l lc rgb "blue" t "total avg",  \
     exp(f(log(x))) lc rgb "red" t sprintf("fit to 8..64 with a = %1.2f",a)                                             # casting to string
unset logscale xy

set title "Communication"
set output "plots/4096_comm.png"
plot file_com u 1:2:4 w filledcu lc rgb "red" t "comm min-max", file_com u 1:3 w l lc rgb "red" t "comm avg"

set output "plots/4096_comm_logxy.png"
set logscale xy
replot
unset logscale xy

# Runtime of Residual Computation
set term png
set output "plots/4096_comm_vs_res.png"
file_res = "out_sweep/4096_t_res.dat"
file_com = "out_sweep/4096_t_comm.dat"
set title "Communication vs Res"
set yrange [1:5000000]
plot file_res u 1:2:4 w filledcu lc rgb "green" t "res min-max" , file_res u 1:3 w l lc rgb "green" t "res avg", \
     file_com u 1:2:4 w filledcu lc rgb "red"   t "comm min-max", file_com u 1:3 w l lc rgb "red"   t "comm avg"
set autoscale y

set output "plots/4096_comm_vs_res_logxy.png"
set logscale xy
replot
unset logscale xy

set title "Reduction"
set output "plots/4096_res.png"
plot file_res u 1:2:4 w filledcu lc rgb "green" t "res min-max", file_res u 1:3 w l lc rgb "green" t "res avg"

set output "plots/4096_res_logxy.png"
set logscale xy
replot
unset logscale xy

# runtime
set title "Runtime for 3 Iters with 2 SOR steps each"
set output "plots/4096_runtime.png"
file_tot = "out_sweep/4096_runtime.dat"
fit [3:] f(x) file_tot using (log($1)):(log($2)) via a,b
plot file_tot u 1:2 w l lc rgb "blue" t "runtime", \
     exp(f(log(x))) lc rgb "red" t sprintf("fit to 8..64 with a = %1.2f",a)                                             # casting to string

set output "plots/4096_runtime_logxy.png"
set logscale xy
replot
unset logscale xy

#speedup
set term png
set output "plots/4096_speedup.png"
set title "Speedup"
plot "out_sweep/4096_speedup.dat" u 1:($3/$2) lc rgb "red" t "speedup"

set term png
set output "plots/theo_speedup.png"
set title "Speedup for 75% parallelizable code"
set ylabel "Speedup"
set xrange [1:100]
plot 1/(1-0.75 + 0.75/x)

set term png
set output "plots/4096_par_eff.png"
set title "Parallel Efficiency"
set ylabel "Efficiency"
set xrange [1:100]
plot "out_sweep/4096_speedup.dat" u 1:($3/$2/$1) lc rgb "red" t "efficiency"
