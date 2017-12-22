set style fill transparent solid 0.25 noborder

set xlabel "h"
set ylabel "vortex length"
unset logscale xy


# Runtime vs Communication
set term png
set output "plots/vortex_length.png"
file_vort = "measurements_step.dat"
set title "Vortex length over Discretization"
f(x) = a*x + b
fit [:] f(x) file_vort using (5.0/$1):($2/21.5*5) via a,b
plot file_vort u (5.0/$1):($2/21.5*5) w l lc rgb "blue" t "vortex", \
f(x) lc rgb "red" t sprintf("fit with a = %1.2f",a)

set output "plots/vortex_length_logxy.png"
set logscale xy
fit [:] f(x) file_vort using (log(5.0/$1)):(log($2/21.5*5)) via a,b
plot file_vort u (5.0/$1):($2/21.5*5) w l lc rgb "blue" t "vortex", \
exp(f(log(x))) lc rgb "red" t sprintf("fit with a = %1.2f",a)
#replot
unset logscale xy