load '~/gnuplot/default.cfg'

set terminal png
set output "plots/err_h.png"
set xlabel "h_y"
set title "Error from analytic solution"
plot 'out_err/discr.dat' u 2:3 w lp t "Error ||.||_2"
