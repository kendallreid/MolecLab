set key autotitle columnhead
set xlabel 'Time'
set ylabel 'Concentration'
set datafile separator ','
plot for [col=2:*] 'concentrations.csv' using 1:col with lines lw 2
