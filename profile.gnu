# calculate the number of points
stats 'test.dat' using 1:2 nooutput

# if you want to have a fixed range for all plots
set xrange [STATS_min_x:STATS_max_x]
set yrange [-1:1]

set terminal pngcairo size 800,400
outtmpl = 'output%07d.png'

do for [i=2:3]{
	set output sprintf(outtmpl, i)
	plot 'test.dat' u 1:i w lp 
	}
set output


