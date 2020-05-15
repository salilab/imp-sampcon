reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10
set border lw 5 lc rgb "#484848"

stats sprintf("%s.ChiSquare_Grid_Stats.txt", sysname) usi 1 prefix "R"
stats sprintf("%s.Sampling_Precision_Stats.txt", sysname) usi 1 prefix "B"

maxx = (int(R_max) - 0  - (10 + int(R_max))%10)
max(a, b) = (a > b ? a : b)
max0=max(maxx, B_max)

set xr [0:max0]
set yr [0:1.1]

set xtics 0, maxx / 4, maxx tc rgb "#484848"
set ytics 0,0.25,1 tc rgb "#484848"
set format y "%.2f"

set encoding iso_8859_1 
set xlabel "Threshold ({\305})" tc rgb "#484848" font "Arial-Bold, 57"
set ylabel "Convergence Criteria" tc rgb "#484848" font "Arial-Bold, 57"

set key above width -2 vertical maxrows 3 tc rgb "#484848"
set linetype 5 dashtype 2 lw 10

set arrow nohead from 0,0.05 to maxx,0.05 lt 5 lw 10 lc rgb "#FF4500" back filled
set arrow nohead from 0,0.10 to maxx,0.10 lt 5 lw 10 lc rgb "#5B6FB5" back filled
set arrow nohead from 0,0.80 to maxx,0.80 lt 5 lw 10 lc rgb "#61B329" back filled

set arrow nohead from B_min,0 to B_max,1.1 lt 5 lw 10 lc rgb "#484848" back filled

set output sprintf("%s.ChiSquare.pdf" , sysname)
plot sprintf("%s.ChiSquare_Grid_Stats.txt", sysname) usi 1:2 w p pt 7 ps 2.5 lc rgb "#FF4500" title "{/Symbol c}^2-test p-value", \
     "" 	            	usi 1:3 w p pt 5 ps 2.5 lc rgb "#5B6FB5" title "Cramer's V", \
     ""				usi 1:($4/100) w p pt 9 ps 2.5 lc rgb "#61B329" title "Clustered population"

set output
