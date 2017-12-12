reset

set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10
set border lw 5 lc rgb "#484848"

stats sprintf("%s.TS.txt", ARG1) usi 1 prefix "A"
stats sprintf("%s.TS.txt", ARG1) usi 2 prefix "B"

minx = (int(A_max/5) - 0   - (500 + int(A_max/5))%500)
maxx = (int(A_max) + 500 + (500 - int(A_max))%500)
set xr [minx:maxx]

miny = (int(B_min) - 0  - (10 + int(B_min))%10)
maxy = (int(B_max) + 10 + (10 - int(B_max))%10)
set yr [miny:maxy]

set xtics minx, (maxx - minx) /4 , maxx tc rgb "#484848"
set ytics miny, (maxy - miny) / 4, maxy tc rgb "#484848"

#set format x "%s"
set format y "%.1f"

set xlabel "Number of Models" tc rgb "#484848" font "Arial-Bold, 57"
set ylabel "Best score" tc rgb "#484848" font "Arial-Bold, 57"
set key tc rgb "#484848"

set linetype 5 dashtype 2 lw 10
set arrow nohead from  minx,B_max to maxx,B_max  lt 5 lc rgb "red"

set output sprintf("%s.TS.pdf", ARG1)
plot sprintf("%s.TS.txt", ARG1) usi 1:2:3 w errorbars lw 2 pt 7 ps 2.5 lc rgb "#484848" notitle
set output

