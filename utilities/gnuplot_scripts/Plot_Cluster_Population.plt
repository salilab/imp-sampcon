reset
set terminal pdfcairo enhanced color font "Arial-Bold,40" size 10,10
set border lw 5 lc rgb "#484848"
set boxwidth 0.5 absolute
set style fill solid 1.00 border rgb "#484848"
unset key 

set encoding iso_8859_1
unset xlabel
set ylabel "Population" tc rgb "#484848" offset 0,0 font "Arial-Bold, 57"
#set y2label "Population (%)" tc rgb "#484848" offset 0,0 font "Arial-Bold, 57"

set pointsize 2
set xtics border in scale 0,0 mirror norotate offset 0.25 left
set xtics  norangelimit
set xtics   ()

stats sprintf("%s.Cluster_Population.txt", sysname) usi 1 prefix "A"
minx = 0
maxx = A_max + 1.75
set xr [minx:maxx]

max(a, b) = (a > b ? a : b)
min(a, b) = (a <= b ? a : b)

stats sprintf("%s.Cluster_Population.txt", sysname) usi 2 prefix "N1"
stats sprintf("%s.Cluster_Population.txt", sysname) usi 3 prefix "N2"

maxy0=max(N1_max, N2_max)
div = ceil(maxy0/2500)
mod=int(maxy0)%2500

ymax=(div)*2500


set yr [0:ymax]
set ytics  0, ymax / 4, ymax nomirror


stats sprintf("%s.Sampling_Precision_Stats.txt", sysname) usi 2 prefix "B2"
stats sprintf("%s.Sampling_Precision_Stats.txt", sysname) usi 3 prefix "B3"

set label 1 sprintf("Number of models:%i, %i\n", N1_sum, N2_sum) at graph 0.95, 0.96 right font 'Arial-Bold, 40' front tc rgb "#484848"
set label 2 sprintf("%i  %i\n", N1_sum, N2_sum) at graph 0.95,0.96 right font 'Arial-Bold, 40.02' front tc rgb "#EB7262"
set label 4 sprintf("%i\n",N2_sum) at graph 0.95,0.96 right font 'Arial-Bold, 40.05' front tc rgb "#000080"

set label 5 sprintf("{/Symbol c}^2-test p-value: %.2f\nCramer's V: %.2f", B2_max, B3_max) at graph 0.95, 0.925 right font 'Arial-Bold, 40' front tc rgb "#484848"

set output sprintf("%s.Cluster_Population.pdf", sysname)
plot sprintf("%s.Cluster_Population.txt", sysname) usi ($1+0.5):2:xtic(sprintf("Cluster %i", $1)) w boxes notitle lc rgb "#EB7262", \
     "" 		       usi ($1+1):3 w boxes notitle lc "#000080"

set output

