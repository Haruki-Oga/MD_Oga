set term eps

set output "temp_profile.eps"

file_liquid = "../result/data/chunk1_3.dat"
file_wall = "../result/data/awk.log"

set x2tics
set xtics nomirror
set yrange [0:45]
set x2range [0.7:1.0]
set arrow 1 from second 0.8269, 0.0 to second 0.8269, 45 nohead

plot file_liquid us 2:1 with lines axis x1y1,\
     file_liquid us 3:1 with lines axis x2y1,\
     file_wall us 3:1 axis x2y1

