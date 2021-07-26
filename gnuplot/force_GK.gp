set term eps

bk = 1.0
temp = 0.8269
area = 3.398E+02

set output "force_GK.eps"
file = "../result/data/force_GK_integ.dat"
plot file us 1:(0.25*($2+$3+$4+$5)/(area*bk*temp)) with lines
