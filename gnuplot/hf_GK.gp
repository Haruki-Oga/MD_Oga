set term eps

bk = 1.0
temp = 0.8269
area = 3.398E+02

set output "hf_GK.eps"
file = "../result/data/hf_GK_integ.dat"
set logscale x
plot file us 1:(0.5*($2+$3)/(area*bk*temp*temp)) with lines
