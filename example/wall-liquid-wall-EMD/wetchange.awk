BEGIN{
    wet = 7.7382E-01
    eta = 0.4
}
FNR!=9 && FNR!=10{
    print $0 > FILENAME
}
FNR==9 || FNR==10{
    print $1, $2, wet*eta, $4, $5 > FILENAME
}
