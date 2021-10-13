BEGIN{
    n = -1
    posz = 0.40130E+002
}
n!=NR{
    print $0
}
$0~"fix_set_pos"{
    n = NR
}
NR==n+1{
    print "6 001 0.0 0.0" posz
}
