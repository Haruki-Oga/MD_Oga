BEGIN{
    set_temp = 0.8269
}
{
    print $4/set_temp,$5/set_temp,$6/set_temp
}
