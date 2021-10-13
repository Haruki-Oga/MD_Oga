#!/bin/bash
awk -f temp_set.awk < result_langevin/avetemp.dat
awk -f temp_set.awk < result_langevin2/avetemp.dat
awk -f temp_set.awk < result_langevin3/avetemp.dat
