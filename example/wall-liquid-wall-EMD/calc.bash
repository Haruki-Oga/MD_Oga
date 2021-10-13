#!/bin/bash

EXE="${HOME}/MD.out"
RUN="mpiexec"
COR="10"

${RUN} -n ${COR} ${EXE} ./calc_pressctrl_relax.dat
${RUN} -n ${COR} ${EXE} ./calc_pressctrl.dat
cat calc_setpos0.dat | sed -e "s/setposz/`awk '{print $6}' ./ave_compv.dat`/g" > calc_setpos.dat
${RUN} -n ${COR} ${EXE} ./calc_setpos.dat
${RUN} -n ${COR} ${EXE} ./calc_relax.dat
${RUN} -n ${COR} ${EXE} ./calc.dat
