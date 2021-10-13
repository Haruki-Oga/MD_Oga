cd $(dirname $0)
NCOR=5
EXE="${HOME}/MD.out"

langevin=(
    "langevin"
    "langevin2"
    "langevin3"
)
    
for i in ${langevin[@]}
do
    echo ${i}
    mpiexec -n ${NCOR} ${EXE} ./calc_${i}.dat
    mpiexec -n ${NCOR} ${EXE} ./calc_${i}.dat
    mkdir -p result_${i}
    ls ./ | grep -v result | xargs -i cp {} result_${i}
done
