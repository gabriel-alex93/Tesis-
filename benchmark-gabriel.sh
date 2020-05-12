#!/bin/sh
if [ "$#" -ne 8 ]; then
    # AQUI PUEDE CAMBIAR 
    echo "run as ./benchmark-gabriel.sh BINARY   REP NC BS STARTN ENDN DN SAMPLES"
    exit;
fi
# AQUI PUEDE CAMBIAR
BINARY=$1 
REP=$2
NC=$3
BS=$4
STARTN=$5
ENDN=$6
DN=$7
SAMPLES=$8
echo ${SAMPLES}
for N in `seq ${STARTN} ${DN} ${ENDN}`;
do
    M=0
    S=0
    # AQUI PUEDE CAMBIAR
    echo "./${BINARY} ${N} 1 ${REP} ${NC} ${BS}"
    for k in `seq 1 ${SAMPLES}`;
    do
        # AQUI PUEDE CAMBIAR
        x=`./${BINARY} ${N} 1 ${REP} ${NC} ${BS}`
        oldM=$M;
        M=$(echo "scale=20;  $M+($x-$M)/$k"           | bc)
        S=$(echo "scale=20;  $S+($x-$M)*($x-${oldM})" | bc)
    done
    echo "done"
    MEAN=$M
    VAR=$(echo "scale=20; $S/(${SAMPLES}-1.0)"  | bc)
    STDEV=$(echo "scale=20; sqrt(${VAR})"       | bc)
    STERR=$(echo "scale=20; ${STDEV}/sqrt(${SAMPLES})" | bc)
    TMEAN1=${MEAN}
    TVAR1=${VAR}
    TSTDEV1=${STDEV}
    TSTERR1=${STERR}

    echo "--->MODE 1: N=${N}   --> (MEAN, VAR, STDEV, STERR) -> (${TMEAN1}[ms], ${TVAR1}, ${TSTDEV1}, ${TSTERR1})"

    M=0
    S=0
    # AQUI PUEDE CAMBIAR
    echo "./${BINARY} ${N} 2 ${REP} ${NC} ${BS}"
    for k in `seq 1 ${SAMPLES}`;
    do
        # AQUI PUEDE CAMBIAR
        x=`./${BINARY} ${N} 2 ${REP} ${NC} ${BS}`
        oldM=$M;
        M=$(echo "scale=20;  $M+($x-$M)/$k"           | bc)
        S=$(echo "scale=20;  $S+($x-$M)*($x-${oldM})" | bc)
    done
    echo "done"
    MEAN=$M
    VAR=$(echo "scale=20; $S/(${SAMPLES}-1.0)"  | bc)
    STDEV=$(echo "scale=20; sqrt(${VAR})"       | bc)
    STERR=$(echo "scale=20; ${STDEV}/sqrt(${SAMPLES})" | bc)
    TMEAN2=${MEAN}
    TVAR2=${VAR}
    TSTDEV2=${STDEV}
    TSTERR2=${STERR}
    # AQUI PUEDE CAMBIAR
    echo "--->MODE 2: N=${N}   --> (MEAN, VAR, STDEV, STERR) -> (${TMEAN2}[ms], ${TVAR2}, ${TSTDEV2}, ${TSTERR2})"
    # AQUI PUEDE CAMBIAR
    echo "${N} ,   ${TMEAN1} , ${TVAR1}  ,${TSTDEV1}  , ${TSTERR1}       ,         ${TMEAN2}     ,   ${TVAR2} , ${TSTDEV2} , ${TSTERR2}" >> data/running-times_NC${NC}_BS${BS}.dat
    echo " "
done 
