#!/bin/bash

# These are exported so be used in envsubst
export TIMESTEP=0.001 
export NUMSTEPS=2
export IO_CHECKSTEPS=200
export IO_INFOSTEPS=200
export GEOMETRY=s-s5.xmt_txt
export PUMIMESH=s-s5_72_pumi.sms
export ORDER=5
SOLVER=../../builds/dist/bin/IncNavierStokesSolver
FLDTOVTK=../../builds/dist/bin/FieldConvert
MODES=(4) #(4 5 6)
for N in ${MODES[@]}; do 
	export NUMMODES=$N
	# echo "Running $N"
	envsubst < template.xml > m${N}.xml
	$SOLVER m${N}.xml |& tee m${N}.log
	# tail m${N}.log
	# echo " done"
	$FLDTOVTK m${N}.xml m${N}.fld m${N}.vtu
	# $FLDTOVTK m${N}.xml m${N}.chk 
	sed -i 's/nan/10.0/g' m${N}.vtu

done