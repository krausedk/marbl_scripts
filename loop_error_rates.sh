#!/bin/bash
#DIRECTIONS#
#this shell needs get_useful_arrays.pl and marbl_v2.0.pl and dependencies M1627_error_likelihoods.txt and the proper snp_database (needs to be referred to in get_useful_arrays.pl)#
#first variable is name of .sam file, second variable is number of CPUs to use#


LENGTH=$(wc -l $1 | awk '{print $1}')
CPUS=$2
LINESPERCPU=$[$LENGTH/$CPUS]
LINESPERCPUROUNDUP=$(echo "$LINESPERCPU" | awk '{printf("%d\n",$1 + 0.5)}')
if [ $(( $LINESPERCPUROUNDUP % 2)) -eq 0 ]
then
	LINEVAR=$[$LINESPERCPUROUNDUP]
else
	LINEVAR=$[$LINESPERCPUROUNDUP+1]
fi
STARTLINE=4
ENDLINE=$[$STARTLINE+$LINEVAR]
while [ $STARTLINE -lt $LENGTH ]
do
	perl gatk_base_recalibration_emulator_v2.pl $1 $1.$STARTLINE.gatk_errorrates.out $STARTLINE $ENDLINE &
	echo "PERL $STARTLINE STARTED"
	STARTLINE=$[$STARTLINE+$LINEVAR+2]
	ENDLINE=$[$ENDLINE+$LINEVAR+2]
done
wait
echo "measure_error_rates_complete"
