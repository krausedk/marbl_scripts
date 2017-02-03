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
	perl get_useful_arrays_v3.pl $1 $1.$STARTLINE.out $STARTLINE $ENDLINE &
	echo "PERL $STARTLINE STARTED"
	STARTLINE=$[$STARTLINE+$LINEVAR+2]
	ENDLINE=$[$ENDLINE+$LINEVAR+2]
done
wait
echo "USEFULARRAYS COMPLETE"
cat $1.*.out > $1_cat.out
echo "FILE_CAT COMPLETE"
perl reconfigure_catfile.pl $1_cat.out $1_cat_reconfigured.out
echo "FILE RECONFIGURATION COMPLETE"
perl marbl_v3.pl $1_cat_reconfigured.out 0.001 0.000001
echo "MARBL COMPLETE!"

