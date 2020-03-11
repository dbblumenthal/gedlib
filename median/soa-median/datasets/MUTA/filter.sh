#!/bin/bash

SIZE=$1
NB_REPET=$2
REPER=$3
ALL_FILES=all_files.cxl
#tmp_$$.cxl

#(for file in $REPER/*.cxl
#do
#    cat $file
#done
#)>$ALL_FILES



declare -a files
files=($(find $REPER -type f -name "*.gxl" -exec  grep  -l "node id=\"100\"" {}  -ne 0 2>/dev/null \;))
#     | while read file
#do
#    muta=`grep $file $ALL_FILES | cut -d= -f3`
    
#done
declare -a mutagen
declare -a nonmutagen

for FILE in "${files[@]}"
do
    file=`basename $FILE`
    muta=`grep $file $ALL_FILES | cut -d= -f3|cut -d\" -f2`
#    echo $file $muta
    if [ $muta == "mutagen" ]
    then
	mutagen+=($file)
    else
	nonmutagen+=($file)
    fi
done
#set -x
NB_SAMPLE=`expr  $SIZE / 2`
i=1
nb_muta=${#mutagen[@]}
nb_nonmuta=${#nonmutagen[@]}

while [ $i -le $NB_REPET ]
do
    MUTAGEN_FILE=file_$i.ds
    >$MUTAGEN_FILE
    j=1
    all_muta=""
    all_nonmuta=""
#    echo "<?xml version=\"1.0\"?>\n<GraphCollection>\n<mutagenicity count=\"$SIZE\">" > $MUTAGEN_FILE
    while [ $j -le $NB_SAMPLE ]
    do
	k1=`expr $RANDOM \% $nb_muta`
	muta=${mutagen[$k1]}
	while echo "$all_muta" | grep -q $muta
	do
	    	k1=`expr $RANDOM \% $nb_muta`
		muta=${mutagen[$k1]}
	done
	all_muta="$all_muta $muta"
	#	echo "<print file=\"$muta\" class=\"mutagen\"/>" >> $MUTAGEN_FILE
	echo "$muta 1"  >> $MUTAGEN_FILE

	k2=`expr $RANDOM \% $nb_nonmuta`
	non_muta=${nonmutagen[$k2]}
	while echo "$all_nonmuta" | grep -q $non_muta
	do
	    	k2=`expr $RANDOM \% $nb_nonmuta`
		non_muta=${nonmutagen[$k2]}
	done
	all_nonmuta="$all_nonmuta $non_muta"
	#	echo "<print file=\"$non_muta\" class=\"nonmutagen\"/>" >> $MUTAGEN_FILE
	echo "$non_muta 0" >> $MUTAGEN_FILE
	j=`expr $j + 1`
    done
#    echo "</mutagenicity></GraphCollection>">> $MUTAGEN_FILE
    i=`expr $i + 1`
done
#rm -f $ALL_FILES
