#!/bin/bash
# Convert text to FASTA file
# Assumes no header

###----- Parameters
IN_FILE='M13mp7_ref.txt'
OUT_FILE='M13mp7_ref.fasta'
REF_HEAD='>M13mp7_ref'
LINE_LEN=74
###-----

touch $OUT_FILE
echo $REF_HEAD >> $OUT_FILE
TXT_LEN=$(wc -c $IN_FILE | cut -d " " -f1)
SEQ=$(cat $IN_FILE)
NEWLINE=$'\n'

for (( i=1; i<$TXT_LEN+1; i++ ));
do
  if [ $i != 1 ] && [ $(( $i%$LINE_LEN )) = 1 ]; then
    echo $NEWLINE >> $OUT_FILE
  fi
  echo -n ${SEQ:$i:1} >> $OUT_FILE
done

echo 'Done!'
