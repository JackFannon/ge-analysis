#!/bin/bash

rm ge_fit_summary.txt
rm diff_default_smeared.txt

RUNLIST="runlist_for_old_ge_2021.txt"
DATA_DIR="/disk02/usr6/sshima/Ge2021/data/LINAC2021/"
OUTPUT_DIR="fit_results/"

cat $RUNLIST | while read line; do
    echo $line | grep -v '^#.*' > /dev/null
    if [ $? -ne 0 ]; then
	continue
    fi
    RUNNO=`echo $line | cut -f1 -d ' '`
    ENERGY=`echo $line | cut -f2 -d ' '`
    XPOS=`echo $line | cut -f4 -d ' '`
    ZPOS=`echo $line | cut -f5 -d ' '`
    DATA_FILE=`echo $line | cut -f6 -d ' '`

    OUTPUT_FILE="run${RUNNO}_${XPOS}_${ZPOS}_e${ENERGY}.ps"
    echo $OUTPUT_FILE

    X_MIN=`echo $line | cut -f7 -d ' '`
    X_MAX=`echo $line | cut -f8 -d ' '`
    E_MIN=`echo $line | cut -f9 -d ' '`
    E_MAX=`echo $line | cut -f10 -d ' '`
    DATA_TYPE=0

    
    # elif [ $ENERGY -eq "3" ]; then
    #    X_MIN=4650
    #    X_MAX=4750
    #    E_MIN=4.0
    #    E_MAX=4.07
    # fi

    echo root -l -b -q "fit_linac_data.C(${DATA_DIR}${DATA_FILE},${DATA_TYPE},${X_MIN},${X_MAX},${E_MIN},${E_MAX},${OUTPUT_DIR}${OUTPUT_FILE},${ENERGY},${XPOS},${ZPOS})"
    root -l -b -q "fit_linac_data.C(\"${DATA_DIR}${DATA_FILE}\",${DATA_TYPE},${X_MIN},${X_MAX},${E_MIN},${E_MAX},\"${OUTPUT_DIR}${OUTPUT_FILE}\",${ENERGY},\"${XPOS}\",\"${ZPOS}\")"

    cd $OUTPUT_DIR
    ps2pdf ${OUTPUT_FILE}
    cd -

    SUMMARY_TABLE="ge_fit_summary.txt"
    P_TMP=`cat "bestfit_momentum_tmp.txt"`
    echo $RUNNO $P_TMP >> $SUMMARY_TABLE
    
done
   
