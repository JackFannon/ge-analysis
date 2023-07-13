#!/bin/bash
H2ROOT=/usr/local/sklib_gcc8/skofl-trunk/bin/h2root
for FILE in `ls hbk/*.hbk`; do
#for FILE in `ls hbk/163[5-9]*.hbk`; do
#for FILE in `ls hbk/164*.hbk`; do
#for FILE in `ls hbk/07000.hbk`; do
    ${H2ROOT} $FILE
done

for FILE2 in `ls hbk/*.root`; do
    mv $FILE2 root/
done
