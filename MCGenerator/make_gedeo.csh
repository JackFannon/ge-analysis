#!/bin/csh -f

set top_dir = `pwd`

foreach ene ( `ls ./ene` )
echo $ene

cat <<! >! script/gedeo$ene.csh
#! /bin/csh -f
hostname
date
cd $top_dir
time ./gedeo card/gedeo$ene.card hbk/$ene.hbk

!

chmod u+x script/gedeo$ene.csh

end
