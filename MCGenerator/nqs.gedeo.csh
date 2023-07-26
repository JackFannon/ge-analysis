#!/bin/csh -f
make gedeo

#set i = 5000
#while ($i <= 5100)

#set i = 5900
#while ($i <= 6100)

#set i = 6800
#while ($i <= 7000)

#set i = 6900
#while ($i <= 7010)

#set i = 8800
#while ($i <= 8900)

#set i = 8900
#while ($i <= 9000)

set i = 13550
while($i <= 13700)

#set i = 13700
#while ($i <= 13800)

#set i = 16300
#while ($i <= 16500)

#set i = 18850
#while ($i <= 19050)

#set i = 8873
#while ($i <= 8873)

set ene = `printf "%05d" $i`
set f_err = /home/jfannon/software/GeAnalysis/MCGenerator/err/$ene.err
set f_out = /home/jfannon/software/GeAnalysis/MCGenerator/log/$ene.log
pjsub -L rscgrp=lowe -e $f_err -o $f_out script/gedeo$ene.csh
#qsub -q ALL -e $f_err -o $f_out script/gedeo$ene.csh

  @   i++
end
