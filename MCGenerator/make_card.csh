#!/bin/csh -f

foreach nsub ( `ls ./ene` )
echo $nsub
cat card/gedeo.card1 ene/$nsub card/gedeo.card2 >! card/gedeo$nsub.card
end
