#!/bin/zsh

export LD_LIBRARY_PATH=`pwd`/../../../src/.libs
export RUBYLIB=`pwd`/../..
# set isotropy information
ISOTROPY=~/code/isotropy/findsym
export ISODATA=~/code/isotropy/


for dir in `/bin/ls -d poscar*`; do
	echo $dir
	cd $dir
	for i in `ls POSCAR-*`;do
		echo -n "  "$i": "
		spglib=`ruby $RUBYLIB/symPoscar.rb -n $i|awk -F"(" '{print $2}'|sed s/\)//`
		isotropy=`ruby $RUBYLIB/poscar2findsym.rb $i|$ISOTROPY|grep Space|awk '{print $3}'`

		if [ $spglib = $isotropy ]; then
			echo 'ok'
		else
			echo 'fail (' $isotropy ')'
		fi
	done
	echo
	cd ..
done

