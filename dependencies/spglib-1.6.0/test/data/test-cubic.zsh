#!/bin/zsh

export LD_LIBRARY_PATH=`pwd`/../../src/.libs
export RUBYLIB=`pwd`/..

findsym=$1
# this is special check using findsym
if [ x$findsym = "xfindsym" ]; then
	export ISODATA=~/tools/isotropy/
fi

for i in `/bin/ls cubic/POSCAR-*`;do
	spg=`ruby ../symPoscar.rb -n $i`
	numspg=`echo $spg|awk -F"(" '{print $2}'|sed s/\)//`
	numposcar=`echo $i|awk -F"/" '{print $2}'|cut -c 8-10|awk '{print $1*1}'`

	if [ x$findsym = "xfindsym" ]; then
		numfindsym=`../poscar2findsym.rb $i | $ISODATA/findsym|grep Space|awk '{print $3}'`
#		spgfindsym=`../poscar2findsym.rb $i | $ISODATA/findsym|grep Space|awk '{print $5}'`
	fi

	if [ $numspg = $numposcar ]; then
		echo -n 'ok'
	else
		echo -n 'fail (' $numposcar ')'
	fi

	if [ x$findsym = "xfindsym" ]; then
		if [ $numspg = $numfindsym ]; then
			echo " ..." $spg '(ok)' 
		else
			echo " ..." $spg ' fail (' $numfindsym ')'
		fi
	else
		echo " ..." $spg 
	fi
done

