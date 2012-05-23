# Makefile for ANBF Python based Cif handling modules

#
package: CifFile.py StarFile.py Parsers documentation
	python setup.py sdist
	python setup.py bdist
#	python setup.py bdist_wininst
#
../PyCifRW.tar: clean package
	(cd ..; tar cvf PyCifRW.tar --exclude tests --exclude CVS --exclude yapps2 --exclude error_reports --exclude old_stuff PyCifRW)
#
%.py : %.nw
	notangle $< > $@
#
documentation: CifFile.nw YappsStarParser.nw StarFile.nw
	noweave -html -index -filter l2h CifFile.nw > CifFile.html
	noweave -html -index -filter l2h StarFile.nw > StarFile.html
	noweave -html -index -filter l2h YappsStarParser.nw > YappsStarParser.html
# 
Parsers: YappsStarParser_DDLm.py YappsStarParser_1_1.py YappsStarParser_1_0.py
	
#
clean: 
	rm -f *.pyc *.g
#
YappsStarParser_1_0.py: YappsStarParser.nw
	notangle -R1.0_syntax YappsStarParser.nw > YappsStarParser_1_0.g
	python ./yapps3/yapps2.py YappsStarParser_1_0.g
#
YappsStarParser_1_1.py: YappsStarParser.nw
	notangle -R1.1_syntax YappsStarParser.nw > YappsStarParser_1_1.g
	python ./yapps3/yapps2.py YappsStarParser_1_1.g
#
YappsStarParser_DDLm.py: YappsStarParser.nw
	notangle -RDDLm_syntax YappsStarParser.nw > YappsStarParser_DDLm.g
	python ./yapps3/yapps2.py YappsStarParser_DDLm.g
#

