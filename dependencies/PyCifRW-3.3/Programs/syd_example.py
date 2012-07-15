# Example python program: as close to ciftbx as we can make it

from CifFile import *     #import definitions
import sys                #to finish early

# open our dictionary
try:
    my_dict = CifDic("cif_core.dic")
except IOError:
    print "Cannot open cif_core.dic"
    sys.exit()

# open our CIF
try:
    my_cif = ValidCifFile(datasource="test.cif",dic=my_dict)#read our CIF file
except IOError:
    print "Cannot open test.cif"
    sys.exit()
except ValidCifError,error_message:
    print "test.cif failed validity checks:"
    print error_message
    sys.exit()
except CifError, error_message:
    print "Syntax error in test.cif:"
    print error_message
    sys.exit()

# get the first blockname
first_block = my_cif.keys()[0]

# store in variable for convenience
my_data_block = my_cif[first_block]

# get some data
cela,siga = get_number_with_esd(my_data_block["_cell_length_a"])#cell dimension
name = my_data_block["_symmetry_cell_setting"]      #cell setting

# get a random data name which is not one of the above
allnames = my_data_block.keys()
allnames.remove("_cell_length_a")
allnames.remove("_symmetry_cell_setting")
data = my_data_block[allnames[0]]

# to print, don't need to check type
print "%s  %s" % (allnames[0],data)

# loop atom sites
names = my_data_block["_atom_site_label"]
xsxs = my_data_block["_atom_site_fract_x"]
as_numbers = map(get_number_with_esd,xsxs)
processed = map(None,names,as_numbers) 
for label, data in processed: 
    print "%s  %d   %d" % (label,data[0],data[1])
