# Example python program: native use

from CifFile import *     #import definitions

# only proceed if CifFile is valid
my_cif = ValidCifFile(datasource="test.cif",dic="cif_core.dic")

# get the first blockname
first_block = my_cif.keys()[0]

# store in variable for convenience
my_data_block = my_cif[first_block]

# get some data
cela = my_data_block["_cell_length_a"] #cell dimension
name = my_data_block["_symmetry_cell_setting"] #cell setting

# get a random data name which may be one of the above
data = my_data_block[allnames[0]]

# to print, don't need to check type
print "%s  %s" % (allnames[0],data)

# loop atom sites
names = my_data_block["_atom_site_label"]
xsxs = my_data_block["_atom_site_fract_x"]
processed = map(None,names,xsxs) 
for label, data in processed: 
    print "%s  %s" % (label,data)
