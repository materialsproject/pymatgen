#
# An example of how to output a subset of looped items.
#
from PyCifRW import CifFile
cf = CifFile.CifFile("martin.cif")["some_stuff"] # open and parse our cif, 
                                      #we want data block named "some_stuff".

#######################
# This section is optional: we check that all our data
# items exist before attempting to access them
needed_items = [
          "_atom_site_label",
          "_atom_site_type_symbol",
          "_atom_site_fract_x",
          "_atom_site_fract_y",
          "_atom_site_fract_z"]

loopitems = cf.GetLoop("_atom_site_label")  #get item names and values
loopkeys = map(lambda a:a[0],loopitems)     #get item names only
if len(filter(lambda a,b=loopkeys:a not in b, needed_items)) != 0:
    print "Error: one or more items missing from atom_site_label loop"
    exit                        
#
# End of optional section
##########################
nb = CifFile.CifBlock()                   # create a new block 
# Now we add data from the original file 
nb.AddCifItem((
         needed_items,
         map(lambda a:cf[a],needed_items)
         ))
df = CifFile.CifFile()                    # create a new cif object
df.NewBlock("changed",nb)             # and add our new block
outfile = open("updated.cif",'w')      #open a file to write to
outfile.write (df.WriteOut(comment="# This file has been updated"))

# Note that the above AddCifItem call expands to be equivalent to
# the call below; we saved our typing energy for more important
# things.
#
# nb.AddCifItem 
#         ((
#         [
#         "_atom_site_label",
#         "_atom_site_type_symbol",
#         "_atom_site_fract_x",
#         "_atom_site_fract_y",
#         "_atom_site_fract_z"]
#        [cf["_atom_site_label"],
#         cf["_atom_site_type_symbol"],
#         cf["_atom_site_fract_x"],
#         cf["_atom_site_fract_y"],
#         cf["_atom_site_fract_z"]]
#         ))
