#!/usr/bin/env python
# This script processes the input from the web form requesting validation
# against one or more CIF dictionaries
import CifFile
import validate_cif
import tempfile

import cgi
import os
import cgitb; cgitb.enable()      # for debugging

# output header
print "Content-Type: text/html\n\n"
            
formdata = cgi.FieldStorage()

# some constants
dic_directory = "/usr/local/lib/cif/"
input_cif = formdata["cif_file"]
if input_cif.file:
    filename = input_cif.filename
else: filename = ""
is_dic = formdata.has_key("is_dic")
input_dics = formdata.getlist("dic_list")

# Save file data to temporary file
tmpfile = tempfile.mkstemp()[1]
jj = open(tmpfile,"w")
jj.write(input_cif.file.read())
jj.close()
try:
    cf = CifFile.ReadCif(tmpfile,scantype="flex")
    os.remove(tmpfile)
except CifFile.StarFile.StarError:
    os.remove(tmpfile)
    import sys,re
    print "<h3>File reading error</h3>"
    print "<p>File %s appears to have one or more syntax errors</p>" % input_cif.filename
    print "<p>The detailed error message is as follows:</p><p><pre>"
    etype,eval,trace = sys.exc_info()
    unsafe_str = str(eval)
    unsafe_str = unsafe_str.replace("&","&amp;")
    unsafe_str = unsafe_str.replace("<","&lt;")
    safe_str = unsafe_str.replace(">","&gt;")
    print safe_str
    print "</pre>"
except:
    os.remove(tmpfile)
    print "Unspecified error reading file %s.  This is most likely a CIF syntax error."
# Now do the validation...
else:
    diclist = map(lambda a:os.path.join(dic_directory,a),input_dics)
    merged_dics = CifFile.merge_dic(diclist)
    validate_cif.output_header(True,filename,input_dics)
    print CifFile.validate_report(CifFile.validate(cf,dic= merged_dics,isdic=is_dic),use_html=True)
    validate_cif.output_footer(True)
