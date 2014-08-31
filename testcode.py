# coding: utf-8
with open("CHANGES") as f:
    contents = f.read()
    
contents
contents.split("-+")
toks = contents.split("-+")
print toks
print toks[0]
re.split("-+", contents)
import re
re.split("-+", contents)
toks = re.split("-+", contents)
print toks[0]
print toks[1]
get_ipython().magic(u'save 1-13 testcode.py')
