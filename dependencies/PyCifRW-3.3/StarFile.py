"""
1.This Software copyright \u00A9 Australian Synchrotron Research Program Inc, ("ASRP").

2.Subject to ensuring that this copyright notice and licence terms
appear on all copies and all modified versions, of PyCIFRW computer
code ("this Software"), a royalty-free non-exclusive licence is hereby
given (i) to use, copy and modify this Software including the use of
reasonable portions of it in other software and (ii) to publish,
bundle and otherwise re-distribute this Software or modified versions
of this Software to third parties, provided that this copyright notice
and terms are clearly shown as applying to all parts of software
derived from this Software on each occasion it is published, bundled
or re-distributed.  You are encouraged to communicate useful
modifications to ASRP for inclusion for future versions.

3.No part of this Software may be sold as a standalone package.

4.If any part of this Software is bundled with Software that is sold,
a free copy of the relevant version of this Software must be made
available through the same distribution channel (be that web server,
tape, CD or otherwise).

5.It is a term of exercise of any of the above royalty free licence
rights that ASRP gives no warranty, undertaking or representation
whatsoever whether express or implied by statute, common law, custom
or otherwise, in respect of this Software or any part of it.  Without
limiting the generality of the preceding sentence, ASRP will not be
liable for any injury, loss or damage (including consequential loss or
damage) or other loss, loss of profits, costs, charges or expenses
however caused which may be suffered, incurred or arise directly or
indirectly in respect of this Software.

6. This Software is not licenced for use in medical applications.
"""

from types import *
from urllib import *         # for arbitrary opening
import re
import copy
class StarList(list):
    pass

# Because DDLm makes a tuple from a tuple...
class StarTuple(tuple):
    def __new__(cls,*arglist):
        return tuple.__new__(cls,arglist)

class StarDict(dict):
    pass

class LoopBlock:
    def __init__(self,data = (), dimension = 0, maxoutlength=2048, wraplength=80, overwrite=True):
        # print 'Creating new loop block, dimension %d' % dimension
        self.block = {}
        self.loops = []
        self.no_packets = 0
        self.item_order = []
        self.lower_keys = []    #for efficiency
        self.comment_list = {}
        self.dimension = dimension
        self.popout = False         #used during load iteration
        self.curitem = -1           #used during iteration
        self.maxoutlength = maxoutlength
        self.wraplength = wraplength
        self.overwrite = overwrite
        if not hasattr(self,'loopclass'):  #in case are derived class
            self.loopclass = LoopBlock  #when making new loops
        self.char_check = re.compile("[][ \n\r\t!%&\(\)*+,./:<=>?@0-9A-Za-z\\\\^`{}\|~\"#$';_-]+",re.M)
        if isinstance(data,(TupleType,ListType)):
            for item in data:
                self.AddLoopItem(item)
        elif isinstance(data,LoopBlock):
            self.block = data.block.copy() 
            self.item_order = data.item_order[:]
            self.lower_keys = data.lower_keys[:]
            self.comment_list = data.comment_list.copy()
            self.dimension = data.dimension
            # loops as well; change loop class 
            for loopno in range(len(data.loops)):
                try:
                    placeholder = self.item_order.index(data.loops[loopno])
                except ValueError:
                    print "Warning: loop %s (%s) in loops, but not in item_order (%s)" % (`data.loops[loopno]`,str(data.loops[loopno]),`self.item_order`)
                    placeholder = -1
                self.item_order.remove(data.loops[loopno])   #gone
                newobject = self.loopclass(data.loops[loopno])
                # print "Recasting and adding loop %s -> %s" % (`data.loops[loopno]`,`newobject`)
                self.insert_loop(newobject,position=placeholder)

    def __str__(self):
        return self.printsection()

    def __setitem__(self,key,value):
        # catch a one member loop, for convenience
        # we assume the key is a string value only
        self.AddLoopItem((key,value))

    def __getitem__(self,key):
        if isinstance(key,IntType):   #return a packet!!
            return self.GetPacket(key)        
        return self.GetLoopItem(key)

    def __delitem__(self,key):
        self.RemoveLoopItem(key)

    def __len__(self):
        blen = len(self.block)
        for aloop in self.loops:
            # print 'Aloop is %s' % `aloop`
            blen = blen + len(aloop)  # also a LoopBlock
        return blen    

    def __nonzero__(self):
        if self.__len__() > 0: return 1
        return 0

    # keys returns all internal keys
    def keys(self):
        thesekeys = self.block.keys()
        for aloop in self.loops:
            thesekeys.extend(aloop.keys())
        return thesekeys

    def values(self):
        ourkeys = self.keys()
        return map(lambda a:self[a],ourkeys)

    def items(self):
        ourkeys = self.keys()
        return map(lambda a,b:(a,b),self.keys(),self.values())

    def has_key(self,key):
        if key.lower() in self.lower_keys:
            return 1
        for aloop in self.loops:
            if aloop.has_key(key): return 1
        return 0

    def get(self,key,default=None):
        if self.has_key(key):
            retval = self.GetLoopItem(key)
        else:
            retval = default
        return retval

    def clear(self):
        self.block = {}
        self.loops = []
        self.item_order = []
        self.lower_keys = []
        self.no_packets = 0

    # doesn't appear to work
    def copy(self):
        newcopy = self.copy.im_class(dimension = self.dimension)
        newcopy.block = self.block.copy()
        newcopy.loops = []
        newcopy.no_packets = self.no_packets
        newcopy.item_order = self.item_order[:]
        newcopy.lower_keys = self.lower_keys[:]
        for loop in self.loops:
            try:
                placeholder = self.item_order.index(loop)
            except ValueError:
                print "Warning: loop %s (%s) in loops, but not in item_order (%s)" % (`loop`,str(loop),`self.item_order`)
                placeholder = -1
            newcopy.item_order.remove(loop)   #gone
            newobject = loop.copy()
            # print "Adding loop %s -> %s" % (`loop`,`newobject`)
            newcopy.insert_loop(newobject,position=placeholder)
        return newcopy

    # this is not appropriate for subloops.  Instead, the loop block
    # should be accessed directly for update
     
    def update(self,adict):
        for key in adict.keys():
            self.AddLoopItem((key,adict[key]))

    def load_iter(self,coords=[]):
        count = 0        #to create packet index 
        while not self.popout:
            # ok, we have a new packet:  append a list to our subloops
            for aloop in self.loops:
                aloop.new_enclosing_packet()
            for iname in self.item_order:
                if isinstance(iname,LoopBlock):       #into a nested loop
                    for subitems in iname.load_iter(coords=coords+[count]):
                        # print 'Yielding %s' % `subitems`
                        yield subitems
                    # print 'End of internal loop'
                else:
                    if self.dimension == 0:
                        # print 'Yielding %s' % `self[iname]`
                        yield self,self[iname]
                    else:
                        backval = self.block[iname]
                        for i in range(len(coords)):
                           # print 'backval, coords: %s, %s' % (`backval`,`coords`)
                           backval = backval[coords[i]]
                        yield self,backval
            count = count + 1      # count packets
        self.popout = False        # reinitialise
        # print 'Finished iterating'
        yield self,'###Blank###'     #this value should never be used

    # an experimental fast iterator for level-1 loops (ie CIF)
    def fast_load_iter(self):
        targets = map(lambda a:self.block[a],self.item_order)
        while targets:
            for target in targets:
                yield self,target

    # Add another list of the required shape to take into account a new outer packet
    def new_enclosing_packet(self):
        if self.dimension > 1:      #otherwise have a top-level list
            for iname in self.keys():  #includes lower levels
                target_list = self[iname]
                for i in range(3,self.dimension): #dim 2 upwards are lists of lists of... 
                    target_list = target_list[-1]
                target_list.append([])
                # print '%s now %s' % (iname,`self[iname]`)

    def recursive_iter(self,dict_so_far={},coord=[]):
        # print "Recursive iter: coord %s, keys %s, dim %d" % (`coord`,`self.block.keys()`,self.dimension)
        my_length = 0
        top_items = self.block.items()
        top_values = self.block.values()       #same order as items
        drill_values = self.block.values()
        for dimup in range(0,self.dimension):  #look higher in the tree
            if len(drill_values)>0:            #this block has values
                drill_values=drill_values[0]   #drill in
            else:
                raise StarError("Malformed loop packet %s" % `top_items[0]`)
        my_length = len(drill_values)
        if self.dimension == 0:                #top level
            for aloop in self.loops:
                for apacket in aloop.recursive_iter():
                    # print "Recursive yielding %s" % `dict(top_items + apacket.items())`
                    prep_yield = StarPacket(top_values+apacket.values())  #straight list
                    for name,value in top_items + apacket.items():
                        setattr(prep_yield,name,value)
                    yield prep_yield
        else:                                  #in some loop
            for i in range(my_length):
                kvpairs = map(lambda a:(a,self.coord_to_group(a,coord)[i]),self.block.keys())
                kvvals = map(lambda a:a[1],kvpairs)   #just values
                # print "Recursive kvpairs at %d: %s" % (i,`kvpairs`)
                if self.loops:
                  for aloop in self.loops:
                    for apacket in aloop.recursive_iter(coord=coord+[i]):
                        # print "Recursive yielding %s" % `dict(kvpairs + apacket.items())`
                        prep_yield = StarPacket(kvvals+apacket.values())
                        for name,value in kvpairs + apacket.items():
                            setattr(prep_yield,name,value)
                        yield prep_yield
                else:           # we're at the bottom of the tree
                    # print "Recursive yielding %s" % `dict(kvpairs)`
                    prep_yield = StarPacket(kvvals)
                    for name,value in kvpairs:
                        setattr(prep_yield,name,value)
                    yield prep_yield

    # small function to use the coordinates. 
    def coord_to_group(self,dataname,coords):
          if not isinstance(dataname,StringType):
             return dataname     # flag inner loop processing
          newm = self[dataname]          # newm must be a list or tuple
          for c in coords:
              # print "Coord_to_group: %s ->" % (`newm`),
              newm = newm[c]
              # print `newm`
          return newm 

    def flat_iterator(self):
        if self.dimension == 0:   
            yield copy.copy(self)
        else:
            my_length = 0
            top_keys = self.block.keys()
            if len(top_keys)>0:
                my_length = len(self.block[top_keys[0]])
            for pack_no in range(my_length):
                yield(self.collapse(pack_no))
            

    def insert_loop(self,newloop,position=-1,audit=True):
        # check that new loop is kosher
        if newloop.dimension != self.dimension + 1:
            raise StarError( 'Insertion of loop of wrong nesting level %d, should be %d' % (newloop.dimension, self.dimension+1))
        self.loops.append(newloop)
        if audit:
            dupes = self.audit()
            if dupes:
                dupenames = map(lambda a:a[0],dupes)
                raise StarError( 'Duplicate names: %s' % `dupenames`)
        if position >= 0:
            self.item_order.insert(position,newloop)
        else:
            self.item_order.append(newloop)
        # print "Insert loop: item_order now" + `self.item_order`

    def remove_loop(self,oldloop):
        # print "Removing %s: item_order %s" % (`oldloop`,self.item_order)
        # print "Length %d" % len(oldloop)
        self.item_order.remove(oldloop)
        self.loops.remove(oldloop)
     
    def AddComment(self,itemname,comment):
        self.comment_list[itemname.lower()] = comment

    def RemoveComment(self,itemname):
        del self.comment_list[itemname.lower()]

    def GetLoopItem(self,itemname):
        # assume case is correct first
        try:
            return self.block[itemname]
        except KeyError:
            for loop in self.loops:
                try:
                    return loop[itemname]
                except KeyError:
                    pass
        if itemname.lower() not in self.lower_keys:
            raise KeyError, 'Item %s not in block' % itemname
        # it is there somewhere, now we need to find it
        real_keys = self.block.keys()
        lower_keys = map(lambda a:a.lower(),self.block.keys()) 
        try:
            k_index = lower_keys.index(itemname.lower())
        except ValueError:
            raise KeyError, 'Item %s not in block' % itemname
        return self.block[real_keys[k_index]]

    def RemoveLoopItem(self,itemname):
        if self.has_key(itemname):
            testkey = itemname.lower()
            real_keys = self.block.keys()
            lower_keys = map(lambda a:a.lower(),real_keys)
            try:
                k_index = lower_keys.index(testkey)
            except ValueError:    #must be in a lower loop
                for aloop in self.loops:
                    if aloop.has_key(itemname):
                        # print "Deleting %s (%s)" % (itemname,aloop[itemname])
                        del aloop[itemname]
                        if len(aloop)==0:  # all gone
                           self.remove_loop(aloop)
                        break
            else:
              del self.block[real_keys[k_index]]
              self.lower_keys.remove(testkey)
              # now remove the key in the order list
              for i in range(len(self.item_order)):
                if isinstance(self.item_order[i],StringType): #may be loop
                    if self.item_order[i].lower()==testkey:
                        del self.item_order[i]
                        break
            if len(self.block)==0:    #no items in loop, length -> 0
                self.no_packets = 0
            return        #no duplicates, no more checking needed

    def AddLoopItem(self,data,precheck=False,maxlength=-1):
        # print "Received data %s" % `data`
        # we accept only tuples, strings and lists!!
        if isinstance(data[0],(TupleType,ListType)):
           # internal loop
           # first we remove any occurences of these datanames in
           # other loops
           for one_item in data[0]:
               if self.has_key(one_item):
                   if not self.overwrite:
                       raise StarError( 'Attempt to insert duplicate item name %s' % data[0])
                   else:
                       del self[one_item]
           newloop = self.loopclass(dimension = self.dimension+1)
           keyvals = zip(data[0],data[1])
           for key,val in keyvals:
               newloop.AddLoopItem((key,val))
           self.insert_loop(newloop)
        elif not isinstance(data[0],StringType):
                  raise TypeError, 'Star datanames are strings only (got %s)' % `data[0]`
        else:
           if data[1] == [] or get_dim(data[1])[0] == self.dimension:
               if not precheck:
                   self.check_data_name(data[0],maxlength)    # make sure no nasty characters   
               # check that we can replace data
               if not self.overwrite:
                   if self.has_key(data[0]):
                       raise StarError( 'Attempt to insert duplicate item name %s' % data[0])
               # now make sure the data is OK type
               regval = self.regularise_data(data[1])
               if not precheck:
                   try:
                       self.check_item_value(regval)
                   except StarError, errmes:
                       raise StarError( "Item name " + data[0] + " " + `errmes`)
               if self.dimension > 0:
                   if self.no_packets <= 0:
                       self.no_packets = len(data[1])  #first item in this loop
                   if len(data[1]) != self.no_packets:
                       raise StarLengthError, 'Not enough values supplied for %s' % (data[0])
               try:
                   oldpos = self.GetItemPosition(data[0])
               except ValueError:
                   oldpos = len(self.item_order)#end of list 
               self.RemoveLoopItem(data[0])     # may be different case, so have to do this
               self.block.update({data[0]:regval})  # trust the data is OK
               self.lower_keys.insert(oldpos,data[0].lower())
               self.item_order.insert(oldpos,data[0])
               #    self.lower_keys.append(data[0].lower())
               #    self.item_order.append(data[0])
                
           else:            #dimension mismatch
               raise StarLengthError, "input data dim %d != required dim %d: %s %s" % (get_dim(data[1])[0],self.dimension,data[0],`data[1]`)

    def check_data_name(self,dataname,maxlength=-1): 
        if maxlength > 0:
            if len(dataname)>maxlength:
                raise StarError( 'Dataname %s exceeds maximum length %d' % (dataname,maxlength))
        if dataname[0]!='_':
            raise StarError( 'Dataname ' + dataname + ' does not begin with _')
        if len (filter (lambda a: ord(a) < 33 or ord(a) > 126, dataname)) > 0:
            raise StarError( 'Dataname ' + dataname + ' contains forbidden characters')

    def check_item_value(self,item):
        test_item = item
        if type(item) != TupleType and type(item) != ListType:
           test_item = [item]         #single item list
        def check_one (it):
            if type(it) == StringType:
                if it=='': return
                me = self.char_check.match(it)            
                if not me:
                    raise StarError( 'Bad character in %s' % it)
                else:
                    if me.span() != (0,len(it)):
                        raise StarError('Data item "' + it + '"... contains forbidden characters')
        map(check_one,test_item)

    def regularise_data(self,dataitem):
        alrighttypes = [IntType, LongType, 
                        FloatType, StringType]
        okmappingtypes = [TupleType, ListType]
        thistype = type(dataitem)
        if thistype in alrighttypes or thistype in okmappingtypes:
            return dataitem
        if isinstance(dataitem,StarTuple) or \
           isinstance(dataitem,StarList) or \
           isinstance(dataitem,StarDict):
            return dataitem
        # so try to make into a list
        try:
            regval = list(dataitem)
        except TypeError, value:
            raise StarError( str(dataitem) + ' is wrong type for data value\n' )
        return regval
        
    def GetLoop(self,keyname):
        if keyname in self.block:        #python 2.2 or above
            return self
        for aloop in self.loops:
            try: 
                return aloop.GetLoop(keyname)
            except KeyError:
                pass
        raise KeyError, 'Item %s does not exist' % keyname

    def GetPacket(self,index):
        thispack = StarPacket([])
        for myitem in self.item_order:
            if isinstance(myitem,LoopBlock):
                pack_list = map(lambda b:myitem[b][index],myitem.item_order)
                # print 'Pack_list -> %s' % `pack_list`
                thispack.append(pack_list)
            elif self.dimension==0:
                thispack.append(self[myitem])
            else:
                thispack.append(self[myitem][index])
                setattr(thispack,myitem,thispack[-1])
        return thispack 

    def AddPacket(self,packet):
        if self.dimension==0:
            raise StarError,"Attempt to add packet to top level block"
        for myitem in self.item_order:
            self[myitem] = list(self[myitem])   #in case we have stored a tuple
            self[myitem].append(packet.__getattribute__(myitem))
        self.no_packets +=1 
            # print "%s now %s" % (myitem,`self[myitem]`)
        
    def RemoveKeyedPacket(self,keyname,keyvalue):
        packet_coord = list(self[keyname]).index(keyvalue)
        loophandle = self.GetLoop(keyname)
        for packet_entry in loophandle.item_order:
            loophandle[packet_entry] = list(loophandle[packet_entry])
            del loophandle[packet_entry][packet_coord]
        self.no_packets -= 1
        
    def GetKeyedPacket(self,keyname,keyvalue):
        #print "Looking for %s in %s" % (keyvalue, self[keyname])
        one_pack= filter(lambda a:getattr(a,keyname)==keyvalue,self)
        if len(one_pack)!=1:
            raise KeyError, "Bad packet key %s = %s: returned %d packets" % (keyname,keyvalue,len(one_pack))
        #print "Keyed packet: %s" % one_pack[0]
        return one_pack[0]

    def GetItemOrder(self):
        return self.item_order[:]

    def ChangeItemOrder(self,itemname,newpos):
        testpos = self.GetItemPosition(itemname)
        del self.item_order[testpos]
        # so we have an object ready for action
        self.item_order.insert(newpos,itemname)

    def GetItemPosition(self,itemname):
        import string
        def low_case(item):
            try:
                return string.lower(item)
            except AttributeError:
                return item
        try:
            testname = string.lower(itemname)
        except AttributeError: 
            testname = itemname
        lowcase_order = map(low_case,self.item_order)
        return lowcase_order.index(testname)

    def collapse(self,packet_no):
        if self.dimension == 0:
            raise StarError( "Attempt to select non-existent packet")
        newlb = LoopBlock(dimension=self.dimension-1)
        for one_item in self.item_order:
            if isinstance(one_item,LoopBlock):
                newlb.insert_loop(one_item.collapse(packet_no))
            else:
                # print "Collapse: %s -> %s" % (one_item,`self[one_item][packet_no]`)
                newlb[one_item] = self[one_item][packet_no] 
        return newlb
        
    def audit(self):
        import sets
        allkeys = self.keys()
        uniquenames = sets.Set(allkeys)
        if len(uniquenames) == len(allkeys): return []
        else:              
            keycount = map(lambda a:(a,allkeys.count(a)),uniquenames)
            return filter(lambda a:a[1]>1,keycount)
        
    def GetLoopNames(self,keyname):
        if keyname in self:
            return self.keys()
        for aloop in self.loops:
            try: 
                return aloop.GetLoopNames(keyname)
            except KeyError:
                pass
        raise KeyError, 'Item does not exist'

    def AddToLoop(self,dataname,loopdata):
        thisloop = self.GetLoop(dataname)
        for itemname,itemvalue in loopdata.items():
            thisloop[itemname] = itemvalue

    def SetOutputLength(self,wraplength=80,maxoutlength=2048):
        if wraplength > maxoutlength:
            raise StarError("Wrap length (requested %d) must be <= Maximum line length (requested %d)" % (wraplength,maxoutlength))
        self.wraplength = wraplength
        self.maxoutlength = maxoutlength
        for loop in self.loops:
            loop.SetOutputLength(wraplength,maxoutlength)

    def printsection(self,instring='',blockstart="",blockend="",indent=0,coord=[]):
        import cStringIO
        import string
        # first make an ordering
        order = self.item_order[:]
        # now do it...
        if not instring:
            outstring = cStringIO.StringIO()       # the returned string
        else:
            outstring = instring
        if not coord:
            coords = [0]*(self.dimension-1)
        else:
            coords = coord
        if(len(coords)<self.dimension-1):
            raise StarError("Not enough block packet coordinates to uniquely define data")
        # print loop delimiter
        outstring.write(blockstart)
        while len(order)>0:
            # print "Order now: " + `order`
            itemname = order.pop(0)
            if self.dimension == 0:            # ie value next to tag
                if not isinstance(itemname,LoopBlock):  #no loop
                   # grab any comment
                   thiscomment = self.comment_list.get(itemname.lower(),'') 
                   itemvalue = self[itemname]
                   if isinstance(itemvalue,StringType):  #need to sanitize
                         thisstring = self._formatstring(itemvalue)
                   else: thisstring = str(itemvalue)
                   # try for a tabstop at 40
                   if len(itemname)<40 and (len(thisstring)-40 < self.wraplength-1):
                       itemname = itemname + ' '*(40-len(itemname))
                   else: itemname = itemname + ' '
                   if len(thisstring) + len(itemname) < (self.wraplength-1):
                         outstring.write('%s%s' % (itemname,thisstring))
                         if thiscomment:
                             if len(thiscomment)+len(thisstring)+len(itemname)< (self.wraplength-3):
                                 outstring.write(' #'+thiscomment)
                   else:
                         outstring.write('%s\n %s' % (itemname, thisstring))
                         if thiscomment:
                             if len(thiscomment)+len(thisstring)<(self.wraplength-3):
                                 outstring.write(' #'+thiscomment)
                             else:
                                 outstring.write('\n#'+thiscomment)
                   outstring.write('\n')
                else:   # we are asked to print an internal loop block
                    #first make sure we have sensible coords.  Length should be one
                    #less than the current dimension
                    outstring.write(' '*indent); outstring.write('loop_\n')
                    itemname.format_names(outstring,indent+2)
                    itemname.format_packets(outstring,coords,indent+2)
            else:   # we are a nested loop
                outstring.write(' '*indent); outstring.write('loop_\n')
                self.format_names(outstring,indent+2)
                self.format_packets(outstring,coords,indent+2)
        if instring: return   #inside a recursion
        else:
            returnstring = outstring.getvalue()
        outstring.close()
        return returnstring

    def format_names(self,outstring,indent=0):
        temp_order = self.item_order[:]
        while len(temp_order)>0:
            itemname = temp_order.pop(0)
            if isinstance(itemname,StringType):  #(not loop)
                outstring.write(' ' * indent) 
                outstring.write(itemname)
                outstring.write("\n")
            else:                                # a loop
                outstring.write(' ' * indent) 
                outstring.write("loop_\n")
                itemname.format_names(outstring,indent+2)
                outstring.write(" stop_\n")

    def format_packets(self,outstring,coordinates,indent=0):
       import cStringIO
       import string
       # get our current group of data
       # print 'Coords: %s' % `coordinates`
       alldata = map(lambda a:self.coord_to_group(a,coordinates),self.item_order)
       # print 'Alldata: %s' % `alldata`
       packet_data = apply(zip,alldata)
       # print 'Packet data: %s' % `packet_data`
       curstring = ''
       for position in range(len(packet_data)):
           for point in range(len(packet_data[position])):
               datapoint = packet_data[position][point]
               packstring = self.format_packet_item(datapoint,indent)
               if len(curstring) + len(packstring)> self.wraplength-2: #past end of line with space
                   curstring = curstring + '\n' + ' '*indent + packstring
               elif curstring == '':
                   curstring = curstring + ' '*indent + packstring
               else:
                   curstring = curstring + ' ' + packstring
           outstring.write(curstring + '\n')     #end of one packet
           curstring = ''
       outstring.write(' ' + curstring + '\n')    #last time through
               
    def format_packet_item(self,pack_item,indent):
        # print 'Formatting %s' % `pack_item`
        curstring = ''
        if isinstance(pack_item,(StringType,IntType,FloatType,LongType,StarTuple,StarList)):
           if isinstance(pack_item,StringType):
               thisstring = self._formatstring(pack_item) #no spaces yet
               if '\n' in thisstring:    #must have semicolon digraph then 
                   curstring = curstring + thisstring
                   curstring = curstring + (' ' * indent)
                   thisstring = ''
           else: 
               thisstring = '%s' % str(pack_item)
           if len(curstring) + len(thisstring)> self.wraplength-2: #past end of line with space
               curstring = curstring + '\n' #add the space
               curstring = curstring + (' ' * indent) + thisstring
           else: 
               curstring = curstring + ' ' + thisstring
        # Now, for each nested loop we call ourselves again
        # After first outputting the current line
        else:               # a nested packet
           if not isinstance(pack_item[0],(ListType,TupleType)):  #base packet
               item_list = pack_item
           else:
               item_list = apply(zip,pack_item)
           for sub_item in item_list:
               curstring = curstring + ' ' + self.format_packet_item(sub_item,indent)
           # stop_ is not issued at the end of each innermost packet
           if isinstance(pack_item[0],(ListType,TupleType)):
               curstring = curstring + ' stop_ '
        return curstring          

    def _formatstring(self,instring):
        import string
        if len(instring)==0: return "''"
        if len(instring)< (self.maxoutlength-2) and '\n' not in instring and not ('"' in instring and '\'' in instring):
            if not ' ' in instring and not '\t' in instring and not '\v' \
              in instring and not '_' in instring and not (instring[0]=="'" or \
                 instring[0]=='"'):                  # no blanks
                return instring
            if not "'" in instring:                                       #use apostrophes
                return "'%s'" % (instring)
            elif not "\"" in instring:
                return '"%s"' % (instring)
        # is a long one or one that needs semicolons due to carriage returns
        outstring = "\n;"
        # if there are returns in the string, try to work with them
        while 1:
            retin = string.find(instring,'\n')+1
            if retin < self.maxoutlength and retin > 0:      # honour this break
                outstring = outstring + instring[:retin]
                instring = instring[retin:]
            elif len(instring)<self.maxoutlength:            # finished
                outstring = outstring + instring + '\n;\n'
                break
            else:                             # find a space
                for letter in range(self.maxoutlength-1,self.wraplength-1,-1): 
                    if instring[letter] in ' \t\f': break
                outstring = outstring + instring[:letter+1]
                outstring = outstring + '\n'
                instring = instring[letter+1:]            
        return outstring



class StarBlock(LoopBlock):
    def __init__(self,*pos_args,**keyword_args):
        LoopBlock.__init__(self,*pos_args,**keyword_args)
        self.saves = BlockCollection(element_class=LoopBlock,type_tag="save")

    def __getitem__(self,key):
        if key == "saves":
            return self.saves
        else:
            return LoopBlock.__getitem__(self,key)

    def __setitem__(self,key,value):
        if key == "saves":
            self.saves[key] = value
        else:
            LoopBlock.__setitem__(self,key,value)

    def clear(self):
        LoopBlock.clear(self)
        self.saves = BlockCollection(element_class=LoopBlock,type_tag="save_")

    def copy(self):
        newblock = LoopBlock.copy(self)
        newblock.saves = self.saves.copy()
        return self.copy.im_class(newblock)   #catch inheritance

    def has_key(self,key):
        if key == "saves": return 1
        else: return LoopBlock.has_key(self,key)
        
    def __str__(self):
        retstr = ''
        for sb in self.saves.keys(): 
            retstr = retstr + '\nsave_%s\n\n' % sb
            self.saves[sb].SetOutputLength(self.wraplength,self.maxoutlength)
            retstr = retstr + str(self.saves[sb])
            retstr = retstr + '\nsave_\n\n'
        return retstr + LoopBlock.__str__(self)


class StarPacket(list):
    pass

class BlockCollection:
    def __init__(self,datasource=None,element_class=StarBlock,type_tag=''):
        self.dictionary = {}
        self.type_tag = type_tag
        self.lower_keys = []              # for efficiency
        self.element_class = element_class
        if isinstance(datasource,(DictType,BlockCollection)):
            for key,value in datasource.items():
                if value.__class__ == element_class:
                    self[key]=value
                else:
                    self[key]= element_class(value)
        self.header_comment = ''
     
    def __str__(self):
        return self.WriteOut()

    def __setitem__(self,key,value):
        if isinstance(value,(self.element_class,DictType)):
            self.NewBlock(key,value,replace=True)
        else: raise TypeError
        self.lower_keys.append(key.lower())

    # due to attempt to get upper/lower case treated as identical
    # we have a bit of cruft here
    def __getitem__(self,key):
        try:
            return self.dictionary[key]
        except KeyError:
            if key.lower() not in self.lower_keys:
                raise KeyError, "No such item: %s" % key
        curr_keys = self.dictionary.keys()
        lower_ordered = map(lambda a:a.lower(),curr_keys)
        keyindex = lower_ordered.index(key.lower())
        return self.dictionary[curr_keys[keyindex]]

    # we have to get an ordered list of the current keys,
    # as we'll have to delete one of them anyway
    def __delitem__(self,key):
        try:
            del self.dictionary[key]
            self.lower_keys.remove(key.lower())
        except KeyError:
            if not self.has_key(key):
                raise KeyError
            curr_keys = self.dictionary.keys()
            lower_ordered = map(lambda a:a.lower(),curr_keys)
            keyindex = lower_ordered.index(key.lower())
            del self.dictionary[curr_keys[keyindex]]
        
    def __len__(self):
        return len(self.dictionary)

    def keys(self):
        return self.dictionary.keys()

    # changes to take case independence into account
    def has_key(self,key):
        if not isinstance(key,StringType): return 0
        if self.dictionary.has_key(key):
           return 1
        if key.lower() in self.lower_keys:
           return 1
        return 0

    def get(self,key,default=None):
        if self.dictionary.has_key(key):
            return self.dictionary[key]
        elif self.has_key(key):     # take account of case
            return self.__getitem__(key)
        else:
            return default

    def clear(self):
        self.dictionary.clear()
        self.lower_keys = []

    def copy(self):   
        newcopy = self.dictionary.copy()
        return BlockCollection('',newcopy)
     
    def update(self,adict):
        for key in adict.keys():
            self.dictionary[key] = adict[key]
        self.lower_keys.extend(map(lambda a:a.lower(),adict.keys()))

    def items(self):
        return self.dictionary.items()

    def first_block(self):
        if self.keys():
            return self[self.keys()[0]]

    def NewBlock(self,blockname,blockcontents=(),replace=False,fix=True):
        if not blockcontents:
            blockcontents = self.element_class()
        elif isinstance(blockcontents,DictType):
            blockcontents = self.element_class(blockcontents)
        if not isinstance(blockcontents,self.element_class):
            raise StarError( 'Block is not of required type %s, is %s' % self.element_class.__name__,blockcontents.__class__.__name__)
        if fix:
            newblockname = re.sub('[  \t]','_',blockname)
        else: newblockname = blockname
        new_lowerbn = newblockname.lower()
        if self.lower_keys.count(new_lowerbn):    #already in CIF
            if not replace:
                raise StarError( "Attempt to replace existing block" + blockname)
            # generate a list of lower-case keys in correct order
            current_keys = self.dictionary.keys()
            blocknames = map(lambda a:a.lower(),current_keys)
            location = blocknames.index(new_lowerbn)
            del self.dictionary[current_keys[location]]
            self.lower_keys.remove(new_lowerbn)
        self.dictionary.update({blockname:blockcontents})
        self.lower_keys.append(new_lowerbn)

    def merge(self,new_bc,mode="strict",single_block=[],
                   idblock="",match_att=[],match_function=None):
        if single_block:
            self.dictionary[single_block[0]].merge(new_bc[single_block[1]],mode,
                                                   match_att=match_att,
                                                   match_function=match_function)
            return None
        base_keys = self.keys()
        block_to_item = base_keys   #default
        new_keys = new_bc.keys()
        if match_att:
            #make a blockname -> item name map
            if match_function:
                block_to_item = map(lambda a:match_function(self[a]),self.keys())
            else:
                block_to_item = map(lambda a:self[a].get(match_att[0],None),self.keys())
            #print `block_to_item`
        for key in new_keys:
            if key == idblock: continue
            basekey = key        #default value
            attval = new_bc[key].get(match_att[0],0)
            for ii in range(len(block_to_item)):  #do this way to get looped names
                thisatt = block_to_item[ii]
                #print "Looking for %s in %s" % (attval,thisatt)
                if attval == thisatt or \
                   (isinstance(thisatt,ListType) and attval in thisatt):
                      basekey = base_keys.pop(ii)
                      block_to_item.remove(thisatt)
                      break
            if not self.dictionary.has_key(basekey) or mode=="replace":
                self.dictionary[basekey] = new_bc[key]
            else:
                if mode=="strict":
                    raise StarError( "In strict merge mode: block %s in old and block %s in new files" % (basekey,key))
                elif mode=="overlay":
                    # print "Merging block %s with %s" % (basekey,key)
                    self.dictionary[basekey].merge(new_bc[key],mode,match_att=match_att)
                else:  
                    raise StarError( "Merge called with unknown mode %s" % mode)

    def get_all(self,item_name):
        raw_values = map(lambda a:self[a].get(item_name),self.dictionary.keys())
        raw_values = filter(lambda a:a != None, raw_values)
        ret_vals = []
        for rv in raw_values:
            if isinstance(rv,ListType):
                for rvv in rv:
                    if rvv not in ret_vals: ret_vals.append(rvv)
            else:
                if rv not in ret_vals: ret_vals.append(rv)
        return ret_vals

    def WriteOut(self,comment='',wraplength=80,maxoutlength=2048):
        import cStringIO
        if not comment:
            comment = self.header_comment
        outstring = cStringIO.StringIO()
        outstring.write(comment)
        for datablock in self.dictionary.keys():
            outstring.write('\n' + self.type_tag +datablock+'\n')
            self.dictionary[datablock].SetOutputLength(wraplength,maxoutlength)
            outstring.write(str(self.dictionary[datablock]))
        returnstring =  outstring.getvalue()
        outstring.close()
        return returnstring


class StarFile(BlockCollection):
    def __init__(self,datasource=None,maxinlength=-1,maxoutlength=0,blocktype=StarBlock,**kwargs):
        BlockCollection.__init__(self,datasource=datasource,element_class=blocktype,type_tag='data_')
        if isinstance(datasource, StarFile):
            self.my_uri = datasource.my_uri
        self.maxinlength = maxinlength      #no restriction
        if maxoutlength == 0:
            self.maxoutlength = 2048 
        else:
            self.maxoutlength = maxoutlength
        if type(datasource) is StringType or hasattr(datasource,"read"):
            newself = ReadStar(datasource,self.maxinlength,**kwargs)
            # print "Reinjecting by calling %s.__init__ with kwargs %s" % (`self.__init__.im_class`,kwargs)
            self.__init__.im_class.__init__(self,datasource=newself,maxoutlength=maxoutlength,**kwargs)
        self.header_comment = \
"""#\\#STAR
##########################################################################
#               STAR Format file 
#               Produced by PySTARRW module
# 
#  This is a STAR file.  STAR is a superset of the CIF file type.  For
#  more information, please refer to International Tables for Crystallography,
#  Volume G, Chapter 2.1
#
##########################################################################
"""
    def set_uri(self,my_uri): self.my_uri = my_uri


class StarError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return '\nStar Format error: '+ self.value 

class StarLengthError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return '\nStar length error: ' + self.value
def ReadStar(filename,maxlength=2048,dest=StarFile(),scantype='standard',grammar='1.1'):
    import string
    if grammar=="1.1":
        import YappsStarParser_1_1 as Y
    elif grammar=="1.0":
        import YappsStarParser_1_0 as Y
    elif grammar=="DDLm":
        import YappsStarParser_DDLm as Y
    if isinstance(filename,basestring):
        filestream = urlopen(filename)
    else:
        filestream = filename   #already opened for us
    my_uri = ""
    if hasattr(filestream,"geturl"): 
        my_uri = filestream.geturl()
    text = filestream.read()
    if isinstance(filename,basestring): #we opened it, we close it
        filestream.close()
    if not text:      # empty file, return empty block
        dest.set_uri(my_uri)
        return dest
    # we recognise ctrl-Z as end of file
    endoffile = text.find('\x1a')
    if endoffile >= 0: 
        text = text[:endoffile]
    split = string.split(text,'\n')
    if maxlength > 0:
        toolong = filter(lambda a:len(a)>maxlength,split)
        if toolong:
            pos = split.index(toolong[0])
            raise StarError( 'Line %d contains more than %d characters' % (pos+1,maxlength))
    try: 
        if scantype == 'standard':
            parser = Y.StarParser(Y.StarParserScanner(text))
        else:
            parser = Y.StarParser(Y.yappsrt.Scanner(None,[],text,scantype='flex'))
        proto_star = getattr(parser,"input")()
    except Y.yappsrt.SyntaxError:
        errorstring = 'Syntax error in input file: last value parsed was %s' % Y.lastval
        errorstring = errorstring + '\nParser status: %s' % `parser._scanner`
        raise StarError( errorstring)
    # duplication check on all blocks
    audit_result = map(lambda a:(a,proto_star[a].audit()),proto_star.keys())
    audit_result = filter(lambda a:len(a[1])>0,audit_result)
    if audit_result:
        raise StarError( 'Duplicate keys as follows: %s' % `audit_result`)
    proto_star.set_uri(my_uri)
    return proto_star

def get_dim(dataitem,current=0,packlen=0):
    zerotypes = [IntType, LongType, 
                    FloatType, StringType]
    if type(dataitem) in zerotypes:
        return current, packlen
    if not dataitem.__class__ == ().__class__ and \
       not dataitem.__class__ == [].__class__:
       return current, packlen
    elif len(dataitem)>0: 
    #    print "Get_dim: %d: %s" % (current,`dataitem`)
        return get_dim(dataitem[0],current+1,len(dataitem))
    else: return current+1,0
    


