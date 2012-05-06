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
import re
import StarFile
import sys
class CifLoopBlock(StarFile.LoopBlock):
    def __init__(self,data=(),dimension=0,**kwargs):
        self.loopclass = CifLoopBlock
        if dimension > 1: 
            raise CifError( 'Attempt to nest loops, loop level %d' % dimension)
        StarFile.LoopBlock.__init__(self,data,dimension=dimension,**kwargs)
        # self.__iter__ = self.recursive_iter

    def __iter__(self):
        return self.recursive_iter()
 
    def AddLoopItem(self,data,precheck=False):
        StarFile.LoopBlock.AddLoopItem(self,data,precheck,maxlength=75)

    def insert_loop(self,newloop,**kwargs):
        if newloop.dimension > 1:
            raise  CifError( 'Attempt to insert inner loop, loop level %d' % dimension)
        StarFile.LoopBlock.insert_loop(self,newloop,**kwargs)

class CifBlock(CifLoopBlock):
    def __init__(self,data = (), strict = 1, maxoutlength=2048,wraplength=80,overwrite=True,dimension=0):
        self.strict = strict
        CifLoopBlock.__init__(self,data=data,dimension=0,maxoutlength=maxoutlength,wraplength=wraplength,overwrite=overwrite)
        if isinstance(data,(StarFile.StarBlock,CifBlock)):
            self.saves = StarFile.BlockCollection(datasource=data["saves"],element_class=CifBlock,type_tag="save")
        else:
            self.saves = StarFile.BlockCollection(element_class=CifLoopBlock,type_tag="save")
        if self.strict:
            self.checklengths()
        self.dictionary = None

    def RemoveCifItem(self,itemname): 
        CifLoopBlock.RemoveLoopItem(self,itemname)

    def __getitem__(self,key):
        if key == "saves":
            return self.saves
        try:     
           rawitem = CifLoopBlock.__getitem__(self,key)
        except KeyError:
           if self.dictionary:
               # send the dictionary the required key and a pointer to us
               rawitem = self.dictionary.derive_item(key,self)
           else:
               raise KeyError, 'No such item: %s' % key
        # we now have an item, we can try to convert it to a number if that is appropriate
        if not self.dictionary or not self.dictionary.has_key(key): return rawitem
        return self.dictionary.change_type(key,rawitem)

    def __setitem__(self,key,value):
        if key == "saves":
            self.saves[key] = value
        else:
            self.AddCifItem((key,value))

    def clear(self):
        CifLoopBlock.clear(self)
        self.saves = StarFile.BlockCollection(element_class=CifLoopBlock,type_tag="save_")

    def copy(self):
        newblock = CifLoopBlock.copy(self)
        newblock.saves = self.saves.copy()
        return self.copy.im_class(newblock)   #catch inheritance

    def has_key(self,key):
        if key == "saves": return 1
        else: return CifLoopBlock.has_key(self,key)

    def __str__(self):
        retstr = ''
        for sb in self.saves.keys():
            retstr = retstr + '\nsave_%s\n\n' % sb
            self.saves[sb].SetOutputLength(self.wraplength,self.maxoutlength)
            retstr = retstr + str(self.saves[sb])
            retstr = retstr + '\nsave_\n\n'
        return retstr + CifLoopBlock.__str__(self)

    # this is not appropriate for save blocks.  Instead, the save block
    # should be accessed directly for update
     
    def update(self,adict):
        loopdone = []
        if not isinstance(adict,CifBlock):
            raise TypeError
        for key in adict.block.keys():
            self.AddCifItem((key,adict[key]))
        for aloop in adict.loops:
            self.insert_loop(aloop,audit=True)

    def AddCifItem(self,data):
        # we accept only tuples, strings and lists!!
        if not (isinstance(data[0],(StringType,TupleType,ListType))):
                  raise TypeError, 'Cif datanames are either a string, tuple or list'
        # single items passed straight through to underlying routine
        # we catch single item loops as well...
        if isinstance(data[0],StringType):
            if isinstance(data[1],(TupleType,ListType)) and not isinstance(data[1],(StarFile.StarList,StarFile.StarTuple)):
                CifLoopBlock.AddLoopItem(self,((data[0],),((data[1],))))
            else:
                CifLoopBlock.AddLoopItem(self,data)
            return
        # otherwise, we unpack one level and send along.  This is different
        # to the StarBlock behaviour, which assumes that any tuples imply an
        # inner loop.
        keyvals = zip(data[0],data[1])
        map(lambda a:CifLoopBlock.AddLoopItem(self,a),keyvals)

    def checklengths(self):
        toolong = filter(lambda a:len(a)>75, self.keys())
        outstring = ""
        for it in toolong: outstring += "\n" + it
        if toolong:
           raise CifError( 'Following data names too long:' + outstring)

    def loopnames(self):
        return map(lambda a:a.keys(),self.loops)

    def assign_dictionary(self,dic):
        if not dic.diclang=="DDLm":
            print "Warning: ignoring dictionary %s" % dic.dic_as_cif.my_uri
            return
        self.dictionary = dic

    def merge(self,new_block,mode="strict",match_att=[],match_function=None,nosaves=False,
                   rel_keys = []):
        # deal with save frames
        if not nosaves:
            self["saves"].merge(new_block["saves"],mode,match_att=match_att,
                                                        match_function=match_function)
        if mode == 'strict':
           for key in new_block.item_order: 
               if self.has_key(key) and key not in match_att:
                  raise CifError( "Identical keys %s in strict merge mode" % key)
               elif key not in match_att:           #no change otherwise
                  if isinstance(key,StringType):
                      self[key] = new_block[key] 
                  else:
                      self.insert_loop(key)
        elif mode == 'replace':
           newkeys = new_block.keys()
           for ma in match_att:
              try:
                   newkeys.remove(ma)        #don't touch the special ones
              except ValueError:
                   pass
           for key in new_block.item_order: 
                  if isinstance(key,StringType):
                      self[key] = new_block[key] 
                  else:
                      self.insert_loop(key)   #assume is a loop
        elif mode == 'overlay': 
           for attribute in new_block.keys():
               if attribute in match_att: continue      #ignore this one
               new_value = new_block[attribute]
               #non-looped items
               if isinstance(new_value,StringType):
                  self[attribute] = new_value 
           these_atts = self.keys()
           for newloop in new_block.loops:              
               newkeys = newloop.keys()
               # note that the following line determines packet item order
               overlaps = filter(lambda a: a in these_atts,newkeys)
               if len(overlaps)< len(newloop):#completely new loop
                  self.insert_loop(newloop)
               elif len(overlaps)==len(newloop):
                  # appending packets 
                  # print "In overlay merge mode, found extra packet items:"
                  # print `overlaps`
                  # get key position
                  loop_keys = filter(lambda a:a in rel_keys,overlaps)
                  try:
                     newkeypos = map(lambda a:newkeys.index(a),loop_keys)
                     newkeypos = newkeypos[0]      #one key per loop for now
                     loop_keys = loop_keys[0] 
                  except (ValueError,IndexError):
                     newkeypos = []
                  overlap_data = map(lambda a:listify(self[a]),overlaps) #old packet data
                  new_data = map(lambda a:new_block[a],overlaps) #new packet data
                  packet_data = transpose(overlap_data)
                  new_p_data = transpose(new_data)
                  # remove any packets for which the keys match between old and new; we
                  # make the arbitrary choice that the old data stays
                  if newkeypos:
                      # get matching values in new list
                      print "Old, new data:\n%s\n%s" % (`overlap_data[newkeypos]`,`new_data[newkeypos]`)
                      key_matches = filter(lambda a:a in overlap_data[newkeypos],new_data[newkeypos])
                      # filter out any new data with these key values
                      new_p_data = filter(lambda a:a[newkeypos] not in key_matches,new_p_data)
                      if new_p_data:
                          new_data = transpose(new_p_data)
                      else: new_data = []
                  # wipe out the old data and enter the new stuff
                  byebyeloop = self.GetLoop(overlaps[0])
                  # print "Removing '%s' with overlaps '%s'" % (`byebyeloop`,`overlaps`)
                  # Note that if, in the original dictionary, overlaps are not
                  # looped, GetLoop will return the block itself.  So we check
                  # for this case...
                  if byebyeloop != self:
                      self.remove_loop(byebyeloop)
                  self.AddCifItem(((overlaps,),(overlap_data,)))  #adding old packets
                  for pd in new_p_data:                             #adding new packets
                     if pd not in packet_data:
                        for i in range(len(overlaps)):
                            #don't do this at home; we are appending
                            #to something in place
                            self[overlaps[i]].append(pd[i]) 
                              

class CifFile(StarFile.StarFile):
    def __init__(self,datasource=None,strict=1,maxinlength=2048,maxoutlength=0,**kwargs):
        StarFile.StarFile.__init__(self,datasource=datasource,maxinlength=maxinlength,maxoutlength=maxoutlength,blocktype=CifBlock,**kwargs)
        self.strict = strict
        self.header_comment = \
"""#\\#CIF1.1
##########################################################################
#               Crystallographic Information Format file 
#               Produced by PyCifRW module
# 
#  This is a CIF file.  CIF has been adopted by the International
#  Union of Crystallography as the standard for data archiving and 
#  transmission.
#
#  For information on this file format, follow the CIF links at
#  http://www.iucr.org
##########################################################################
"""
    def NewBlock(self,blockname,*nkwargs,**kwargs):
       if len(blockname)>75:
           raise CifError , 'Blockname %s is longer than 75 characters' % blockname
       else:
           StarFile.StarFile.NewBlock(self,blockname,*nkwargs,**kwargs)


class CifError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return '\nCif Format error: '+ self.value 

class ValidCifError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return '\nCif Validity error: ' + self.value

class CifDic(StarFile.BlockCollection):
    def __init__(self,dic,do_minimum=False,grammar='1.1'):
        self.do_minimum = do_minimum
        self.dic_as_cif = dic
        self.template_cache = {}    #for DDLm imports
        self.ddlm_functions = {}    #for DDLm functions
        self.switch_numpy(False)    #no Numpy arrays returned 
        if isinstance(dic,StringType):
            self.dic_as_cif = CifFile(dic,grammar=grammar)
        (self.dicname,self.diclang,self.defdata) = self.dic_determine(self.dic_as_cif)
        StarFile.BlockCollection.__init__(self,element_class=CifBlock,datasource=self.defdata) 
        self.scopes_mandatory = {"dictionary":[],"category":[],"item":[]}
        self.scopes_naughty = {"dictionary":[],"category":[],"item":[]}
        # rename and expand out definitions using "_name" in DDL dictionaries
        if self.diclang == "DDL1":
            self.DDL1_normalise()   #this removes any non-definition entries
            self.ddl1_cat_load()
        elif self.diclang == "DDL2":
            self.DDL2_normalise()   #iron out some DDL2 tricky bits
        elif self.diclang == "DDLm":
            self.ddlm_normalise()
            self.ddlm_import()      #recursively calls this routine
            if not self.do_minimum:
                print "Doing full dictionary initialisation" 
                self.ddlm_parse_valid() #extract validity information from data block
                self.transform_drel()   #parse the drel functions
                self.add_drel_funcs()   #put the drel functions into the namespace
        self.add_category_info()
        # initialise type information
        self.typedic={}
        self.primdic = {}   #typecode<->primitive type translation
        self.add_type_info()
        self.item_validation_funs = [
            self.validate_item_type,
            self.validate_item_esd,
            self.validate_item_enum,   # functions which check conformance
            self.validate_enum_range,
            self.validate_looping]
        self.loop_validation_funs = [
            self.validate_loop_membership,
            self.validate_loop_key,
            self.validate_loop_references]    # functions checking loop values
        self.global_validation_funs = [
            self.validate_exclusion,
            self.validate_parent,
            self.validate_child,
            self.validate_dependents,
            self.validate_uniqueness] # where we need to look at other values
        self.block_validation_funs = [  # where only a full block will do
            self.validate_mandatory_category]
        self.global_remove_validation_funs = [
            self.validate_remove_parent_child] # removal is quicker with special checks
        self.optimize = False        # default value
        self.done_parents = []
        self.done_children = []
        self.done_keys = []
        # debug
        # j = open("dic_debug","w")
        # j.write(self.__str__())
        # j.close()

    def dic_determine(self,cifdic):
        if cifdic.has_key("on_this_dictionary"): 
            self.master_key = "on_this_dictionary"
            self.type_spec = "_type"
            self.enum_spec = "_enumeration"
            self.cat_spec = "_category"
            self.esd_spec = "_type_conditions"
            self.must_loop_spec = "_list"
            self.must_exist_spec = "_list_mandatory"
            self.list_ref_spec = "_list_reference"
            self.unique_spec = "_list_uniqueness"
            self.child_spec = "_list_link_child"
            self.parent_spec = "_list_link_parent"
            self.related_func = "_related_function"
            self.related_item = "_related_item"
            self.primitive_type = "_type"
            self.dep_spec = "xxx"
            self.cat_list = []   #to save searching all the time
            name = cifdic["on_this_dictionary"]["_dictionary_name"]
            version = cifdic["on_this_dictionary"]["_dictionary_version"]
            return (name+version,"DDL1",cifdic)
        elif len(cifdic.keys()) == 1:              # DDL2/DDLm
            self.master_key = cifdic.keys()[0]      
            name = cifdic[self.master_key]["_dictionary.title"]
            version = cifdic[self.master_key]["_dictionary.version"]
            if name != self.master_key:
                print "Warning: DDL2 blockname %s not equal to dictionary name %s" % (self.master_key,name)
            if cifdic[self.master_key].has_key("_dictionary.class"):   #DDLm
                self.unique_spec = "_category_key.generic"
                return(name+version,"DDLm",cifdic[self.master_key]["saves"]) 
            #otherwise DDL2
            self.type_spec = "_item_type.code" 
            self.enum_spec = "_item_enumeration.value"
            self.esd_spec = "_item_type_conditions.code"
            self.cat_spec = "_item.category_id" 
            self.loop_spec = "there_is_no_loop_spec!"
            self.must_loop_spec = "xxx"
            self.must_exist_spec = "_item.mandatory_code"
            self.child_spec = "_item_linked.child_name"
            self.parent_spec = "_item_linked.parent_name"
            self.related_func = "_item_related.function_code"
            self.related_item = "_item_related.related_name"
            self.unique_spec = "_category_key.name"
            self.list_ref_spec = "xxx"
            self.primitive_type = "_type"
            self.dep_spec = "_item_dependent.dependent_name"
            return (name+version,"DDL2",cifdic[self.master_key]["saves"])
        else:
            raise CifError, "Unable to determine dictionary DDL version"
        
    def DDL1_normalise(self):
        # add default type information in DDL2 style
        # initial types and constructs
        base_types = ["char","numb","null"]
        prim_types = base_types[:] 
        base_constructs = [".*",
            '(-?(([0-9]*[.][0-9]+)|([0-9]+)[.]?)([(][0-9]+[)])?([eEdD][+-]?[0-9]+)?)|\?|\.',
            "\"\" "]
        for key,value in self.dictionary.items():
           if value.has_key("_name"):
               real_name = value["_name"]
               if type(real_name) is ListType:        #looped values
                   for looped_name in real_name:
                       new_value = value.copy()
                       new_value["_name"] = looped_name  #only looped name
                       self.dictionary[looped_name] = new_value
               else: self.dictionary[real_name] = value
           # delete the old one
           del self.dictionary[key]
        # loop again to normalise the contents of each definition
        for key,value in self.dictionary.items():
           # deal with a missing _list, _type_conditions
           if not value.has_key("_list"): value["_list"] = 'no'
           if not value.has_key("_type_conditions"): value["_type_conditions"] = 'none'
           # deal with enumeration ranges
           if value.has_key("_enumeration_range"):
               max,min = self.getmaxmin(value["_enumeration_range"])
               if min == ".":
                   self.dictionary[key].AddLoopItem((("_item_range.maximum","_item_range.minimum"),((max,max),(max,min))))
               elif max == ".":
                   self.dictionary[key].AddLoopItem((("_item_range.maximum","_item_range.minimum"),((max,min),(min,min))))
               else:
                   self.dictionary[key].AddLoopItem((("_item_range.maximum","_item_range.minimum"),((max,max,min),(max,min,min))))
           #add any type construct information
           if value.has_key("_type_construct"):
               base_types.append(value["_name"]+"_type")   #ie dataname_type
               base_constructs.append(value["_type_construct"]+"$")
               prim_types.append(value["_type"])     #keep a record
               value["_type"] = base_types[-1]   #the new type name
               
        
        #make categories conform with ddl2
        #note that we must remove everything from the last underscore
           if value["_category"] == "category_overview":
                last_under = value["_name"].rindex("_")
                catid = value["_name"][1:last_under]
                value["_category.id"] = catid  #remove square bracks
                if catid not in self.cat_list: self.cat_list.append(catid) 
        # we now add any missing categories before filling in the rest of the
        # information
        for key,value in self.dictionary.items():
            if self[key].has_key("_category"):
                if self[key]["_category"] not in self.cat_list:
                    # rogue category, add it in
                    newcat = self[key]["_category"]
                    fake_name = "_" + newcat + "_[]" 
                    newcatdata = CifBlock()
                    newcatdata["_category"] = "category_overview"
                    newcatdata["_category.id"] = newcat
                    newcatdata["_type"] = "null"
                    self[fake_name] = newcatdata
                    self.cat_list.append(newcat)
        # write out the type information in DDL2 style
        self.dic_as_cif[self.master_key].AddLoopItem((
            ("_item_type_list.code","_item_type_list.construct",
              "_item_type_list.primitive_code"),
            (base_types,base_constructs,prim_types)
            ))
     
    def DDL2_normalise(self):
       listed_defs = filter(lambda a:isinstance(self[a].get('_item.name'),ListType),self.keys()) 
       # now filter out all the single element lists!
       dodgy_defs = filter(lambda a:len(self[a]['_item.name']) > 1, listed_defs)
       for item_def in dodgy_defs:
          # print "DDL2 norm: processing %s" % item_def
          thisdef = self[item_def]
          packet_no = thisdef['_item.name'].index(item_def)
          realcat = thisdef['_item.category_id'][packet_no] 
          realmand = thisdef['_item.mandatory_code'][packet_no]
          # first add in all the missing categories
          # we don't replace the entry in the list corresponding to the
          # current item, as that would wipe out the information we want
          for child_no in range(len(thisdef['_item.name'])):
              if child_no == packet_no: continue
              child_name = thisdef['_item.name'][child_no]
              child_cat = thisdef['_item.category_id'][child_no]
              child_mand = thisdef['_item.mandatory_code'][child_no]
              if not self.has_key(child_name):
                  self[child_name] = CifBlock()
                  self[child_name]['_item.name'] = child_name
              self[child_name]['_item.category_id'] = child_cat
              self[child_name]['_item.mandatory_code'] = child_mand
          self[item_def]['_item.name'] = item_def
          self[item_def]['_item.category_id'] = realcat
          self[item_def]['_item.mandatory_code'] = realmand
       # go through any _item_linked tables
       dodgy_defs = filter(lambda a:isinstance(self[a].get('_item_linked.child_name'),ListType),self.keys()) 
       dodgy_defs = filter(lambda a:len(self[a]['_item_linked.child_name']) > 1, dodgy_defs)
       for item_def in dodgy_defs:
          thisdef = self[item_def]
          child_list = thisdef.get('_item_linked.child_name',[])
          parents = thisdef.get('_item_linked.parent_name',[])
          # zap the parents, they will confuse us!!
          del thisdef['_item_linked.parent_name']
          if isinstance(child_list,StringType):
              self[child_list]['_item_linked.parent_name'] = parents
              self[parents]['_item_linked.child_name'] = child_list 
          else:
              # for each parent, find the list of children.
              family = map(None,parents,child_list)
              notmychildren = family
              while len(notmychildren):
                  # get all children of first entry
                  mychildren = filter(lambda a:a[0]==notmychildren[0][0],family)
                  # print "Parent %s: %d children" % (notmychildren[0][0],len(mychildren))
                  for parent,child in mychildren:   #parent is the same for all
                      self[child]['_item_linked.parent_name'] = parent
                  # put all the children into the parent
                  try:
                      del self[mychildren[0][0]]['_item_linked.child_name']
                  except ValueError: pass
                  self[mychildren[0][0]]['_item_linked.child_name'] = map(lambda a:a[1],mychildren)
                  # now make a new,smaller list
                  notmychildren = filter(lambda a:a[0]!=mychildren[0][0],notmychildren)
       # now flatten any single element lists
       single_defs = filter(lambda a:len(self[a]['_item.name'])==1,listed_defs)
       for flat_def in single_defs:
           flat_keys = self[flat_def].GetLoop('_item.name').keys()
           for flat_key in flat_keys: self[flat_def][flat_key] = self[flat_def][flat_key][0]
       # now deal with the multiple lists
       # next we do aliases
       all_aliases = filter(lambda a:self[a].has_key('_item_aliases.alias_name'),self.keys()) 
       for aliased in all_aliases:
          my_aliases = listify(self[aliased]['_item_aliases.alias_name'])
          for alias in my_aliases:
              self[alias] = self[aliased].copy()   #we are going to delete stuff...
              del self[alias]["_item_aliases.alias_name"]
 
    def ddlm_normalise(self):
        for key,value in self.dictionary.items():
           if value.has_key("_name.category_id"):
               real_name = "_" + value["_name.category_id"] + "." + value["_name.object_id"]
               self[real_name] = value
               # delete the old one
               del self[key]
        
    def ddlm_parse_valid(self):
        if not self.dic_as_cif[self.master_key].has_key("_dictionary_valid.scope"):
            return
        for scope_pack in self.dic_as_cif[self.master_key].GetLoop("_dictionary_valid.scope"):
            scope = getattr(scope_pack,"_dictionary_valid.scope")
            valid_info = getattr(scope_pack,"_dictionary_valid.attributes")
            valid_info = valid_info.split()
            for i in range(0,len(valid_info),2): 
                if valid_info[i]=="+":
                   self.scopes_mandatory[scope.lower()].append(valid_info[i+1].lower())
                elif valid_info[i]=="!":
                   self.scopes_naughty[scope.lower()].append(valid_info[i+1].lower())

    def ddlm_import(self):
        import urllib
        #first check the outermost datablocks.  Note we expect our dREL
        #machinery to create _import_list.id only if the individual values are available
        #For this to happen, we need the ddl.dic to have been assigned
        try:
            to_be_imported = self.dic_as_cif[self.master_key]["_import_list.id"]
        except KeyError:
            pass
        else:
            # deal with foreshortened import blocks
            for import_target in to_be_imported:
                if len(import_target)==3:                     #default values have been left off
                    import_target.append('Exit')
                    import_target.append('Exit')
            for scope,dict_block,file_loc,on_dupl,on_miss in to_be_imported: 
                scope = scope.lower()                         #work around capitalisation in draft dics
                if scope == 'att' or scope == 'sta' or scope == 'val':
                    print 'Improper import directive at top level in %s: ignored' % self.master.key
                    continue 
                # resolve URI  
                full_uri = self.resolve_path(file_loc)
                dic_as_cif = CifFile(urllib.urlopen(full_uri),grammar="DDLm")
                import_from = CifDic(dic_as_cif,do_minimum=True)  #this will recurse internal imports
                # and now merge these definitions
                if scope == "dic":
                    self.get_whole_dict(import_from,on_dupl,on_miss)
                elif scope=="cat":
                    self.get_one_cat(import_from,dict_block,on_dupl,on_miss)
                elif scope=="grp":
                    self.get_one_cat_with_children(import_from,dict_block,on_dupl,on_miss)
                elif scope=="itm":  #not clear what to do if category is missing
                    self.add_one_defn(import_from,dict_block,on_dupl,on_miss)
            # it will never happen again...
            del self.dic_as_cif[self.master_key]["_import_list.id"]
        # next we resolve per-definition imports
        for one_def in self.keys():
            try: 
                to_be_imported = self[one_def]["_import_list.id"]
            except KeyError:
                pass
            else:
                if len(to_be_imported) == 5 and len(to_be_imported[0])!=5:
                    #catch an error in earlier versions of the dictionaries where
                    #the outer brackets were missing
                    to_be_imported = [to_be_imported]
                # deal with foreshortened import blocks
                for import_target in to_be_imported:
                    if len(import_target)==3:                     #default values have been left off
                        import_target.append('Exit')
                        import_target.append('Exit')
                for scope,block,file_loc,on_dupl,on_miss in to_be_imported: 
                    scope = scope.lower()                         #work around capitalisation in draft dics
                    if scope == 'dic' or scope == 'cat' or scope == 'grp' or scope == "itm":
                        print 'Improper import directive at definition level in %s: ignored' % self.master.key
                        continue 
                    full_uri = self.resolve_path(file_loc)
                    if full_uri not in self.template_cache:
                        dic_as_cif = CifFile(urllib.urlopen(full_uri),grammar="DDLm")
                        self.template_cache[full_uri] = CifDic(dic_as_cif,do_minimum=True)  #this will recurse internal imports
                        print 'Added %s to cached dictionaries' % full_uri
                    import_from = self.template_cache[full_uri]
                    if scope == 'att': 
                        self.import_attributes(one_def,import_from,block,on_dupl,on_miss) 
                    elif scope == 'sta': 
                        self.import_loop(one_def,import_from,block,'_enumeration_set.state',on_miss)
                    elif scope == 'val': 
                        self.import_loop(one_def,import_from,block,'_enumeration_default.value',on_miss)
                    else:
                        raise CifError, "Unrecognised import scope %s" % scope
                # remove the import attribute 
                del self[one_def]["_import_list.id"]   
                    
    def resolve_path(self,file_loc):
        import urlparse
        url_comps = urlparse.urlparse(file_loc)
        if url_comps[0]: return file_loc    #already full URI
        new_url = urlparse.urljoin(self.dic_as_cif.my_uri,file_loc)
        print "Transformed %s to %s for import " % (file_loc,new_url)
        return new_url
        
    def get_whole_dict(self,source_dict,on_dupl,on_miss):
        print "Cat_map: `%s`" % source_dict.cat_map.values()
        for source_cat in source_dict.cat_map.values():
            self.get_one_cat(source_dict,source_cat,on_dupl,on_miss)
       
    def get_one_cat(self,source_dict,source_cat,on_dupl,on_miss):
        ext_cat = source_dict.get(source_cat,"")
        this_cat = self.get(source_cat,"")
        print "Adding category %s" % source_cat
        if not ext_cat:
            if on_miss == "Ignore":
               pass
            else:
               raise CifError, "Missing category %s" % source_cat 
        else:
            all_ext_defns = source_dict.keys()
            cat_list = filter(lambda a:source_dict[a].get("_name.category_id","").lower()==source_cat.lower(),
                               all_ext_defns) 
            print "Items: %s" % `cat_list`
            if this_cat:     # The category block itself is duplicated
                if on_dupl=="Ignore":
                    pass
                elif on_dupl == "Exit":
                    raise CifError, "Duplicate category %s" % source_cat
                else: 
                    self[source_cat] = ext_cat
            else:
                self[source_cat] = ext_cat
            # now do all member definitions
            for cat_defn in cat_list:
                self.add_one_defn(source_dict,cat_defn,on_dupl)

    def add_one_defn(self,source_dict,cat_defn,on_dupl):
        if self.has_key(cat_defn):
           if on_dupl == "Ignore": pass
           elif on_dupl == "Exit": 
                   raise CifError, "Duplicate definition %s" % cat_defn
           else: self[cat_defn] = source_dict[cat_defn]
        else: self[cat_defn] = source_dict[cat_defn]
        print "    "+cat_defn
        
    def get_one_cat_with_children(self,source_dict,source_cat,on_dupl,on_miss):
        self.get_one_cat(source_dict,source_cat,on_dupl,on_miss)
        child_cats = filter(lambda a:source_dict[a]["_category.parent_id"]==source_dict[source_cat]["_definition.id"],source_dict.cat_map.values())
        for child_cat in child_cats: self.get_one_cat(source_dict,child_cat,on_dupl,on_miss) 

    def import_attributes(self,mykey,source_dict,source_def,on_dupl,on_miss):
        # process missing 
        if not source_dict.has_key(source_def): 
            if on_miss == 'Exit':
                raise CifError, 'Missing definition for import %s' % source_def
            else: return          #nothing else to do
        # now do the import
        print 'Adding attributes from %s to %s' % (source_def,mykey)
        self[mykey].merge(source_dict[source_def],mode='replace',match_att= \
              ['_definition.id','_name.category_id','_name.object_id'])

    def import_loop(self,mykey,source_dict,source_def,loop_name,on_miss):
        # process imssing
        if not source_dict.has_key(source_def): 
            if on_miss == 'Exit':
                raise CifError, 'Missing definition for import %s' % source_def
            else: return          #nothing else to do
        print 'Adding %s attributes from %s to %s' % (loop_name,source_def,mykey)
        state_loop = source_dict[source_def].GetLoop(loop_name)
        self[mykey].insert_loop(state_loop) 
       

    def ddl1_cat_load(self):
        deflist = self.keys()       #slight optimization
        cat_mand_dic = {}
        cat_unique_dic = {}
        # a function to extract any necessary information from each definition
        def get_cat_info(single_def):
            if self[single_def].get(self.must_exist_spec)=='yes':
                thiscat = self[single_def]["_category"]
                curval = cat_mand_dic.get(thiscat,[])
                curval.append(single_def)
                cat_mand_dic[thiscat] = curval
            # now the unique items...
            # cif_core.dic throws us a curly one: the value of list_uniqueness is
            # not the same as the defined item for publ_body_label, so we have
            # to collect both together.  We assume a non-listed entry, which
            # is true for all current (May 2005) ddl1 dictionaries.
            if self[single_def].get(self.unique_spec,None)!=None:
                thiscat = self[single_def]["_category"]
                new_unique = self[single_def][self.unique_spec]
                uis = cat_unique_dic.get(thiscat,[])
                if single_def not in uis: uis.append(single_def)
                if new_unique not in uis: uis.append(new_unique)
                cat_unique_dic[thiscat] = uis
            
        map(get_cat_info,deflist)       # apply the above function
        for cat in cat_mand_dic.keys():
            cat_entry = self.get_ddl1_entry(cat)
            self[cat_entry]["_category_mandatory.name"] = cat_mand_dic[cat]
        for cat in cat_unique_dic.keys():
            cat_entry = self.get_ddl1_entry(cat)
            self[cat_entry]["_category_key.name"] = cat_unique_dic[cat]

    # A helper function get find the entry corresponding to a given category name:
    # yes, in DDL1 the actual name is different in the category block due to the
    # addition of square brackets which may or may not contain stuff.

    def get_ddl1_entry(self,cat_name):
        chop_len = len(cat_name) 
        possibles = filter(lambda a:a[1:chop_len+3]==cat_name+"_[",self.keys())
        if len(possibles) > 1 or possibles == []:
            raise ValidCifError, "Category name %s can't be matched to category entry" % cat_name
        else:
            return possibles[0]

    def add_type_info(self):
        if self.dic_as_cif[self.master_key].has_key("_item_type_list.construct"): 
            types = self.dic_as_cif[self.master_key]["_item_type_list.code"]
            prim_types = self.dic_as_cif[self.master_key]["_item_type_list.primitive_code"]
            constructs = map(lambda a: a + "$", self.dic_as_cif[self.master_key]["_item_type_list.construct"])
            # add in \r wherever we see \n, and change \{ to \\{
            def regex_fiddle(mm_regex):
                brack_match = r"((.*\[.+)(\\{)(.*\].*))" 
                ret_match = r"((.*\[.+)(\\n)(.*\].*))" 
                fixed_regexp = mm_regex[:]  #copy
                # fix the brackets
                bm = re.match(brack_match,mm_regex)
                if bm != None: 
                    fixed_regexp = bm.expand(r"\2\\\\{\4")
                # fix missing \r
                rm = re.match(ret_match,fixed_regexp)
                if rm != None:
                    fixed_regexp = rm.expand(r"\2\3\\r\4")    
                #print "Regexp %s becomes %s" % (mm_regex,fixed_regexp)
                return fixed_regexp
            constructs = map(regex_fiddle,constructs)
            packed_up = map(None,types,constructs)
            for typecode,construct in packed_up:
                self.typedic[typecode] = re.compile(construct,re.MULTILINE|re.DOTALL)
            # now make a primitive <-> type construct mapping
            packed_up = map(None,types,prim_types)
            for typecode,primtype in packed_up:
                self.primdic[typecode] = primtype

    def add_category_info(self):
        if self.diclang == "DDLm":
            categories = filter(lambda a:self[a].get("_definition.scope","Item")=="Category",self.keys())
            category_ids = map(lambda a:self[a]["_definition.id"],categories)


        else:
            categories = filter(lambda a:self[a].has_key("_category.id"),self.keys())
            # get the category id
            category_ids = map(lambda a:self[a]["_category.id"],categories)

        # match ids and entries in the dictionary
        catpairs = map(None,category_ids,categories)
        self.cat_map = {}
        for catid,cat in catpairs:self.cat_map[catid] = cat

    def names_in_cat(self,cat):
        nameblocks = filter(lambda a:self[a].get("_name.category_id","").lower()
                             ==cat.lower(),self.keys())
        return map(lambda a:"_" + self[a]["_name.category_id"]+"." + self[a]["_name.object_id"],nameblocks)
        
    def get_key_pack(self,category,value,data):
        keyname = self[category][self.unique_spec]
        onepack = data.GetPackKey(keyname,value)
        return onepack
     
    def get_number_with_esd(numstring):
        import string
        numb_re = '((-?(([0-9]*[.]([0-9]+))|([0-9]+)[.]?))([(][0-9]+[)])?([eEdD][+-]?[0-9]+)?)|(\?)|(\.)' 
        our_match = re.match(numb_re,numstring)
        if our_match:
            a,base_num,b,c,dad,dbd,esd,exp,q,dot = our_match.groups()
        #    print "Debug: %s -> %s" % (numstring, `our_match.groups()`)
        else:
            return None,None
        if dot or q: return None,None     #a dot or question mark
        if exp:          #has exponent 
           exp = string.replace(exp,"d","e")     # mop up old fashioned numbers
           exp = string.replace(exp,"D","e")
           base_num = base_num + exp
        #print "Debug: have %s for base_num from %s" % (base_num,numstring)
        base_num = float(base_num)
        # work out esd, if present.
        if esd:
            esd = float(esd[1:-1])    # no brackets
            if dad:                   # decimal point + digits
                esd = esd * (10 ** (-1* len(dad)))
            if exp:
                esd = esd * (10 ** (float(exp[1:])))
        return base_num,esd

    def getmaxmin(self,rangeexp):
        regexp = '(-?(([0-9]*[.]([0-9]+))|([0-9]+)[.]?)([eEdD][+-]?[0-9]+)?)*' 
        regexp = regexp + ":" + regexp
        regexp = re.match(regexp,rangeexp)
        try:
            minimum = regexp.group(1)
            maximum = regexp.group(7)
        except AttributeError:
            print "Can't match %s" % rangeexp
        if minimum == None: minimum = "." 
        else: minimum = float(minimum)
        if maximum == None: maximum = "." 
        else: maximum = float(maximum)
        return maximum,minimum

    def transform_drel(self):
        import drel_yacc
        parser = drel_yacc.parser
        my_namespace = self.keys()
        my_namespace = dict(map(None,my_namespace,my_namespace)) 
        parser.loopable_cats = filter(lambda a:self[a].get("_definition.class","Datum")=="List",self.keys())
        parser.loopable_cats = map(lambda a:self[a]["_definition.id"],parser.loopable_cats) 
        parser.listable_items = filter(lambda a:"*" in self[a].get("_type.dimension",""),self.keys()) 
        derivable_list = filter(lambda a:self[a].has_key("_method.expression") and self[a].get("_definition.scope","")!='Category' and self[a].get("_name.category_id","")!= "function",self.keys())
        for derivable in derivable_list:
            parser.target_id = derivable
            # reset the list of visible names for parser
            parser.special_id = [my_namespace]
            # reset list of looped with statements
            parser.withtable = {}
            print "Target id: %s" % derivable
            drel_expr = self[derivable]["_method.expression"]
            if isinstance(drel_expr,ListType):
               drel_expr = drel_expr[0] 
            print "Transforming %s" % drel_expr
            # List categories are treated differently...
            pyth_meth = parser.parse(drel_expr,debug=True)
            self[derivable]["_loop_categories"] = pyth_meth[1].keys()
            self[derivable]["_method.expression"] = drel_yacc.make_func(pyth_meth,"pyfunc",None) 
            print "Final result:\n " + self[derivable]["_method.expression"]

    def add_drel_funcs(self):
        import drel_yacc
        funclist = filter(lambda a:self[a].get("_name.category_id","")=='function',self.keys())
        funcnames = map(lambda a:self[a]["_name.object_id"],funclist)
        funcbodys = map(lambda a:self[a]["_method.expression"],funclist)
        # create executable python code...
        parser = drel_yacc.parser
        for funcname,funcbody in zip(funcnames,funcbodys):
            parser.target_id = funcname
            parser.special_id = [{}]   #first element is always global namespace of dictionary
            parser.withtable = {}
            res,ww = parser.parse(funcbody[0])
            print 'dREL library function ->\n' + res
            global_table = globals()
            global_table.update(self.ddlm_functions)
            exec res in global_table    #add to namespace
        print "All functions -> " + `self.ddlm_functions`
 
    def switch_numpy(self,to_val):
        if to_val:
            self.recursive_numerify = self.numpy_numerify
        else:
            self.recursive_numerify = self.normal_numerify

    def derive_item(self,key,cifdata,store_value = False):
        # store any default value in case we have a problem
        def_val = self[key].get("_enumeration.default","")
        def_index_val = self[key].get("_enumeration.def_index_id","")
        the_func = self[key].get("_method.expression","")  
        if def_val and not the_func : return def_val
        if def_index_val and not the_func:            #derive a default value
            index_vals = self[key]["_enumeration_default.index"]
            val_to_index = cifdata[def_index_val]     #what we are keying on
            # Handle loops
            if isinstance(val_to_index,ListType):
                keypos = map(lambda a:index_vals.index(a),val_to_index)
                result = map(lambda a:self[key]["_enumeration_default.value"][a] ,keypos)
            else:
                keypos = index_vals.index(val_to_index)   #value error if no such value available
                result = self[key]["_enumeration_default.value"][keypos] 
            print "Indexed on %s to get %s for %s" % (def_index_val,`result`,`val_to_index`)
            return result
        # read it in
        the_category = self[key]["_name.category_id"]
        the_type = self[the_category]["_definition.class"]
        global_table = globals()
        global_table.update(self.ddlm_functions)
        exec the_func in global_table,locals() #will access dREL functions, puts "pyfunc" in scope
        print 'Executing following function'
        print the_func
        print 'With following loop categories:' + `self[key].get("_loop_categories","")`
        # print 'in following global environment: ' + `global_table` 
        if self[key].get("_loop_categories",""): 
           loop_category = self[key]["_loop_categories"][0]
           loop_names = self.names_in_cat(loop_category)
           no_of_packs = len(cifdata[loop_names[0]])
           packlist = []
           for pack_index in range(no_of_packs):
              packlist.append(pyfunc(self,cifdata,pack_index))
        # now try to insert the new information into the right place
        # find if items of this category already appear...
           if store_value: 
               cat_names = filter(lambda a:self[a].get["_name.category_id",None]==the_category,self.keys())
               has_cat_names = filter(lambda a:cifdata.has_key(a),cat_names)
               if len(has_cat_names)>0:
                  target_loop = cifdata.GetLoop(has_cat_names[0])
                  target_loop[key] = packlist      #lengths must match or else!!
               else:
                  cifdata[key] = packlist  
           return packlist
        else:              # No looped categories
           return pyfunc(self,cifdata)

    def change_type(self,itemname,inval):
        import numpy
        # we need to iterate over the structure description.  For now we deal only with
        # Single and Array containers, with types that are a simple specification
        item_type = self[itemname]["_type.contents"]
        item_container = self[itemname]["_type.container"]
        isnumeric = (item_type == "Real" or \
           item_type == "Float" or \
           item_type == "Count" or \
           item_type == "Integer" or \
           item_type == "Digit")
        if not isnumeric: return inval   # we don't attempt any changes
        # even for a 'Single' container, it may be looped 
        # print 'Changing type for %s' % `inval`
        if StarFile.get_dim(inval)[0] == 0:
            if item_container == 'Single': return float_with_esd(inval)
            if item_container == 'Array':
                return self.recursive_numerify(inval) 
        else:
            if item_container == 'Single': return map(float_with_esd,inval)
            if item_container == 'Array': return map(self.recursive_numerify,inval)
     
    # A utility function to recursively make all atomic values numeric
    # All embedded values will be either StarTuples or StarLists
    def normal_numerify(self,valarray): 
        # print 'Recursive evaluation of %s' % `valarray`
        if isinstance(valarray,StarFile.StarTuple):
            return StarFile.StarTuple(map(self.recursive_numerify,valarray))
        if isinstance(valarray,StarFile.StarList):
            return StarFile.StarList(map(self.recursive_numerify,valarray))
        if isinstance(valarray,(StringType,IntType,LongType)):
            return float_with_esd(valarray)
        else:
            return valarray    #assume is OK
                
    # Identical to the above except that a numpy array is returned.  We
    # do the normal_numerify call in order to perform the float conversion.
    #
    def numpy_numerify(self,valarray):
        import numpy
        return numpy.array(self.normal_numerify(valarray))

    def validate_item_type(self,item_name,item_value):
        def mymatch(m,a):  
            res = m.match(a)
            if res != None: return res.group() 
            else: return ""
        target_type = self[item_name].get(self.type_spec) 
        if target_type == None:          # e.g. a category definition
            return {"result":True}                  # not restricted in any way
        matchexpr = self.typedic[target_type]
        item_values = listify(item_value)
        #for item in item_values:
            #print "Type match " + item_name + " " + item + ":",
        #skip dots and question marks
        check_all = filter(lambda a: a !="." and a != "?",item_values)
        check_all = filter(lambda a: mymatch(matchexpr,a) != a, check_all)
        if len(check_all)>0: return {"result":False,"bad_values":check_all}
        else: return {"result":True}

    def validate_item_esd(self,item_name,item_value):
        if self[item_name].get(self.primitive_type) != 'numb':
            return {"result":None}
        can_esd = self[item_name].get(self.esd_spec,"none") == "esd" 
        if can_esd: return {"result":True}         #must be OK!
        item_values = listify(item_value)
        check_all = filter(lambda a: get_number_with_esd(a)[1] != None, item_values)
        if len(check_all)>0: return {"result":False,"bad_values":check_all}
        return {"result":True}

    def validate_enum_range(self,item_name,item_value):
        if not self[item_name].has_key("_item_range.minimum") and \
           not self[item_name].has_key("_item_range.maximum"):
            return {"result":None}
        minvals = self[item_name].get("_item_range.minimum",default = ["."])
        maxvals = self[item_name].get("_item_range.maximum",default = ["."])
        def makefloat(a):
            if a == ".": return a
            else: return float(a)
        maxvals = map(makefloat, maxvals)
        minvals = map(makefloat, minvals)
        rangelist = map(None,minvals,maxvals)
        item_values = listify(item_value)
        def map_check(rangelist,item_value):
            if item_value == "?" or item_value == ".": return True
            iv,esd = get_number_with_esd(item_value)
            if iv==None: return None  #shouldn't happen as is numb type
            for lower,upper in rangelist:
                #check the minima
                if lower == ".": lower = iv - 1
                if upper == ".": upper = iv + 1
                if iv > lower and iv < upper: return True
                if upper == lower and iv == upper: return True
            # debug
            # print "Value %s fails range check %d < x < %d" % (item_value,lower,upper)
            return False
        check_all = filter(lambda a,b=rangelist: map_check(b,a) != True, item_values)
        if len(check_all)>0: return {"result":False,"bad_values":check_all}
        else: return {"result":True}
                
    def validate_item_enum(self,item_name,item_value):
        try: 
            enum_list = self[item_name][self.enum_spec][:]
        except KeyError:
            return {"result":None}
        enum_list.append(".")   #default value
        enum_list.append("?")   #unknown
        item_values = listify(item_value)
        #print "Enum check: %s in %s" % (`item_values`,`enum_list`)
        check_all = filter(lambda a: a not in enum_list,item_values)
        if len(check_all)>0: return {"result":False,"bad_values":check_all}
        else: return {"result":True}

    def validate_looping(self,item_name,item_value):
        try:
            must_loop = self[item_name][self.must_loop_spec]
        except KeyError:
            return {"result":None}
        if must_loop == 'yes' and isinstance(item_value,StringType): # not looped
            return {"result":False}      #this could be triggered
        if must_loop == 'no' and not isinstance(item_value,StringType): 
            return {"result":False}
        return {"result":True}


    def validate_loop_membership(self,loop_names):
        try:
            categories = map(lambda a:self[a][self.cat_spec],loop_names)
        except KeyError:       #category is mandatory
            raise ValidCifError( "%s missing from dictionary %s for item in loop containing %s" % (self.cat_spec,self.dicname,loop_names[0]))
        bad_items =  filter(lambda a:a != categories[0],categories)
        if len(bad_items)>0:
            return {"result":False,"bad_items":bad_items}
        else: return {"result":True}

    def validate_loop_key(self,loop_names):
        category = self[loop_names[0]][self.cat_spec]
        # find any unique values which must be present 
        entry_name = self.cat_map[category]
        key_spec = self[entry_name].get("_category_mandatory.name",[])
        for names_to_check in key_spec:
            if isinstance(names_to_check,StringType):   #only one
                names_to_check = [names_to_check]
            for loop_key in names_to_check:
                if loop_key not in loop_names: 
                    #is this one of those dang implicit items?
                    if self[loop_key].get(self.must_exist_spec,None) == "implicit":
                        continue          #it is virtually there...
                    alternates = self.get_alternates(loop_key)
                    if alternates == []: 
                        return {"result":False,"bad_items":loop_key}
                    for alt_names in alternates:
                        alt = filter(lambda a:a in loop_names,alt_names)
                        if len(alt) == 0: 
                            return {"result":False,"bad_items":loop_key}  # no alternates   
        return {"result":True}
        
    def validate_loop_references(self,loop_names):
        must_haves = map(lambda a:self[a].get(self.list_ref_spec,None),loop_names)
        must_haves = filter(lambda a:a != None,must_haves)
        # build a flat list.  For efficiency we don't remove duplicates,as
        # we expect no more than the order of 10 or 20 looped names.
        def flat_func(a,b): 
            if isinstance(b,StringType): 
               a.append(b)       #single name
            else:
               a.extend(b)       #list of names
            return a
        flat_mh = reduce(flat_func,must_haves,[])
        group_mh = filter(lambda a:a[-1]=="_",flat_mh)
        single_mh = filter(lambda a:a[-1]!="_",flat_mh)
        res = filter(lambda a: a not in loop_names,single_mh)
        def check_gr(s_item, name_list):
            nl = map(lambda a:a[:len(s_item)],name_list)
            if s_item in nl: return True
            return False
        res_g = filter(lambda a:check_gr(a,loop_names),group_mh)
        if len(res) == 0 and len(res_g) == 0: return {"result":True}
        # construct alternate list
        alternates = map(lambda a: (a,self.get_alternates(a)),res)
        alternates = filter(lambda a:a[1] != [], alternates)
        # next two lines purely for error reporting
        missing_alts = filter(lambda a: a[1] == [], alternates)
        missing_alts = map(lambda a:a[0],missing_alts)
        if len(alternates) != len(res): 
           return {"result":False,"bad_items":missing_alts}   #short cut; at least one
                                                       #doesn't have an altern
        #loop over alternates
        for orig_name,alt_names in alternates:
             alt = filter(lambda a:a in loop_names,alt_names)
             if len(alt) == 0: return {"result":False,"bad_items":orig_name}# no alternates   
        return {"result":True}        #found alternates
             
    def get_alternates(self,main_name,exclusive_only=False):
        alternates = self[main_name].get(self.related_func,None)
        alt_names = []
        if alternates != None: 
            alt_names =  self[main_name].get(self.related_item,None)
            if isinstance(alt_names,StringType): 
                alt_names = [alt_names]
                alternates = [alternates]
            together = map(None,alt_names,alternates)
            if exclusive_only:
                alt_names = filter(lambda a:a[1]=="alternate_exclusive" \
                                             or a[1]=="replace", together)
            else:
                alt_names = filter(lambda a:a[1]=="alternate" or a[1]=="replace",together)
            alt_names = map(lambda a:a[0],alt_names)
        # now do the alias thing
        alias_names = listify(self[main_name].get("_item_aliases.alias_name",[]))
        alt_names.extend(alias_names)
        # print "Alternates for %s: %s" % (main_name,`alt_names`)
        return alt_names
        

    def validate_exclusion(self,item_name,item_value,whole_block,provisional_items={},globals={}):
       alternates = map(lambda a:a.lower(),self.get_alternates(item_name,exclusive_only=True))
       item_name_list = map(lambda a:a.lower(),whole_block.keys())
       item_name_list.extend(map(lambda a:a.lower(),provisional_items.keys()))
       item_name_list.extend(map(lambda a:a.lower(),globals.keys()))
       bad = filter(lambda a:a in item_name_list,alternates)
       if len(bad)>0:
           print "Bad: %s, alternates %s" % (`bad`,`alternates`)
           return {"result":False,"bad_items":bad}
       else: return {"result":True}

    # validate that parent exists and contains matching values
    def validate_parent(self,item_name,item_value,whole_block,provisional_items={},globals={}):
        parent_item = self[item_name].get(self.parent_spec)
        if not parent_item: return {"result":None}   #no parent specified
        if isinstance(parent_item,ListType): 
            parent_item = parent_item[0]
        if self.optimize:
            if parent_item in self.done_parents:
                return {"result":None}
            else: 
                self.done_parents.append(parent_item)
                print "Done parents %s" % `self.done_parents`
        # initialise parent/child values
        if isinstance(item_value,StringType):
            child_values = [item_value]
        else: child_values = item_value[:]    #copy for safety
        # track down the parent
        # print "Looking for %s parent item %s in %s" % (item_name,parent_item,`whole_block`)
        # if globals contains the parent values, we are doing a DDL2 dictionary, and so 
        # we have collected all parent values into the global block - so no need to search
        # for them elsewhere. 
        # print "Looking for %s" % `parent_item`
        parent_values = globals.get(parent_item)
        if not parent_values:
            parent_values = provisional_items.get(parent_item,whole_block.get(parent_item))
        if not parent_values:  
            # go for alternates
            namespace = whole_block.keys()
            namespace.extend(provisional_items.keys())
            namespace.extend(globals.keys())
            alt_names = filter_present(self.get_alternates(parent_item),namespace)
            if len(alt_names) == 0:
                if len(filter(lambda a:a != "." and a != "?",child_values))>0:
                    return {"result":False,"parent":parent_item}#no parent available -> error
                else:
                    return {"result":None}       #maybe True is more appropriate??
            parent_item = alt_names[0]           #should never be more than one?? 
            parent_values = provisional_items.get(parent_item,whole_block.get(parent_item))
            if not parent_values:   # check global block
                parent_values = globals.get(parent_item)
        if isinstance(parent_values,StringType):
            parent_values = [parent_values]   
        #print "Checking parent %s against %s, values %s/%s" % (parent_item,
        #                                          item_name,`parent_values`,`child_values`)
        missing = self.check_parent_child(parent_values,child_values)
        if len(missing) > 0:
            return {"result":False,"bad_values":missing,"parent":parent_item}
        return {"result":True}

    def validate_child(self,item_name,item_value,whole_block,provisional_items={},globals={}):
        try:
            child_items = self[item_name][self.child_spec][:]  #copy
        except KeyError:
            return {"result":None}    #not relevant
        # special case for dictionaries  -> we check parents of children only
        if globals.has_key(item_name):  #dictionary so skip
            return {"result":None}
        if isinstance(child_items,StringType): # only one child
            child_items = [child_items]
        if isinstance(item_value,StringType): # single value
            parent_values = [item_value]
        else: parent_values = item_value[:]
        # expand child list with list of alternates
        for child_item in child_items[:]:
            child_items.extend(self.get_alternates(child_item))
        # now loop over the children
        for child_item in child_items:
            if self.optimize:
                if child_item in self.done_children:
                    return {"result":None}
                else: 
                    self.done_children.append(child_item)
                    print "Done children %s" % `self.done_children`
            if provisional_items.has_key(child_item):
                child_values = provisional_items[child_item][:]
            elif whole_block.has_key(child_item):
                child_values = whole_block[child_item][:]
            else:  continue 
            if isinstance(child_values,StringType):
                child_values = [child_values]
            #    print "Checking child %s against %s, values %s/%s" % (child_item,
            #                                          item_name,`child_values`,`parent_values`)
            missing = self.check_parent_child(parent_values,child_values)
            if len(missing)>0:
                return {"result":False,"bad_values":missing,"child":child_item}
        return {"result":True}       #could mean that no child items present
           
    #a generic checker: all child vals should appear in parent_vals
    def check_parent_child(self,parent_vals,child_vals):
        # shield ourselves from dots and question marks
        pv = parent_vals[:]
        pv.extend([".","?"])
        res =  filter(lambda a:a not in pv,child_vals)
        #print "Missing: %s" % res
        return res

    def validate_remove_parent_child(self,item_name,whole_block):
        try:
            child_items = self[item_name][self.child_spec]
        except KeyError:
            return {"result":None}
        if isinstance(child_items,StringType): # only one child
            child_items = [child_items]
        for child_item in child_items:
            if whole_block.has_key(child_item): 
                return {"result":False,"child":child_item}
        return {"result":True}
         
    def validate_dependents(self,item_name,item_value,whole_block,prov={},globals={}):
        try:
            dep_items = self[item_name][self.dep_spec][:]
        except KeyError:
            return {"result":None}    #not relevant
        if isinstance(dep_items,StringType):
            dep_items = [dep_items]
        actual_names = whole_block.keys()
        actual_names.extend(prov.keys())
        actual_names.extend(globals.keys())
        missing = filter(lambda a:a not in actual_names,dep_items)
        if len(missing) > 0:
            alternates = map(lambda a:[self.get_alternates(a),a],missing)
            # compact way to get a list of alternative items which are 
            # present
            have_check = map(lambda b:[filter_present(b[0],actual_names),
                                       b[1]],alternates) 
            have_check = filter(lambda a:len(a[0])==0,have_check)
            if len(have_check) > 0:
                have_check = map(lambda a:a[1],have_check)
                return {"result":False,"bad_items":have_check}
        return {"result":True}
        
    def validate_uniqueness(self,item_name,item_value,whole_block,provisional_items={},
                                                                  globals={}):
        category = self[item_name].get(self.cat_spec)
        if category == None:
            print "No category found for %s" % item_name
            return {"result":None}
        # print "Category %s for item %s" % (`category`,item_name)
        catentry = self.cat_map[category]
        # we make a copy in the following as we will be removing stuff later!
        unique_i = self[catentry].get("_category_key.name",[])[:]
        if isinstance(unique_i,StringType):
            unique_i = [unique_i]
        if item_name not in unique_i:       #no need to verify
            return {"result":None}
        if isinstance(item_value,StringType):  #not looped
            return {"result":None}
        # print "Checking %s -> %s -> %s ->Unique: " % (item_name,category,catentry) + `unique_i`
        # check that we can't optimize by not doing this check
        if self.optimize:
            if unique_i in self.done_keys:
                return {"result":None}
            else:
                self.done_keys.append(unique_i)
        val_list = []
        # get the matching data from any other data items
        unique_i.remove(item_name)
        other_data = []
        if len(unique_i) > 0:            # i.e. do have others to think about
           for other_name in unique_i:
           # we look for the value first in the provisional dict, then the main block
           # the logic being that anything in the provisional dict overrides the
           # main block
               if provisional_items.has_key(other_name):
                   other_data.append(provisional_items[other_name]) 
               elif whole_block.has_key(other_name):
                   other_data.append(whole_block[other_name])
               elif self[other_name].get(self.must_exist_spec)=="implicit":
                   other_data.append([item_name]*len(item_value))  #placeholder
               else:
                   return {"result":False,"bad_items":other_name}#missing data name
        # ok, so we go through all of our values
        # this works by comparing lists of strings to one other, and
        # so could be fooled if you think that '1.' and '1' are 
        # identical
        for i in range(len(item_value)):
            #print "Value no. %d" % i ,
            this_entry = item_value[i]
            for j in range(len(other_data)):
                this_entry = " ".join([this_entry,other_data[j][i]]) 
            #print "Looking for %s in %s: " % (`this_entry`,`val_list`)
            if this_entry in val_list: 
                return {"result":False,"bad_values":this_entry}
            val_list.append(this_entry)
        return {"result":True}


    def validate_mandatory_category(self,whole_block,globals={},fake_mand=False):
        if fake_mand:
            return {"result":True}
        mand_cats = filter(lambda a:self[a].get("_category.mandatory_code","no")=="yes",
                    self.keys())
        # map to actual ids
        catlist = self.cat_map.items()
        # print "Mandatory categories - %s" % `mand_cats`
        all_keys = whole_block.keys() #non-save block keys
        if globals:         #
            all_keys.extend(globals.abs_all_keys)
        for mand_cat in mand_cats:
            cat_id = filter(lambda a:a[1]==mand_cat,catlist)[0][0]
            no_of_items = len(filter(lambda a:self[a].get(self.cat_spec)==cat_id,
                                 all_keys))
            if no_of_items == 0:
                return {"result":False,"bad_items":cat_id}
        return {"result":True}

    def find_prob_cats(self,whole_block):
        mand_cats = filter(lambda a:self[a].get("_category.mandatory_code","no")=="yes",
                    self.keys())
        # map to actual ids
        catlist = self.cat_map.items()
        # find missing categories
        wbs = whole_block["saves"]
        abs_all_keys = whole_block.keys()
        abs_all_keys.extend(reduce(lambda a,b:a+(wbs[b].keys()),wbs.keys(),[]))
        prob_cats = []
        for mand_cat in mand_cats:
            cat_id = filter(lambda a:a[1]==mand_cat,catlist)[0][0]
            
            if len(filter(lambda a:self[a].get(self.cat_spec)==cat_id,abs_all_keys))==0:
                prob_cats.append(cat_id)
        if len(prob_cats) > 0:
            return (False,{'whole_block':[('validate_mandatory_category',{"result":False,"bad_items":problem_cats})]})
        else:
            return (True,{})


    def run_item_validation(self,item_name,item_value):
        return {item_name:map(lambda f:(f.__name__,f(item_name,item_value)),self.item_validation_funs)}

    def run_loop_validation(self,loop_names):
        return {loop_names[0]:map(lambda f:(f.__name__,f(loop_names)),self.loop_validation_funs)}

    def run_global_validation(self,item_name,item_value,data_block,provisional_items={},globals={}):
        results = map(lambda f:(f.__name__,f(item_name,item_value,data_block,provisional_items,globals)),self.global_validation_funs)
        return {item_name:results}

    def run_block_validation(self,whole_block,globals={},fake_mand=False):
        results = map(lambda f:(f.__name__,f(whole_block,globals,fake_mand)),self.block_validation_funs)
        # fix up the return values
        return {"whole_block":results}

    def optimize_on(self):
        self.optimize = True
        self.done_keys = []
        self.done_children = []
        self.done_parents = []

    def optimize_off(self):
        self.optimize = False
        self.done_keys = []
        self.done_children = []
        self.done_parents = []


class ValidCifBlock(CifBlock):
    def __init__(self,dic = None, diclist=[], mergemode = "replace",*args,**kwords):
        CifBlock.__init__(self,*args,**kwords)    
        if dic and diclist:
            print "Warning: diclist argument ignored when initialising ValidCifBlock"
        if isinstance(dic,CifDic):
            self.fulldic = dic
        else:
            raise TypeError( "ValidCifBlock passed non-CifDic type in dic argument")
        if len(diclist)==0 and not dic:
            raise ValidCifError( "At least one dictionary must be specified")
        if diclist and not dic:
            self.fulldic = merge_dic(diclist,mergemode)
        if not self.run_data_checks()[0]:
            raise ValidCifError( self.report())

    def run_data_checks(self,verbose=False):
        self.v_result = {}
        self.fulldic.optimize_on()
        for dataname in self.keys():
            update_value(self.v_result,self.fulldic.run_item_validation(dataname,self[dataname]))
            update_value(self.v_result,self.fulldic.run_global_validation(dataname,self[dataname],self))
        for loop in self.loops:
            update_value(self.v_result,self.fulldic.run_loop_validation(loop.keys()))
        # now run block-level checks
        update_value(self.v_result,self.fulldic.run_block_validation(self))
        # return false and list of baddies if anything didn't match
        self.fulldic.optimize_off()
        for test_key in self.v_result.keys():
            #print "%s: %s" % (test_key,`self.v_result[test_key]`)
            self.v_result[test_key] = filter(lambda a:a[1]["result"]==False,self.v_result[test_key])
            if len(self.v_result[test_key]) == 0: 
                del self.v_result[test_key]
        isvalid = len(self.v_result)==0
        #if not isvalid:
        #    print "Baddies:" + `self.v_result`
        return isvalid,self.v_result

    def single_item_check(self,item_name,item_value):
        #self.match_single_item(item_name)
        if not self.fulldic.has_key(item_name):
            result = {item_name:[]}
        else:
            result = self.fulldic.run_item_validation(item_name,item_value)
        baddies = filter(lambda a:a[1]["result"]==False, result[item_name])
        # if even one false one is found, this should trigger
        isvalid = (len(baddies) == 0)
        # if not isvalid: print "Failures for %s:" % item_name + `baddies`
        return isvalid,baddies

    def loop_item_check(self,loop_names):
        in_dic_names = filter(lambda a:self.fulldic.has_key(a),loop_names)
        if len(in_dic_names)==0:
            result = {loop_names[0]:[]}
        else:
            result = self.fulldic.run_loop_validation(in_dic_names)
        baddies = filter(lambda a:a[1]["result"]==False,result[in_dic_names[0]])
        # if even one false one is found, this should trigger
        isvalid = (len(baddies) == 0)
        # if not isvalid: print "Failures for %s:" % `loop_names` + `baddies`
        return isvalid,baddies

    def global_item_check(self,item_name,item_value,provisional_items={}):
        if not self.fulldic.has_key(item_name):
            result = {item_name:[]}
        else:
            result = self.fulldic.run_global_validation(item_name,
               item_value,self,provisional_items = provisional_items)
        baddies = filter(lambda a:a[1]["result"]==False,result[item_name])
        # if even one false one is found, this should trigger
        isvalid = (len(baddies) == 0)
        # if not isvalid: print "Failures for %s:" % item_name + `baddies`
        return isvalid,baddies

    def remove_global_item_check(self,item_name):
        if not self.fulldic.has_key(item_name):
            result = {item_name:[]}
        else:
            result = self.fulldic.run_remove_global_validation(item_name,self,False)
        baddies = filter(lambda a:a[1]["result"]==False,result[item_name])
        # if even one false one is found, this should trigger
        isvalid = (len(baddies) == 0)
        # if not isvalid: print "Failures for %s:" % item_name + `baddies`
        return isvalid,baddies

    def AddToLoop(self,dataname,loopdata):
        # single item checks
        paired_data = loopdata.items()
        for name,value in paired_data:
            valid,problems = self.single_item_check(name,value) 
            self.report_if_invalid(valid,problems)
        # loop item checks; merge with current loop
        found = 0
        for aloop in self.block["loops"]:
            if aloop.has_key(dataname):
                loopnames = aloop.keys()
                for new_name in loopdata.keys():
                    if new_name not in loopnames: loopnames.append(new_name)
                valid,problems = self.looped_item_check(loopnames)
                self.report_if_invalid(valid,problems)
        prov_dict = loopdata.copy()
        for name,value in paired_data: 
            del prov_dict[name]   # remove temporarily
            valid,problems = self.global_item_check(name,value,prov_dict)
            prov_dict[name] = value  # add back in
            self.report_if_invalid(valid,problems)
        CifBlock.AddToLoop(self,dataname,loopdata)
 
    def AddCifItem(self,data):
        if isinstance(data[0],StringType):   # single item
            valid,problems = self.single_item_check(data[0],data[1])
            self.report_if_invalid(valid,problems,data[0])
            valid,problems = self.global_item_check(data[0],data[1])
            self.report_if_invalid(valid,problems,data[0])
        elif isinstance(data[0],TupleType) or isinstance(data[0],ListType):
            paired_data = map(None,data[0],data[1])
            for name,value in paired_data:
                valid,problems = self.single_item_check(name,value) 
                self.report_if_invalid(valid,problems,name)
            valid,problems = self.loop_item_check(data[0])
            self.report_if_invalid(valid,problems,data[0])
            prov_dict = {}            # for storing temporary items
            for name,value in paired_data: prov_dict[name]=value
            for name,value in paired_data: 
                del prov_dict[name]   # remove temporarily
                valid,problems = self.global_item_check(name,value,prov_dict)
                prov_dict[name] = value  # add back in
                self.report_if_invalid(valid,problems,name)
        CifBlock.AddCifItem(self,data)

    # utility function
    def report_if_invalid(self,valid,bad_list,data_name):
        if not valid:
            error_string = reduce(lambda a,b: a + "," + b[0], bad_list, "") 
            error_string = `data_name` + " fails following validity checks: "  + error_string
            raise ValidCifError( error_string)

    def __delitem__(self,key):
        # we don't need to run single item checks; we do need to run loop and
        # global checks.
        if self.has_key(key):
            try: 
                loop_items = self.GetLoop(key)
            except TypeError:
                loop_items = []
            if loop_items:             #need to check loop conformance
                loop_names = map(lambda a:a[0],loop_items)
                loop_names = filter(lambda a:a != key,loop_names)
                valid,problems = self.loop_item_check(loop_names)
                self.report_if_invalid(valid,problems)
            valid,problems = self.remove_global_item_check(key)
            self.report_if_invalid(valid,problems)
        self.RemoveCifItem(key)


    def report(self):
       import cStringIO
       outstr = cStringIO.StringIO()
       outstr.write( "Validation results\n")
       outstr.write( "------------------\n")
       print "%d invalid items found\n" % len(self.v_result)
       for item_name,val_func_list in self.v_result.items():
           outstr.write("%s fails following tests:\n" % item_name)
           for val_func in val_func_list:
               outstr.write("\t%s\n")
       return outstr.getvalue()


class ValidCifFile(CifFile):
    def __init__(self,dic=None,diclist=[],mergemode="replace",*args,**kwargs):
        if not diclist and not dic and not hasattr(self,'bigdic'):
            raise ValidCifError( "At least one dictionary is required to create a ValidCifFile object")
        if not dic and diclist:     #merge here for speed
            self.bigdic = merge_dic(diclist,mergemode)
        elif dic and not diclist:
            self.bigdic = dic
        CifFile.__init__(self,*args,**kwargs)
        #for blockname in self.keys():
    #       self.dictionary[blockname]=ValidCifBlock(data=self.dictionary[blockname],dic=self.bigdic)

    def NewBlock(self,blockname,blockcontents,**kwargs):
        CifFile.NewBlock(self,blockname,blockcontents,**kwargs)
        # dictionary[blockname] is now a CifBlock object.  We
        # turn it into a ValidCifBlock object
        self.dictionary[blockname] = ValidCifBlock(dic=self.bigdic,
                                         data=self.dictionary[blockname])


class ValidationResult:
    """Represents validation result"""
    def __init__(self,results):
        """results is return value of validate function"""
        self.valid_result, self.no_matches = results

    def report(self,use_html):
        """Return string with human-readable description of validation result"""
        return validate_report((self.valid_result, self.no_matches),use_html)

    def is_valid(self,block_name=None):
        """Return True for valid CIF file, otherwise False"""
        if block_name is not None:
            block_names = [block_name]
        else:
            block_names = self.valid_result.iterkeys()
        for block_name in block_names:
            if not self.valid_result[block_name] == (True,{}):
                valid = False
                break
            else:
                valid = True
        return valid
    
    def has_no_match_items(self,block_name=None):
        """Return true if some items are not found in dictionary"""
        if block_name is not None:
            block_names = [block_name]
        else:
            block_names = self.no_matches.iter_keys() 
        for block_name in block_names:
            if self.no_matches[block_name]:
                has_no_match_items = True
                break
            else:
                has_no_match_items = False
        return has_no_match_items
    

        
def validate(ciffile,dic = "", diclist=[],mergemode="replace",isdic=False,fake_mand=True):
    check_file = CifFile(ciffile)
    if not dic:
        fulldic = merge_dic(diclist,mergemode)
    else:
        fulldic = dic
    no_matches = {}
    valid_result = {}
    if isdic:          #assume one block only
        blockname = check_file.keys()[0]
        check_bc = check_file[blockname]["saves"]
        check_globals = check_file[blockname] 
        # collect a list of parents for speed
        poss_parents = fulldic.get_all("_item_linked.parent_name")
        for parent in poss_parents:
            curr_parent = listify(check_globals.get(parent,[]))
            new_vals = check_bc.get_all(parent)
            new_vals.extend(curr_parent)
            if len(new_vals)>0:
                check_globals[parent] = new_vals
                # print "Added %s (len %d)" % (parent,len(check_globals[parent]))
        # next dictionary problem: the main DDL2 dictionary has what
        # I would characterise as a mandatory_category problem, but
        # in order to gloss over it, we allow a different 
        # interpretation, which requires only a single check for one
        # block.
        if fake_mand:
            valid_result[blockname] = fulldic.find_prob_cats(check_globals)
            no_matches[blockname] = filter(lambda a:not fulldic.has_key(a),check_globals.keys())
    else:
        check_bc = check_file
        check_globals = CifBlock()   #empty
    for block in check_bc.keys(): 
        #print "Validating block %s" % block 
        no_matches[block] = filter(lambda a:not fulldic.has_key(a),check_bc[block].keys())
        # remove non-matching items
        # print "Not matched: " + `no_matches[block]`
        for nogood in no_matches[block]:
             del check_bc[block][nogood]
        valid_result[block] = run_data_checks(check_bc[block],fulldic,globals=check_globals,fake_mand=fake_mand)
    return valid_result,no_matches

def validate_report(val_result,use_html=False):
    import cStringIO
    valid_result,no_matches = val_result
    outstr = cStringIO.StringIO()
    if use_html:
        outstr.write("<h2>Validation results</h2>")
    else:
        outstr.write( "Validation results\n")
        outstr.write( "------------------\n")
    if len(valid_result) > 10:  
        suppress_valid = True         #don't clutter with valid messages
        if use_html:
           outstr.write("<p>For brevity, valid blocks are not reported in the output.</p>")
    else:
        suppress_valid = False
    for block in valid_result.keys():
        block_result = valid_result[block]
        if block_result[0]:
            out_line = "Block '%s' is VALID" % block
        else:
            out_line = "Block '%s' is INVALID" % block
        if use_html:
            if (block_result[0] and (not suppress_valid or len(no_matches[block])>0)) or not block_result[0]:
                outstr.write( "<h3>%s</h3><p>" % out_line)
        else:
                outstr.write( "\n %s\n" % out_line)
        if len(no_matches[block])!= 0:
            if use_html:
                outstr.write( "<p>The following items were not found in the dictionary")
                outstr.write(" (note that this does not invalidate the data block):</p>")
                outstr.write("<p><table>\n")
                map(lambda it:outstr.write("<tr><td>%s</td></tr>" % it),no_matches[block])
                outstr.write("</table>\n")
            else:
                outstr.write( "\n The following items were not found in the dictionary:\n")
                outstr.write("Note that this does not invalidate the data block\n")
                map(lambda it:outstr.write("%s\n" % it),no_matches[block])
        # now organise our results by type of error, not data item...
        error_type_dic = {}
        for error_item, error_list in block_result[1].items():
            for func_name,bad_result in error_list:
                bad_result.update({"item_name":error_item})
                try:
                    error_type_dic[func_name].append(bad_result)
                except KeyError:
                    error_type_dic[func_name] = [bad_result]
        # make a table of test name, test message
        info_table = {\
        'validate_item_type':\
            "The following data items had badly formed values",
        'validate_item_esd':\
            "The following data items should not have esds appended",
        'validate_enum_range':\
            "The following data items have values outside permitted range",
        'validate_item_enum':\
            "The following data items have values outside permitted set",
        'validate_looping':\
            "The following data items violate looping constraints",
        'validate_loop_membership':\
            "The following looped data names are of different categories to the first looped data name",
        'validate_loop_key':\
            "A required dataname for this category is missing from the loop\n containing the dataname",
        'validate_loop_references':\
            "A dataname required by the item is missing from the loop",
        'validate_parent':\
            "A parent dataname is missing or contains different values",
        'validate_child':\
            "A child dataname contains different values to the parent",
        'validate_uniqueness':\
            "One or more data items do not take unique values",
        'validate_dependents':\
            "A dataname required by the item is missing from the data block",
        'validate_exclusion': \
            "Both dataname and exclusive alternates or aliases are present in data block",
        'validate_mandatory_category':\
            "A required category is missing from this block"}

        for test_name,test_results in error_type_dic.items():
           if use_html:
               outstr.write(html_error_report(test_name,info_table[test_name],test_results)) 
           else:
               outstr.write(error_report(test_name,info_table[test_name],test_results)) 
               outstr.write("\n\n")
    return outstr.getvalue()
         
# A function to lay out a single error report.  We are passed
# the name of the error (one of our validation functions), the
# explanation to print out, and a dictionary with the error 
# information.  We print no more than 50 characters of the item

def error_report(error_name,error_explanation,error_dics):
   retstring = "\n\n " + error_explanation + ":\n\n"
   headstring = "%-32s" % "Item name"
   bodystring = ""
   if error_dics[0].has_key("bad_values"):
      headstring += "%-20s" % "Bad value(s)"
   if error_dics[0].has_key("bad_items"):
      headstring += "%-20s" % "Bad dataname(s)"
   if error_dics[0].has_key("child"):
      headstring += "%-20s" % "Child"
   if error_dics[0].has_key("parent"):
      headstring += "%-20s" % "Parent" 
   headstring +="\n"
   for error in error_dics:
      bodystring += "\n%-32s" % error["item_name"]
      if error.has_key("bad_values"):
          out_vals = map(lambda a:a[:50],error["bad_values"])
          bodystring += "%-20s" % out_vals 
      if error.has_key("bad_items"):
          bodystring += "%-20s" % error["bad_items"]
      if error.has_key("child"):
          bodystring += "%-20s" % error["child"]
      if error.has_key("parent"):
          bodystring += "%-20s" % error["parent"]
   return retstring + headstring + bodystring 

#  This lays out an HTML error report

def html_error_report(error_name,error_explanation,error_dics,annotate=[]):
   retstring = "<h4>" + error_explanation + ":</h4>"
   retstring = retstring + "<table cellpadding=5><tr>"
   headstring = "<th>Item name</th>"
   bodystring = ""
   if error_dics[0].has_key("bad_values"):
      headstring += "<th>Bad value(s)</th>"
   if error_dics[0].has_key("bad_items"):
      headstring += "<th>Bad dataname(s)</th>"
   if error_dics[0].has_key("child"):
      headstring += "<th>Child</th>"
   if error_dics[0].has_key("parent"):
      headstring += "<th>Parent</th>" 
   headstring +="</tr>\n"
   for error in error_dics:
      bodystring += "<tr><td><tt>%s</tt></td>" % error["item_name"]
      if error.has_key("bad_values"):
          bodystring += "<td>%s</td>" % error["bad_values"]
      if error.has_key("bad_items"):
          bodystring += "<td><tt>%s</tt></td>" % error["bad_items"]
      if error.has_key("child"):
          bodystring += "<td><tt>%s</tt></td>" % error["child"]
      if error.has_key("parent"):
          bodystring += "<td><tt>%s</tt></td>" % error["parent"]
      bodystring += "</tr>\n"
   return retstring + headstring + bodystring + "</table>\n"

def run_data_checks(check_block,fulldic,globals={},fake_mand=False):
    v_result = {}
    for key in check_block.keys():
        update_value(v_result, fulldic.run_item_validation(key,check_block[key]))
        update_value(v_result, fulldic.run_global_validation(key,check_block[key],check_block,globals=globals))
    for loop in check_block.loops:
        update_value(v_result, fulldic.run_loop_validation(loop.keys()))
    update_value(v_result,fulldic.run_block_validation(check_block,globals=globals,fake_mand=fake_mand))
    # return false and list of baddies if anything didn't match
    for test_key in v_result.keys():
        v_result[test_key] = filter(lambda a:a[1]["result"]==False,v_result[test_key])
        if len(v_result[test_key]) == 0: 
            del v_result[test_key]
    # if even one false one is found, this should trigger
    # print "Baddies:" + `v_result`
    isvalid = len(v_result)==0
    return isvalid,v_result
    

def get_number_with_esd(numstring):
    import string
    numb_re = '((-?(([0-9]*[.]([0-9]+))|([0-9]+)[.]?))([(][0-9]+[)])?([eEdD][+-]?[0-9]+)?)|(\?)|(\.)' 
    our_match = re.match(numb_re,numstring)
    if our_match:
        a,base_num,b,c,dad,dbd,esd,exp,q,dot = our_match.groups()
    #    print "Debug: %s -> %s" % (numstring, `our_match.groups()`)
    else:
        return None,None
    if dot or q: return None,None     #a dot or question mark
    if exp:          #has exponent 
       exp = string.replace(exp,"d","e")     # mop up old fashioned numbers
       exp = string.replace(exp,"D","e")
       base_num = base_num + exp
    #print "Debug: have %s for base_num from %s" % (base_num,numstring)
    base_num = float(base_num)
    # work out esd, if present.
    if esd:
        esd = float(esd[1:-1])    # no brackets
        if dad:                   # decimal point + digits
            esd = esd * (10 ** (-1* len(dad)))
        if exp:
            esd = esd * (10 ** (float(exp[1:])))
    return base_num,esd

def float_with_esd(inval):
    if isinstance(inval,StringType):
        j = inval.find("(")
        if j>=0:  return float(inval[:j])
    return float(inval)
        
    
                
# A utility function to append to item values rather than replace them
def update_value(base_dict,new_items):
    for new_key in new_items.keys():
        if base_dict.has_key(new_key):
            base_dict[new_key].extend(new_items[new_key])
        else:
            base_dict[new_key] = new_items[new_key]

#Transpose the list of lists passed to us
def transpose(base_list):
    new_lofl = []
    full_length = len(base_list)
    opt_range = range(full_length)
    for i in range(len(base_list[0])):
       new_packet = [] 
       for j in opt_range:
          new_packet.append(base_list[j][i])
       new_lofl.append(new_packet)
    return new_lofl

# listify strings - used surprisingly often
def listify(item):
    if isinstance(item,StringType): return [item]
    else: return item

# given a list of search items, return a list of items 
# actually contained in the given data block
def filter_present(namelist,datablocknames):
    return filter(lambda a:a in datablocknames,namelist)

# merge ddl dictionaries.  We should be passed filenames or CifFile
# objects
def merge_dic(diclist,mergemode="replace",ddlspec=None):
    dic_as_cif_list = []
    for dic in diclist:
        if not isinstance(dic,CifFile) and \
           not isinstance(dic,StringType):
               raise TypeError, "Require list of CifFile names/objects for dictionary merging"
        if not isinstance(dic,CifFile): dic_as_cif_list.append(CifFile(dic))
        else: dic_as_cif_list.append(dic)
    # we now merge left to right
    basedic = dic_as_cif_list[0]
    if basedic.has_key("on_this_dictionary"):   #DDL1 style only
        for dic in dic_as_cif_list[1:]:
           basedic.merge(dic,mode=mergemode,match_att=["_name"])
    elif len(basedic.keys()) == 1:                     #One block: DDL2 style
        old_block = basedic[basedic.keys()[0]]
        for dic in dic_as_cif_list[1:]:
           new_block = dic[dic.keys()[0]]
           basedic.merge(dic,mode=mergemode,
                         single_block=[basedic.keys()[0],dic.keys()[0]],
                         match_att=["_item.name"],match_function=find_parent)
    return CifDic(basedic)

def find_parent(ddl2_def):
    if not ddl2_def.has_key("_item.name"):
       return None 
    if isinstance(ddl2_def["_item.name"],StringType):
        return ddl2_def["_item.name"]
    if not ddl2_def.has_key("_item_linked.child_name"):
        raise CifError("Asked to find parent in block with no child_names")
    if not ddl2_def.has_key("_item_linked.parent_name"):
        raise CifError("Asked to find parent in block with no parent_names")
    result = filter(lambda a:a not in ddl2_def["_item_linked.child_name"],ddl2_def["_item.name"]) 
    if len(result)>1 or len(result)==0:
        raise CifError("Unable to find single unique parent data item")
    return result[0]


def ReadCif(filename,strict=1,maxlength=2048,scantype="standard",grammar="1.1"):
    proto_cif = StarFile.ReadStar(filename,maxlength,scantype=scantype,grammar=grammar)
    # convert to CifFile
    proto_cif = CifFile(proto_cif)
    # check for nested loops
    for bname,bvalue in proto_cif.items():
        nests = filter(lambda a:len(a.loops)>0,bvalue.loops)
        if len(nests) > 0:
            raise CifError( "Block %s contains nested loops")
        # check for save frame references (not yet implemented in PySTARRW)
        # check for global blocks (not yet implemented in PySTARRW)
    return proto_cif


