#!/usr/bin/env python
# -*- encoding: iso-8859-1 -*-
"""
Generate workflows view

Author: Dan Gunter <dkgunter@lbl.gov>
Created: 14 October 2011
"""
__rcsid__ = "$Id$"

import datetime
import pymongo

# States
CR, BL, RN, ST, KL, AB, CM = (
    "created", "building", "running",
    "stopped", "killed", "abyss", "completed")

# State transitions (links from child to parent)
parents = { BL : CR,
            RN : BL,
            KL : RN,
            AB : RN,
            CM : RN,
            ST : RN }

# Maximum time since a task was updated
# before we consider it dead.
MAX_TASK_IDLE = 60*60*8

class Workflow:
    """Graph of states.

    The tree is represented internally as a nested sequence.
    For example, if the tree is:
      a
     /\ \
    b  c f
       /\
      d  e

   The list will be: (a (b c (d e) f))
    """
    def __init__(self, root="start"):
        """Create with a root node.
        """
        self._tree = [root]
        self._nodes = 1
        
    def __len__(self):
        return self._nodes
    
    def add(self, value):
        """Add a sub-tree, pointed to by root.
        """
        self._nodes += len(value)
        self._tree.append(value)

    def add_subwf(self, state, suffix=""):
        """Add a subwf, off the root, for the given current state.
        """
        subwf = nest(path(state, suffix=suffix))
        self.add(subwf)

    def chain_subwf(self, state, item, suffix=""):
        """Add a subwf, off the i-th item, for the given
        current state.
        """
        subwf = nest(path(state, suffix=suffix))
        self._nodes += len(subwf)
        item.append(subwf)

    def __add__(self, other):
        self.add(other._tree)
        return(self)
        
    def as_dot(self):
        """Represent for graphviz.

        Returns:
           (str) Graphviz commands
        """
        def _format(parent, child, links):
             if parent:
                 links.append("{0} -> {1};".format(parent, child))
        links = [ ]
        traverse(None, self._tree, _format, (links,))
        result = ("digraph workflow {\n"
                  "rankdir=LR;\n" +
                  '\n'.join(links) +
                  "}\n")
        return result
    
def build_xtl_workflow(xtl_id,  engines, tasks):
    """Build the crystal workflow usin the engine and task collections.
    """
    # get the engine IDs associated with this crystal
    ee = list(engines.find({'crystal_id':xtl_id},
                            fields=('engine_id', 'state', 'failed'),
                            sort=[('created_at', 1)]))
    if not ee:
        return None # no workflow to build
    # init workflow obj
    xtl_idstr = "_{0:d}".format(xtl_id)
    wf = Workflow(root = "Crystal" + xtl_idstr)
    # for each engine-id, retrieve associated info from task table
    # and append as a sub-workflow to the workflow
    prev, cur = None, None
    has_multi = len(ee) > 1
    for i, item in enumerate(ee):
        eid = item['engine_id']
        task = tasks.find_one({'task_id':eid})
        if not task:
            raise ValueError("No task for engine_id {0:d}".format(eid))
        state = determine_state(item['state'], task)
        if has_multi:
            sfx = "_{0:d}".format(i+1)
        else:
            sfx = ""
        sfx += xtl_idstr
        cur = subwf(state, suffix=sfx)
        if prev:
            prev.append(cur) # add to root of wf
        else:
            wf.add(cur) # chain off of previous subwf
        prev = cur
    # that's it, return the workflow
    return wf

def determine_state(engine_state, task):
    if engine_state == CR:
        return CR
    if engine_state == RN:
        if not task['was_successful']:
            return KL
        elif updated_sec(task) > MAX_TASK_IDLE:
            return AB
        else:
            return RN
    if engine_state == CM:
        if task['was_stopped']:
            return ST
        else:
            return CM
    raise ValueError("Unexpected engine state: " + engine_state)

def updated_sec(task):
    now = datetime.datetime.now()
    return (now - task['updated_at']).seconds

def subwf(state, **kw):
    return nest(path(state, **kw))

def path(state, suffix=''):
    """Make list of all states up to and including input state.

    Args:
      state - (string) name of a state.
    Returns:
      list of all states up to and including input state
      
    For example:
    >>> path('foo')
    ['foo']
    >>> path('killed')
    ['created', 'building', 'running', 'killed']
    """
    p = parents.get(state, None)
    state_name = (state, state + suffix)[bool(suffix)]
    if p is None:
        return [state_name]
    else:
        return path(p, suffix=suffix) + [state_name]
        
def nest(states):
    """Nest a list of strings, so that it fits the tree structure
    expected in Workflow.

    Args:
       states - list of string
    Returns:
       nested list with same strings
       
    For example,
    >>> nest(["start", CR, BL, RN])
    ['start', ['created', ['building', ['running']]]]
    """
    if len(states) == 1:
        return states
    else:
        return [states[0], nest(states[1:])]

def traverse(prev, cur, func, args):
    """LR depth-first traversal of a list (iterable) representing a tree.
   """
    prev_terminal = None
    for node in cur:
        if isinstance(node, str):
            func(prev, node, *args)
            prev_terminal = node
        else:
            if not prev_terminal:
                raise ValueError("bad tree structure")
            traverse(prev_terminal, node, func, args)

# Run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()

# Run another test
if __name__ == '__main__':
    import sys
    conf = open("wf.conf")
    user, passwd = conf.readline().split() 
    db1 = pymongo.Connection('materialsdb.lbl.gov').test_automation_v_0_1
    db1.authenticate(user, passwd)
    engines_coll = db1.engines
    db2 = pymongo.Connection('materials.nersc.gov').mg_core_v_0_5_prod
    user, passwd = conf.readline().split() 
    db2.authenticate(user, passwd)
    tasks_coll = db2.tasks
    master_wf = Workflow('Crystals')
    xtl = 1
    while len(master_wf) < 10:        
        try:
            wf = build_xtl_workflow(xtl, engines_coll, tasks_coll)
        except ValueError, err:
            print("# ERROR: {0}".format(err))
            wf = None
        if wf:
            master_wf = master_wf + wf
        else:
            print("# Did not find crystal '{0:d}'".format(xtl))
        xtl += 1
    print(master_wf.as_dot())
#
# 
#

# Local Variables:
# mode:python
# End:
