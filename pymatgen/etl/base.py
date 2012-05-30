"""
Interfaces and utilities to perform
Extract-Transform-Load operations 
between MongoDB collections.
"""
__author__ = "Dan Gunter"
__copyright__ = "Copyright 2012, The Materials Project"
__maintainer__ = "Dan Gunter"
__email__ = "dkgunter@lbl.gov"
__date__ = "29 May 2012"
__rcsid__ = "$Id$"

import importlib
from StringIO import StringIO
import yaml
import warnings
try:
    import pymongo
except ImportError:
    pymongo = None
    warnings.warn("Failed to import 'pymongo'. "
    "Database operations will fail with NoneType error.")

class ETLError(Exception):
    def __init__(self, src=None, tgt=None, msg=None, base_exc="None"):
        if src is None:
            s = "Failed: {m}. Base exception: {e}".format(m=msg, e=base_exc)
        else:
            s = "Failed for source {s} -> target {t}. " \
                "Base exception: {e}".format(s=src, t=tgt, e=base_exc)
        Exception.__init__(self, s)

class ETLBase(object):
    """Base class for extract-transform-load.
    """
    def __init__(self, src=None, tgt=None):
        self.src, self.tgt = src, tgt
        
    def extract_transform_load(self):
        """Subclasses must override this to actually
        perform the operation.
        
        Returns: None
        Raises: Any Exception
        """
        return None
        
class ETLRunner:
    """Working from a YAML configuration file,
    perform arbitrary extract-transform-load (ETL) operations
    from one or more source collections to a target collection.
    """
    # Configuration layout with defaults.
    # If default is 'None' then the value must be provided in the
    # configuration file.
    CONF = {
        "sources" : {
            "collection" : None,
            "module" : None,
            "class" : "ETL",
            "param" : { }
        },
        "target": {
            "host" : "localhost",
            "port" : 27017,
            "user" : "",
            "password" : "",
            "db" : None,
            "collection" : None
        },
    }

    def __init__(self, conf):
        """Configure with YAML from file or string, or
        a pre-configured dictionary.
        
        Raises: ValueError, yaml.YAMLError
        """
        if hasattr(conf, "read"):
            self._conf = yaml.load(conf)
        elif hasattr(conf, "get"):
            self._conf = conf
        else:
            self._conf = yaml.load(StringIO(conf))
        self.target, self.sources = None, [ ]
        for section in self.CONF:
            if not section in self._conf:
                raise ValueError("Missing section: {}".format(section))
            if section == "target":
                contents = self._conf[section]
                self.target = self._get_values(section, contents)
            elif section == "sources":
                for contents in self._conf[section]:
                    values = self._get_values(section, contents)
                    self.sources.append(values)
            else:
                raise ValueError("Unknown section: {}".format(section))

    # Helpers for __init__
    
    def _get_values(self, section, data):
        v = { }
        for key, default in self.CONF[section].iteritems():
            uval = data.get(key, default)
            if uval is None:
                raise ValueError("Missing key: {sec}.{key}".format(
                sec=section, key=key))
            v[key] = uval
        return v
                  
    def __len__(self):
        """Number of source collections"""
        return len(self.sources)
        
    # Run
    
    def run(self):
        """Run all the ETL operations.
        
        Raises: ETLError
        Returns: Number run
        """
        dbconn = self._connect()
        target = dbconn[self.target["collection"]]
        n = 0
        for s in self.sources:
            source = dbconn[s["collection"]]
            try:
                etl_mod = self._load_module(s["module"])
            except ImportError, err:
                raise ETLError(src=source, tgt=target, base_exc=err)                
            etl_cls = getattr(etl_mod, s["class"])
            etl_param = s["param"]
            etl = etl_cls(src=source, tgt=target, **etl_param)
            try:
                etl.extract_transform_load()
            except ETLError:
                raise
            except Exception, err:
                raise ETLError(src=source, tgt=target, base_exc=err)
            n += 1
        return n
    # helpers for run()

    def _load_module(self, mod_name):
        """Dynamically load a Python module.
        
        Returns: module object
        Raises: ImportError
        """
        if mod_name[0] == ".":
            # perform relative imports from pymatgen package
            mod = importlib.import_module(mod_name, "pymatgen")
        else:
            mod = importlib.import_module(mod_name)
        return mod

    def _connect(self):
        """Connect and authorize to MongoDB
        
        Returns: DB connection obj
        Raises: ETLError
        """
        cfg = self.target
        host, port, db = cfg["host"], cfg["port"], cfg["db"]
        if cfg["user"]:
            user, passwd = cfg["user"], cfg["password"]
            uri = ("mongodb://{u}:{w}@{h}:{p}/{d}"
                   .format(u=user, w=passwd, h=host, p=port, d=db))
        else:
            uri = "mongodb://{h}:{p}".format(h=host, p=port)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                conn = pymongo.Connection(uri)
        except pymongo.errors.ConnectionFailure, err:
            raise ETLError(msg="Connect to {h}:{p}/{d}"
                            .format(h=host, p=port, d=db, base_exc=err))
        dbconn = conn[db]
        return dbconn