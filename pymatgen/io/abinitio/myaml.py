import yaml

_HAS_CYAML = False

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    import warnings
    msg = ("Cannot import LibYAML bindings. "
           "If you want to use LibYAML bindings, which are much faster than the pure Python version,\n" +
           "you need to download and install LibYAML. See http://pyyaml.org/wiki/PyYAMLDocumentation\n")
    #warnings.warn(msg)
    _HAS_CYAML = False
    from yaml import Loader, Dumper


def has_cyaml():
    """True if we have the C implementation"""
    return _HAS_CYAML


def load(stream, Loader=Loader):
    return yaml.load(stream, Loader=Loader)


def dump(data, Dumper=Dumper):
    return yaml.dump(data, Dumper=Dumper)
