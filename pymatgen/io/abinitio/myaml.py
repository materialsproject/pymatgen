import yaml

HAS_CYAML = True
HAS_CYAML = False
from yaml import Loader, Dumper

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    import warnings
    warnings.warn("Cannot import LibYAML bindings. If you want to use LibYAML bindings, which are much faster than the pure Python version,\n" +
                  "you need to download and install LibYAML. See http://pyyaml.org/wiki/PyYAMLDocumentation\n")
    HAS_CYAML = False
    from yaml import Loader, Dumper

def has_cyaml():
    return HAS_CYAML


def load(stream, Loader=Loader):
    return yaml.load(stream, Loader=Loader)


def dump(data, Dumper=Dumper):
    return yaml.dump(data, Dumper=Dumper)
