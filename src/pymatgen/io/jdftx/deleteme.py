from __future__ import annotations


class SuperClass:
    def __init__(self):
        self.a = 1


class SubClass(SuperClass):
    def __init__(self):
        super().__init__()
        self.b = 2


class SubSubClass(SubClass):
    def __init__(self):
        super().__init__()
        self.c = 3


sc = SuperClass()
subc = SubClass()
subsubc = SubSubClass()

print(isinstance(sc, SuperClass))
print(isinstance(subc, SuperClass))
print(isinstance(subsubc, SuperClass))
