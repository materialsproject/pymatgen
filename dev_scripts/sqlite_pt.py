import json
import ruamel.yaml as yaml
import sqlite3
from monty.serialization import loadfn


def main():
    with open("periodic_table.yaml") as f:
        contents = yaml.load(f)

    keys = list(contents["Li"].keys())
    definitions = [("symbol", "text")]
    for k in keys:
        if k in ["Atomic mass", "Atomic radius", "X"]:
            dtype = "real"
        elif k in ["Atomic no", "Mendeleev no"]:
            dtype = "int"
        else:
            dtype = "text"
        definitions.append((k, dtype))
    deftext = ", ".join(["%s %s" % (k.replace(" ", "_"), v) for k, v in definitions])

    conn = sqlite3.connect('periodic_table.db')
    c = conn.cursor()
    c.execute('CREATE TABLE elements (%s)' % deftext)

    for el, d in contents.items():
        data = ["'%s'" % el]
        for k, dtype in definitions:
            if k != "symbol":
                if dtype == "text":
                    if k in d:
                        data.append("'%s'" % d[k])
                    else:
                        data.append("NULL")
                else:
                    val = d.get(k, 'NULL')
                    val = "NULL" if val == "no data" else val
                    data.append("%s" % val)

        datatext = ", ".join(data)
        print(datatext)
        c.execute("INSERT INTO elements VALUES (%s)" % datatext)

    # Save (commit) the changes
    conn.commit()

    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()

def load_sql():
    conn = sqlite3.connect('periodic_table.db')
    c = conn.cursor()
    c.execute("SELECT * FROM elements")
    return c.fetchall()


def load_json():
    with open("periodic_table.json") as f:
        return json.load(f)

def load_yaml():
    with open("periodic_table.yaml") as f:
        return yaml.load(f)

def cload_yaml():
    with open("periodic_table.yaml") as f:
        return yaml.load(f, Loader=yaml.CSafeLoader)

def monty_load_yaml():
    return loadfn("periodic_table.yaml")

if __name__ == "__main__":
    import timeit
    n = 20
    print(timeit.timeit("load_sql()", setup="from __main__ import load_sql", number=n))
    print(timeit.timeit("load_json()", setup="from __main__ import load_json", number=n))
    print(timeit.timeit("cload_yaml()", setup="from __main__ import cload_yaml", number=n))
    print(timeit.timeit("monty_load_yaml()", setup="from __main__ import monty_load_yaml", number=n))
