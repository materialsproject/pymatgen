from pymatgen.ext.matproj import MPRester, TaskType

material_ids = ["mp-32800", "mp-23494"]
task_types = [TaskType.GGA_OPT, TaskType.GGA_UNIFORM]
file_patterns = ["vasprun*", "OUTCAR*"]
with MPRester("BsTwIdAypx0Far6t") as mpr:
    meta, urls = mpr.get_download_info(
        material_ids, task_types=task_types, file_patterns=file_patterns
    )

print(meta)
print(urls)
