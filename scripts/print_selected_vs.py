# blender script that prints selectes vertices' ids

import bpy
import bmesh

obj = bpy.context.edit_object
me = obj.data

bm = bmesh.from_edit_mesh(me)

l = []

for v in bm.verts:
    if v.select:
        l.append(v.index)

txt = ""
for i in l:
    txt += str(i) + ", "
    
print(txt)