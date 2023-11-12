
reinitialize
bg white
space cmyk

load 7ar6.pdb

as cartoon
set cartoon_rect_length, 1
set cartoon_oval_length, 1
set stick_radius, 0.25
set sphere_scale, 0.5
set transparency, 0.0
set surface_color, grey

set ambient, 0.5

symexp sym, 7ar6, (7ar6), 0.1
set_name sym01000000, sym
remove not alt ''+A
dss

hide ////HOH
hide ////CL
hide ////DMS
#hide cartoon
show surface
remove hydrogens

#color yellow, 7ar6
set surface_color, orange, 7ar6

#color blue, sym
set surface_color, skyblue, sym

color atomic, (not elem C)

deselect
set ray_trace_mode, 1

set_view (\
    -0.367271572,    0.102026783,   -0.924502552,\
     0.033059150,   -0.991902709,   -0.122598901,\
    -0.929527640,   -0.075589187,    0.360924602,\
     0.000000000,    0.000000000, -267.076599121,\
     0.000003815,    0.407508850,    0.000011444,\
   234.994918823,  299.158386230,  -20.000000000 )

ray 1024,1024
save dimer-cartoon.png

hide everything, 7ar6

ray 1024,1024
save monomerA-cartoon.png

hide everything, sym
show surface, 7ar6

ray 1024,1024
save monomerB-cartoon.png