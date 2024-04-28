
reinitialize
bg white
space cmyk

load 7ar6.pdb

as cartoon
set cartoon_rect_length, 1
set cartoon_oval_length, 1
set stick_radius, 0.25
set sphere_scale, 0.5
set transparency, 0.7
set surface_color, grey

symexp sym, 7ar6, (7ar6), 0.1
set_name sym01000000, sym
remove not alt ''+A
dss

hide ////HOH
hide ////CL
hide ////DMS
remove hydrogens

color grey70, 7ar6
show surface, 7ar6

select active_site, /7ar6///41+49+143+144+145+163+164+165+166+167+187+188+189+190+191+192
color orange, active_site
show spheres, active_site

select hotspots, /sym///214+256+284
show sticks, hotspots
color yellow, hotspots

# these dont have cov values
hide cartoon, ////304-306

# draw the covariance onto the protein
inFile = open("Ca_covariance_matrix.txt", "r")
stored.newB = []
for line in inFile.readlines(): stored.newB.append(float(line))
stored.newB2 = []
for x in stored.newB: stored.newB2.append(x)
inFile.close()

# alter 7ar6, b=0.0
# alter 7ar6 and n. CA, b=stored.newB.pop(0)
# cmd.spectrum("b", "magenta_yellow", "7ar6 and n. CA", 0.000001, 0.001)

alter sym, b=0.0
alter sym and n. CA, b=stored.newB2.pop(0)
cmd.spectrum("b", "blue_yellow", "sym", 0.0001, 0.001)

color yellow, hotspots
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

ray 2048,2048
save covariance_overview.png

set_view (\
     0.221945718,   -0.849299788,    0.478984177,\
     0.791667342,    0.443723142,    0.419942230,\
    -0.569197357,    0.285993934,    0.770851195,\
    -0.000116262,   -0.000166513, -146.785186768,\
    -5.780332088,   12.013784409,   15.735269547,\
   116.568534851,  176.786315918,  -20.000000000 )

ray 2048,1024
save covariance_dimer_interface.png
