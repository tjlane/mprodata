
reinitialize
bg white

load 7ar6.pdb

as cartoon
set cartoon_rect_length, 1
set cartoon_oval_length, 1
set stick_radius, 0.25
set sphere_scale, 0.5
set transparency, 0.5
set surface_color, grey

symexp sym, 7ar6, (7ar6), 0.1
set_name sym01000000, sym
remove not alt ''+A
dss

hide ////HOH
hide ////CL
hide ////DMS
remove hydrogens

color grey, 7ar6
show surface, 7ar6

select active_site, /7ar6///41+49+143+144+145+163+164+165+166+167+187+188+189+190+191+192
color red, active_site
show spheres, active_site

select hotspots, /sym///214+256+284
show sticks, hotspots
#color yellow, hotspots

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
cmd.spectrum("b", "magenta_yellow", "sym", 0.000001, 0.001)

color atomic, (not elem C)


set_view (\
    -0.555719018,    0.156507254,   -0.816507816,\
    -0.011010773,   -0.983419061,   -0.181006759,\
    -0.831299365,   -0.091598041,    0.548227012,\
     0.000000000,    0.000000000, -200.993118286,\
     0.000003815,    0.407508850,    0.000011444,\
   158.464508057,  243.521728516,  -20.000000000 )

deselect
set ray_trace_mode, 1
#ray
