
reinitialize
bg white
space cmyk

fetch 7z3u

as cartoon
set cartoon_rect_length, 1
set cartoon_oval_length, 1
set stick_radius, 0.25
set sphere_scale, 0.5
set transparency, 0.0
set surface_color, grey

#set ambient, 0.5
#dss

color lightorange

hide ////HOH
hide ////CL
hide ////DMS
#hide cartoon
remove hydrogens

select f140loop, ////139-145
show sticks, f140loop
show sticks, ////41+117

select lig, ////RN2`401
show sticks, lig
color grey70, lig

distance /7z3u/A/A/GLY`143/N, /7z3u/C/A/RN2`401/O26`B
distance /7z3u/A/A/SER`144/N, /7z3u/C/A/RN2`401/O26`B
distance /7z3u/A/A/SER`144/OG, /7z3u/C/A/RN2`401/O26`B
distance /7z3u/A/A/CYS`145/N, /7z3u/C/A/RN2`401/O26`B

color atomic, (not elem C)

deselect
set ray_trace_mode, 1

set_view (\
     0.958532691,   -0.191656783,    0.210895166,\
    -0.281202018,   -0.756093383,    0.590968132,\
     0.046193030,   -0.625773132,   -0.778639674,\
     0.000109069,   -0.000018097,  -38.771678925,\
   -12.720992088,   12.731163025,    8.608248711,\
    31.667379379,   45.826698303,  -20.000000000 )

#ray 1024,1024
#save monomerA-cartoon.png
