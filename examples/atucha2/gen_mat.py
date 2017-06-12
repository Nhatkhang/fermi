#!/usr/bin/env python3

n=451
nz=20
directions=["nw","ne","se","sw"];
print(n*4*20+1," 2")
for i in range(n):
  for j in directions:
    for t in range(nz):
      print( "\"c",i+1, "-",j,"-t",t,"\" 1.500  0.400  0.01012  0.080032  0.000  0.135  1.000  0.000   0.000  0.020" , sep='')

    
print("\"refl_rad\" 1.500  0.400  0.01012  0.020032  1.000 ")
print("\"plenum_sup\" 1.500  0.400  0.01012  0.020032  1.000 ")
print("\"plenum_inf\" 1.500  0.400  0.01012  0.020032  1.000 ")
print("\"plenum_inf_inconel\" 1.500  0.400  0.01012  0.020032  1.000")



