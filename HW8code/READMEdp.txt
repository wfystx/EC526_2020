The double pendulum code ask for parameter to start and print tracet of
file with these unique paramters.

Look at the "PlotTrace.txt" to see the formate:

e.g.

Trace_th1=1.00_th2=2.00_v1=2.00_v2=-2.00_tmax= 1000.00_dt=0.001000_tsample=0.010000.dat


=====  Physics equation in dbl_pendulum.cpp ++++ 

/*********************************************************
Program which simulates the motion of a double pendulum
   by solving numerically the Lagrange equations with the
   Runge-Kutta algorithm.  This program writes on file the 
   evolution data. 

see beautiful page: https://en.wikipedia.org/wiki/Double_pendulum

The equation for  first figure  with m1 = m2 = m and L1 = L2 = l
has  Kinetic energy is:   K = v1^2/2m + v2^2/2m

The petential energy is for hights of masses:

h1 = 1 - cos(th1) &  h2 = 2 -  cos(th1) - cos(th2)

is:   V = m g l (1-  cos(th1)) +  m g l (2-  cos(th1) + cos(th2))

g = 10 for approximately the acceleration of gravity

                         *** WARNING ***  
The code is actually  written for the more complicated case
of solid bars in the second figure. Don't worry about this.
The  results are exactly the same in the end as far as the dymamics!

******************************************************************/

