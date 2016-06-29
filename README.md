# ProjectMun
Applying physics to a mission to the Mun in Kerbal Space Program. Most of the information will be on my blog,
but all the relevant programs will be kept here.

###Mission Planning

The algorithm for finding a nice lunar trajectory is complicated, so I chose to
implement it in C++ rather than kOS language. I have one program to work out the
orbits:

* *munplanner.cpp* - This program executes a version of Richard Battin's method
for finding a lunar free return trajectory. The method is taken from his book,
*An Introduction to the Mathematics and Methods of Astrodynamics, Revised Edition*.

###Autopilots

I will also be using kOS programs to perform all the necessary maneuvers once all the calculations are finished.
The kOS mod allows rockets in the game to programmed, allowing much more precise control than I could manage myself.

* *tothemun.ks* - Takes the output from munplanner.cpp and executes a series of burns
to enter into the orbit munplanner chooses.
