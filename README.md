# ProjectMun
Applying physics to a mission to the Mun in Kerbal Space Program. Most of the information will be on my blog,
but all the relevant programs will be kept here.

###Mission Planning

The algorithm for finding a nice lunar trajectory is complicated, so I chose to
implement it in a Julia notebook. I like Julia for programs for with lots of math,
and the notebook format makes it easy to inspect and visualize function outputs.

* *MunBook.ipynb* - A Jupyter-style Julia notebook. I use the package IJulia to 
display and edit it. It roughly follows the method for determining lunar orbits
found in *Fundamentals of Astrodynamics* by Bate et al. 


* *C++/munplanner.cpp* - This program executes a version of Richard Battin's method
for finding a lunar free return trajectory. The method is taken from his book,
*An Introduction to the Mathematics and Methods of Astrodynamics, Revised Edition*.
Before I started the Julia version, this was my main orbit solving program.

###Autopilots

I will also be using kOS programs to perform all the necessary maneuvers once all the calculations are finished.
The kOS mod allows rockets in the game to programmed, allowing much more precise control than I could manage myself.

* *tothemun.ks* - Takes the output from one of the orbit-solving programs and executes a series of burns
to enter the corresponding orbit.
