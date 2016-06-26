declare parameter a.
declare parameter e.      //may be unnecessary
declare parameter angle.  //angle of periapsis relative to +x
declare parameter timeTMI.   //time of insertion burn

//current version assumes circular orbit with 0 inclination

set period to ORBIT:PERIOD.
set TMIwait to timeTMI-time:seconds.
set deltaT to TMIwait%period.

//work out math for delta-v for timing burn
//use vis-viva eqn and orbital period equation

periV = sqrt(muKerbin(2/ORBIT:periapsis - 1/ORBIT:A)).
newV = sqrt(muKerbin*(2/ORBIT:periapsishiehgt - (muKerbin*period^2/4/constant:PI^2)^(1/3))).

deltaV1 = newV-periV.
set node1 to node(timeTMI-period-deltaT,0,0,deltaV1).
add node1.

//-------Timing burn----------

WARPTO(timeTMI-period-deltaT-30).

lock steering to ship:prograde.

//work out burn time eqn
list engines in a.
set eng to a[0].
set f to eng:fuelflow.
set i to eng:isp.

set m0 to SHIP:MASS.
set m2 to m0*e^(isp/deltaV1).
set m1 to sqrt(m0*m2).

set t0 to (m0-m1)/f.
set t1 to (m1-m2)/f.

eng:activate().
wait until time:seconds >= timeTMI-period-deltaT-t0.
lock throttle to 1.
wait until time:seconds >= timeTMI-period-deltaT+t1.
lock throttle to 0.

unlock steering.
wait 5.
//remove node1.

//-------TMI burn-----------

WARPTO(timeTMI-60).
lock steering to ship:prograde.


[x] 2. Approach desired periapsis angle
[ ] 3. Perform burn so that craft returns to periapsis at t_a-tfl
[ ] 4. Wait until periapsis
[ ] 5. Burn so that velocity matches calculated value
