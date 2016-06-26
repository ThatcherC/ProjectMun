declare parameter a.
declare parameter e.      //may be unnecessary
declare parameter angle.  //angle of periapsis relative to +x
declare parameter timeTMI.   //time of insertion burn

//current version assumes circular orbit with 0 inclination

function mod{     //computes n%d
  parameter n.
  parameter d.

  if n-d<0{
    return n.
  }
  return mod(n-d,d).
}

//wait until we reach correct angle
lock shipangle to ship:orbit:lan+ship:orbit:argumentofperiapsis+ship:orbit:trueanomaly.
set waitangle to mod(angle-mod(shipangle,360)+360,360).
set waitangle to waitangle*constant:degtorad.

print "Warping to proper angle...".
warpto(time:seconds + waitangle/ship:velocity:orbit:mag*(ship:altitude+kerbin:radius) - 10).
print "...Done".
wait 2.

wait until mod(shipangle,360)>=angle.
print "Reached proper angle.".

set period to ship:orbit:period.
set TMIwait to timeTMI-time:seconds.
set deltaT to mod(TMIwait,period).

//work out math for delta-v for timing burn
//use vis-viva eqn and orbital period equation

set periV to sqrt(kerbin:mu*(2/(ship:orbit:periapsis+kerbin:radius) - 1/ship:orbit:semimajoraxis)).

set newV to sqrt(kerbin:mu*(2/(ship:orbit:periapsis+kerbin:radius)-(kerbin:mu*(period+deltaT)^2/4/constant:PI^2)^(-1/3))).

print "Creating timing burn node...".
set deltaV1 to newV-periV.
set node1 to node(timeTMI-period-deltaT,0,0,deltaV1).
add node1.

//-------Timing burn----------

print "Warping to timing burn...".
WARPTO(timeTMI-period-deltaT-30).
print "...Done".
print "".
lock steering to ship:prograde.

//work out burn time eqn
list engines in a.
set eng to a[0].
set i to eng:isp*9.81.
set f to eng:maxthrust/i.

set m0 to SHIP:MASS.
set m2 to m0*e^(i/deltaV1).
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


//[x] 2. Approach desired periapsis angle
//[ ] 3. Perform burn so that craft returns to periapsis at t_a-tfl
//[ ] 4. Wait until periapsis
//[ ] 5. Burn so that velocity matches calculated value
