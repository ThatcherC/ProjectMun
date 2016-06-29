declare parameter a.
declare parameter angle.  //angle of periapsis relative to +x
declare parameter timeTMI.   //time of insertion burn

run once nextnode.

function mod{     //computes n%d
  parameter n.
  parameter d.

  if n-d<0{
    return n.
  }
  return mod(n-d,d).
}

//current version assumes circular orbit with 0 inclination

set angle to angle*constant:radtodeg.

//wait until we reach correct angle
lock shipangle to ship:orbit:lan+ship:orbit:argumentofperiapsis+ship:orbit:trueanomaly.     //in degrees
set waitangle to mod(angle-mod(shipangle,360)+360,360).                                     //in degrees
set waitangle to waitangle*constant:degtorad.                                               //in radians

print "Warping to proper angle (waiting " + (waitangle/ship:velocity:orbit:mag*(ship:altitude+kerbin:radius) - 15) +"s)".

warpto(time:seconds + waitangle/ship:velocity:orbit:mag*(ship:altitude+kerbin:radius) - 10).
wait waitangle/ship:velocity:orbit:mag*(ship:altitude+kerbin:radius) - 10.

print "...Done".

print "Desired angle: "+angle.
print "Current angle: "+mod(shipangle,360).
wait until mod(shipangle,360)>=angle.
print "Reached proper angle.".
print "Desired angle: "+angle.
print "Current angle: "+mod(shipangle,360).


set period to ship:orbit:period.
set TMIwait to timeTMI-time:seconds.
set deltaT to mod(TMIwait,period).

//Use vis-viva eqn and orbital period eqn to work out timing burn and transfer burn
set periV to sqrt(kerbin:mu*(2/(ship:orbit:periapsis+kerbin:radius) - 1/ship:orbit:semimajoraxis)).
set timingV to sqrt(kerbin:mu*(2/(ship:orbit:periapsis+kerbin:radius)-(kerbin:mu*(period+deltaT)^2/4/constant:PI^2)^(-1/3))).
set transferV to sqrt(kerbin:mu*(2/(ship:orbit:periapsis+kerbin:radius) - 1/a)).

print "Creating timing burn node...".
set deltaV1 to timingV-periV.
set deltaV2 to transferV-timingV.
set node1 to node(timeTMI-period-deltaT,0,0,deltaV1).
set node2 to node(timeTMI,0,0,deltaV2).
add node1.
add node2.

//-------Timing burn----------

print "Warping to timing burn...".
WARPTO(timeTMI-period-deltaT-30).
wait until time:seconds>(timeTMI-period-deltaT-30).
print "...Done".

nextnode().



//-------TMI burn-----------

WARPTO(timeTMI-60).

nextnode().


//[x] 2. Approach desired periapsis angle
//[ ] 3. Perform burn so that craft returns to periapsis at t_a-tfl
//[ ] 4. Wait until periapsis
//[ ] 5. Burn so that velocity matches calculated value
