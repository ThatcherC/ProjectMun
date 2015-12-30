//Creates a maneuver burn at the desired angle from mun prograde

declare parameter angle.

lock munangle to 270-(ship:geoposition:lng + Mun:rotationangle + vang(solarprimevector,Kerbin:position-Mun:position)).

set angle to angle*3.1415926/180.

set a1 to munangle.
set time1 to time:seconds.
wait 1.
set a2 to munagle.
set time2 to time:seconds.

set w = (a2-a1)/(t2-t1).
print w.

set dp to angle-a2.
print dp.

//make sure this happens in the future not the past
if dp<0 {
set dp to dp+360.
}.

//divide by rate of phase change to get time
set t to dp/w.

//create the node with the proper delta-v and add it to the map
set mynode to NODE(t2+t,0,0,100).
add mynode.
