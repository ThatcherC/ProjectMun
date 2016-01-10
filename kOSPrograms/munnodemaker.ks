//Creates a maneuver burn at the desired angle from mun prograde

declare parameter angle.

declare parameter speed.

lock munangle to arctan2(mun:position:x,mun:position:z)-arctan2((Kerbin:position-mun:position):x,(Kerbin:position-mun:position):z)+90.

set angle to angle*180/3.1415926.

set a1 to munangle.
set time1 to time:seconds.
wait 1.
set a2 to munangle.
set time2 to time:seconds.

set w to (a2-a1)/(time2-time1).
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
set mynode to NODE(time2+t,0,0,speed).
add mynode.
