//Creates a maneuver node to impact the Mun

declare parameter angle.

//ShipPos and MunPos are the vector positions of the ship and the Mun
//relative to Kerbin
lock ShipPos to -Kerbin:position.
lock MunPos to Mun:position-Kerbin:position.

//This gets the current angle between the ship and the Mun
set phase to arctan2(ShipPos:x,ShipPos:z)-arctan2(MunPos:x,MunPos:z).

//dp is difference between current and desired phase
set dp to (phase-angle).

//make sure this happens in the future not the past
if dp<0 {
set dp to dp+360.
}.

//divide by rate of phase change (from phasechange.ks) to get time
set t to dp/.1894.

print t.

//create the node with the proper delta-v and add it to the map
set mynode to NODE(time:seconds+t,0,0,856.3).
add mynode.
