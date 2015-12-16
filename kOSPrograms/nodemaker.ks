//Creates a maneuver node to impact the Mun

lock ShipPos to -Kerbin:position.
lock MunPos to Mun:position-Kerbin:position.

set phase to arctan2(ShipPos:x,ShipPos:z)-arctan2(MunPos:x,MunPos:z).

//replace 62.56 with calculated phase angle
set dp to (phase-180+62.56).
if dp<0 {
set dp to dp+360.
}.

set t to dp/.1894.

print t.

set mynode to NODE(time:seconds+t,0,0,856.3).

add mynode.
