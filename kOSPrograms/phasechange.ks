lock ShipPos to -Kerbin:position.
lock MunPos to Mun:position-Kerbin:position.

set phase1 to arctan2(ShipPos:x,ShipPos:z)-arctan2(MunPos:x,MunPos:z).
set time1 to time:seconds.
wait 1.
set phase2 to arctan2(ShipPos:x,ShipPos:z)-arctan2(MunPos:x,MunPos:z).
set time2 to time:seconds.

print (phase2-phase1)-(time2-time1).