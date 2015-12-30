lock munangle to arctan2(mun:position:x,mun:position:z)-arctan2((Kerbin:position-mun:position):x,(Kerbin:position-mun:position):z)+90.

set a1 to munangle.
set time1 to time:seconds.
wait 1.
set a2 to munangle.
set time2 to time:seconds.

print a2.
print (a2-a1)/(time2-time1).
