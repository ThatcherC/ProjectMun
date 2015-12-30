lock munangle to 270-(ship:geoposition:lng + Mun:rotationangle + vang(solarprimevector,Kerbin:position-Mun:position)).

set a1 to munangle.
set time1 to time:seconds.
wait 1.
set a2 to munagle.
set time2 to time:seconds.

print (a2-a1)/(t2-t1).
