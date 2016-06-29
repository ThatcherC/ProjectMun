global nodeNumber = 1.

function nextnode{
  print "------ Executing Node #"+nodeNumber +"------".
  set nd to nextnode().

  lock steering to nd:deltav.
  set dv to nd:deltav:mag.
  print "       Delta-V: "+dv.

  //work out burn time eqn
  list engines in englist.
  set eng to englist[0].

  eng:activate().
  print "       Engine Isp: "+eng:isp.

  set i to eng:isp*9.81.
  set f to eng:maxthrust/i.

  set m0 to SHIP:MASS.
  set m2 to m0*constant:e^(-dv/i).
  set m1 to sqrt(m0*m2).

  set t0 to (m0-m1)/f.
  set t1 to (m1-m2)/f.

  print "       tO = " + t0.
  print "       Total Burn = " + (t0+t1).

  wait nd:eta-t0.
  lock throttle to 1.
  wait t0+t1.
  lock throttle to 0.

  unlock steering.
  wait 5.
  remove node1.
}
