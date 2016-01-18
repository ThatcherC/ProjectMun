#set up
#set terminal png size 1500,1125 crop font "Arial" 18
unset colorbox
unset ztics
unset key

set xtics out 10
set ytics out .1

set xrange [840:900]
set yrange [.8:1.8]
set hidden3d front
set zrange [-1e19:2e19]

set xlabel "Periapsis speed (m/s)" offset 0,3
set ylabel "Periapsis angle\n(radians)" offset -1,0

set view 0,180,1,1
set isosamples 100

#constants
GMk = 3.5316e12
GMm = 6.5138398e10
rsoi = 2429559.1
dm = 12000000.0


#orbital constants
rpm = 210000.0				#Mun Periapsis
rpk = 635000.0  			#Kerbin's Periapsis

#functions
vmtheta(vburn) = vburn*rpm/rsoi		#Tangential velocity @SOI in Mun frame

vmr_sqr(vburn) = vburn*vburn*(1-rpm*rpm/(rsoi*rsoi)) + 2*GMm*(1/rsoi - 1/rpm)

E(v) = v*v/2 - GMm/rpm

vmr_sqr3(vburn) = 2*(E(vburn) + GMm/rsoi - vburn*vburn*rpm*rpm/(rsoi*rsoi))

vmr(vburn) = sqrt(vmr_sqr3(vburn))	#Radial velocity @SOI in Mun frame

#rx: altitude of craft upon entry to Kerbin system
rx(thetam) = sqrt(dm*dm + rsoi*rsoi + 2*dm*rsoi*sin(thetam))

T(thetam) = 542.5 * rx(thetam)/dm 	#Speed of Kerbin entry point

thetam(vp,thetap) = thetap + def_angle(vp) + .21465
def_angle(vp) = 3.14159 - acos ( -(vp*vp*rpm*rpm - GMm*rsoi)/(vp*vp*rpm*rsoi - GMm*rsoi))


#kerbin equations

#-------One option for vt_sqr
#vt_sqr(vp,thetap) = vmr_sqr(vp) + vmtheta(vp)**2 + T(thetam(vp,thetap))**2 + \
#	4*vmr(vp)*vmtheta(vp)*cos(thetam(vp,thetap))*sin(thetam(vp,thetap))+ \
#	2*T(thetam(vp,thetap))*(vmr(vp)*cos(thetam(vp,thetap)) +  \
#	vmtheta(vp)*sin(thetam(vp,thetap)))

vtheta(vp,thetap) = vmr(vp)*cos(thetam(vp,thetap)) - \
		    vmtheta(vp)*sin(thetam(vp,thetap)) + \
		    T(thetam(vp,thetap))

#-------Another option, more modifiable
vr(vp,thetap) = vmr(vp)*sin(thetam(vp,thetap)) + \
		    vmtheta(vp)*cos(thetam(vp,thetap))

vt_sqr2(vp,thetap) = vtheta(vp,thetap)**2 + vr(vp,thetap)**2

#THE PLOT
splot (vt_sqr2(x,y)-2*GMk/rx(thetam(x,y)))*rpk*rpk + 2*GMk*rpk - vtheta(x,y)**2*rx(thetam(x,y))**2, \
0 with pm3d title ""
