#set up
set xrange [730:900]
set yrange [.6:2.6]
set hidden3d
set zrange [-2e19:2e19]

set view 0,180,1,1
set isosamples 80

#constants
GMk = 3.5316e12
GMm = 6.5138398e10
rsoi = 2429559.1
dm = 12000000


#orbital constants
rpm = 210000				#Mun Periapsis
rpk = 635000  				#Kerbin's Periapsis

#functions
vmtheta(vburn) = vburn*rpm/rsoi		#Tangential velocity @SOI in Mun frame

vmr_sqr(vburn) = vburn*vburn*(1-rpm*rpm/(rsoi*rsoi)) + 2*GMm*(1/rsoi - 1/rpm)
vmr(vburn) = sqrt(vmr_sqr(vburn))	#Radial velocity @SOI in Mun frame

#rx: altitude of craft upon entry to Kerbin system
rx(thetam) = sqrt(dm*dm + rsoi*rsoi - 2*dm*rsoi*sin(thetam))

T(thetam) = 542.5 * rx(thetam)/dm 	#Speed of Kerbin entry point

thetam(vp,thetap) = thetap + def_angle(vp)
def_angle(vp) = acos ( -(vp*vp*rpm*rpm - GMm*rsoi)/(vp*vp*rpm*rsoi - GMm*rsoi))


#kerbin equations

vt_sqr(vp,thetap) = vmr_sqr(vp) + vmtheta(vp)**2 + T(thetam(vp,thetap))**2 + \
	4*vmr(vp)*vmtheta(vp)*cos(thetam(vp,thetap))*sin(thetam(vp,thetap))+ \
	2*T(thetam(vp,thetap))*(vmr(vp)*cos(thetam(vp,thetap)) +  \
	vmtheta(vp)*sin(thetam(vp,thetap)))

vtheta(vp,thetap) = vmr(vp)*cos(thetam(vp,thetap)) + \
		    vmtheta(vp)*sin(thetam(vp,thetap)) + \
		    T(thetam(vp,thetap))



#THE PLOT
splot (vt_sqr(x,y)-2*GMk/rx(thetam(x,y)))*rpk*rpk + 2*GMk*rpk - vtheta(x,y)**2*rx(thetam(x,y))**2, 0
