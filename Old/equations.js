//
//

//variables
var GM_Mun = 6.5138398e10;	//m^3 s^-2
var SOI_Mun = 2429559.1;  	//meters


//Angular Velocity
//d(theta)/dt = l/(mu*r^2)
//l is angular momentum, mu is reduced mass, r is distance from center of body

function w(v0,r0,r){
	return v0*r0/(r*r);
}

//Radial Velocity
//dr/dt = sqrt(2*(E + GM/r - l^2/2r^2)

function radial_velocity_mun(v0,r0,r){
	var energy = v0*v0/2 - GM_Mun/r0;
	return Math.sqrt(2*(energy + GM_Mun/r - v0*v0*r0*r0/(2*r*r)));
}

