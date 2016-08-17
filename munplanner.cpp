#include <iostream>
#include <stdint.h>
#include <cmath>
#include "stdio.h"
#include "stdlib.h"
#include <string>

#include <vmmlib/vector.hpp>
#include <vmmlib/math.hpp>

#include <defs.cpp>
#include <basicfunctions.cpp>

// g++ munplanner.cpp -I ./ -o munplanner -std=c++0x

using namespace std;

double desired_rm;      //mun periapsis
double desired_ta = 0;      //time of arrival
double desired_il = 0;
double desired_ir = 0;
double desired_rl;      //Kerbin periapsis outbound
double desired_rr;      //Kerbin periapsis inbound



void setParamters(){
  int toa[5];

  //equivalences to seconds for years, days, hours, minutes, and seconds
  int toSeconds[5] = {9203400, 6*3600,3600,60,1};

  int p = 0;

  printf("Time of arrival (y:d:h:m:s):  ");
  scanf("%d:%d:%d:%d:%d", &toa[0], &toa[1], &toa[2], &toa[3], &toa[4]);

  printf("Initial Kerbin periapsis (m): ");
  scanf("%d",&p);
  desired_rl = (double)p+kerbinRadius;

  printf("Kerbin return periapsis (m):  ");
  scanf("%d",&p);
  desired_rr = (double)p+kerbinRadius;

  printf("Mun periapsis (m):            ");
  scanf("%d", &p);
  desired_rm = (double)p+munRadius;


  toa[0] = toa[0]-1;
  toa[1] = toa[1]-1;
  for(int i = 0; i < 5; i++){
    desired_ta += toSeconds[i] * toa[i];
  }

  //desired_ta = 0;
  printf("\n----------------------------\n");
}


//Calculates the shape (a, e) of an orbit that satisfies the arguments
//Works for both inbound and outbound
Orbit getOrbit(double rl, vmml::vector<3,double> Rt, double theta){
  Orbit out;

  out.a = semimajor(rl,Rt.length(), theta);
  out.e = 1-rl/out.a;                     //eccentricity
  out.p = out.a*(1-out.e*out.e);

  //find angle between between Rt and Mun's center, and add that to munAngle
  vmml::vector<3,double> munPos = munPosition(desired_ta);
  double smallAngle = acos( Rt.dot(munPos) / (Rt.length()*munAltitude) );
  out.aop = munAngle(desired_ta)+smallAngle-theta;;

  return out;
}

//This function is very likely known to work for outbound
//Returns the position and velocity conditions at the Munar patch point
//Assumes an outbound trajectory
MunIntercept getIntercept(double rl, vmml::vector<3,double> Rt, double toa, double angle, int situation){
  vmml::vector<3,double> Rtm;
  vmml::vector<3,double> Vt;
  vmml::vector<3,double> Vtm;

  MunIntercept intercept;
  Orbit traj;

  //Step 1:
  //Choose an r_T and a t_F, then find the trajectory
  Rtm = Rt - munPosition(toa);

  double rt = Rt.length();

  traj = getOrbit(rl,Rt,angle);

  //solved for all necessary orbital parameters

  //Step 2:
  //Find R_tm and V_tm - already have Rtm
  //Comes from IttMaMoA pg. 445

  //should change angle for inbound orbits
  angle = angle - situation*(2*angle);

  Vt = 1.0/rt * (sqrt(muKerbin/traj.p)*traj.e*sin(angle)*Rt + sqrt(muKerbin*traj.p)/rt * i.cross(Rt) );
  Vtm = Vt-munVelocity(toa);

  intercept.Rt = Rt;
  intercept.Rtm = Rtm;
  intercept.Vt = Vt;
  intercept.Vtm = Vtm;
  intercept.toa = toa;
  //is this tfl the same as tf? test
  intercept.tfl = tfl(rl,rt,angle);
  intercept.theta = angle;

  return intercept;
}

//Returns an orbit from Kerbin to Mun ending in the Mun's center
Orbit findOrbit(double periapsis, double toa, double angle, int situation){
  MunIntercept I1;
  vmml::vector<3,double> Ra;

  //Step 3:
  //Iterate over angles to guess a good starting intercept angle
  const double RAguessingThreshold = 20000.0*5.0/3.0*1000.0 * munRadius/1737000.0 *0.5;//times .5 for extra good guess

  double phi;
  //check one quarter of SOI for outbound, check other half for inbound
  for(phi = 0+situation*1.57; phi < 1.57+situation*1.57; phi += .1){
    //Intercept point relative to Mun
    vmml::vector<3,double> RtGuess (-munSOI*sin( munAngle(toa)+phi ), munSOI*cos( munAngle(toa)+phi ), 0 );

    RtGuess = RtGuess + munPosition(toa);


    I1 = getIntercept(periapsis, RtGuess, toa, angle, situation);

    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
    //printf("Phi: %f  R_a: %f\n", phi, Ra.length());
    if(Ra.length()<RAguessingThreshold){
      break;
    }
  }

  //TODO: This might be wayyy to many iterations
  for(int x=0; x<14; x++){
    //printf("%d: Ra = %f\n",x,Ra.length());
    I1.Rt = munPosition(toa) - munSOI * I1.Vtm/I1.Vtm.length();
    I1 = getIntercept(periapsis, I1.Rt, toa, angle, situation);
    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
  }

  Orbit O1 = getOrbit(desired_rl, I1.Rt, angle);
  O1.intercept = I1;
  return O1;
}



double getMunSOItime(double v){
  double a_h = 1.0/(v*v/muMun - 2.0/munSOI);
  //What is rm? desired_rm?? Check section 9.1
  double sineNu = 1.0/(1+desired_rm/a_h);
  double H = acosh((1+munSOI/a_h)*sineNu);
  double t = 2*sqrt(a_h*a_h*a_h/muMun)*(sinh(H)/sineNu - H);
  return t;
}

double getMunRm(MunIntercept munBound, MunIntercept earthBound){
  double nu = asin( munBound.Vtm.cross(earthBound.Vtm).length() / munBound.Vtm.squared_length())/2;
  //printf("Deflection angle: %f\n", nu);
  double a_h = 1.0/(munBound.Vtm.squared_length()/muMun - 2.0/munSOI);
  //printf("A_h: %f\n", a_h);
  return a_h * (1.0/sin(nu) - 1);
}

//Not sure if this is the proper way to do this....
double getMunRa(MunIntercept munBound, MunIntercept earthBound){
  double nu = asin( munBound.Vtm.cross(earthBound.Vtm).length() / munBound.Vtm.squared_length())/2;
  //printf("Deflection angle: %f\n", nu);
  double a_h = 1.0/(munBound.Vtm.squared_length()/muMun - 2.0/munSOI);
  //printf("A_h: %f\n", a_h);
  return a_h * 1.0/tan(nu);
}

//Calculates orbits that satisfy the laws of physics and returns the Mun periapsis
Orbit findConsistentOrbits(double thetaFL){
  Orbit O1;
  Orbit O2;

  O1 = findOrbit(desired_rl, desired_ta, thetaFL, OUTBOUND);

  //Step 4: Get ToF through Mun SOI
  //Need Vtm and Rtm
  double tSOI = getMunSOItime(O1.intercept.Vtm.length());
  //printf("Mun SOI time: %f\n", tSOI);

  //Step 6: Vary t_fr so that Munar entry and exit velocities match
  double thetaFR = 3.1;
  for(int x = 0; x< 10; x++){

    //Step 5: find return orbit
    O2 = findOrbit(desired_rr, desired_ta+tSOI, thetaFR, INBOUND);
    //printf("VTM outbound: %f   VTM inbound: %f\n", O1.intercept.Vtm.length(),O2.intercept.Vtm.length());

    if(abs(O1.intercept.Vtm.length()-O2.intercept.Vtm.length())<20){
      break;
    }

    thetaFR -= .1;
  }

  //printf("Iteratively matching entry and exit speeds:\n");
  for(int x = 0; x < 4; x++){
    double deriv = (findOrbit(desired_rr, desired_ta+tSOI, thetaFR+0.001, INBOUND).intercept.Vtm.length()-O2.intercept.Vtm.length())/0.001;

    thetaFR = thetaFR + (O1.intercept.Vtm.length() - O2.intercept.Vtm.length())/deriv;
    O2 = findOrbit(desired_rr, desired_ta+tSOI, thetaFR, INBOUND);

    //printf("VTM outbound: %f   VTM inbound: %f\n", O1.intercept.Vtm.length(),O2.intercept.Vtm.length());
  }
  return O2;
  //return getMunRm(O1.intercept, O2.intercept);
}


double getOutboundTFL(Orbit orbit, double angle){
  double e = orbit.e;
  double E = acos( (e + cos(angle)) / (1+e*cos(angle)) );

  //Here there was a bug! I squashed it. the '-' used to be a '*'!
  return sqrt(orbit.a*orbit.a*orbit.a/muKerbin)*(E - e*sin(E));
}

void printVector(vmml::vector<3,double> v){
  printf("x: %f   y: %f\n", v[0], v[1]);
}

//X_y, Xy : vector
//x_y, xy : magnitude of same vector


int main(){
  //Set up
  printf("\n\nMun Planner v0.1, (June 17, 2016)\n----------------------------\n");
  setParamters();

  Orbit O1;
  Orbit O2;

  //Estimate for burnout-Mun angle
  double thetaFL = 2.9;
  //double thetaFR = 3.1;
  double RMthreshold = 50000;

  //Step 7: Vary t_fl (and repeat step 6) so that r_m matches desired value
  //Limit on this is currently optimized-ish to 100km, 35km, 10km target
  for(int c = 0; c< 30; c++){
    O1 = findOrbit(desired_rl, desired_ta, thetaFL, OUTBOUND);
    double rm = getMunRm(O1.intercept, findConsistentOrbits(thetaFL).intercept);

    if(abs(rm-desired_rm) < RMthreshold){
      printf("ThetaFL: %f      Calculated RM: %f\n", thetaFL, rm);
      break;
    }

    thetaFL -= 0.01;
  }

  //Time to do a Newton's method on thetaFL
  for(int x = 0; x < 4; x++){
    O1 = findOrbit(desired_rl, desired_ta, thetaFL, OUTBOUND);
    O2 = findConsistentOrbits(thetaFL);
    double rm = getMunRm(O1.intercept, O2.intercept);

    double rm2 = getMunRm(O1.intercept, findConsistentOrbits(thetaFL+0.001).intercept);

    double deriv = (rm2-rm)/0.001;

    thetaFL = thetaFL + (desired_rm-rm)/deriv;
    // = findOrbit(desired_rr, desired_ta+tSOI, thetaFR, INBOUND);

    //printf("RM Desired: %f   RM Calculated: %f\n", rm, desired_rm);
  }

  O1 = findOrbit(desired_rl, desired_ta, thetaFL, OUTBOUND);
  O2 = findConsistentOrbits(thetaFL);
  printf("Ra: %f", getMunRa(O1.intercept, O2.intercept));


  /*

  printf("\n-------Results:--------\n");

  printf("a: %f\n", outbound.a);
  printf("e: %f\n", outbound.e);
  printf("AoP: %f (%f°)\n", outbound.aop, outbound.aop*180/pi);
  printf("Outbound TFL: %f\n", getOutboundTFL(outbound,thetaFL));
  printf("v: %f\n",sqrt(muKerbin * (2/desired_rl-1/outbound.a)));
  printf("\nrun tothemun(%f, %f, %f).\n\n",outbound.a,outbound.aop,desired_ta-getOutboundTFL(outbound,thetaFL));

  printf("\n");*/
  return 0;
}
