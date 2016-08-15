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

  //printf("Time of arrival (y:d:h:m:s): ");
  //scanf("%d:%d:%d:%d:%d", &toa[0], &toa[1], &toa[2], &toa[3], &toa[4]);

  printf("Initial Kerbin periapsis (m): ");
  scanf("%d",&p);
  desired_rl = (double)p+kerbinRadius;

  printf("Kerbin return periapsis (m): ");
  scanf("%d",&p);
  desired_rr = (double)p+kerbinRadius;

  printf("Mun periapsis (m): ");
  scanf("%d", &p);
  desired_rm = (double)p+munRadius;

  /*
  toa[0] = toa[0]-1;
  toa[1] = toa[1]-1;
  for(int i = 0; i < 5; i++){
    desired_ta += toSeconds[i] * toa[i];
  }
  */
  desired_ta = 0;
  printf("\n----------------------------\n");
}


//Calculates the shape (a, e) of an orbit that satisfies the arguments
//Works for both inbound and outbound
Orbit getOrbit(double rl, vmml::vector<3,double> Rt, double theta){
  Orbit out;

  out.a = semimajor(rl,Rt.length(), theta);
  out.e = 1-rl/out.a;                     //eccentricity
  out.p = out.a*(1-out.e*out.e);

  out.aop = atan2(Rt[1],Rt[0])-theta;

  return out;
}
/*
//Returns a fully specified orbit with the time of periapsis and burn
//Only works for zero inclination orbits.
Orbit getOrbit(double _rl, MunIntercept intercept){
  Orbit out;

  out.a = 1.0/(-intercept.Vt.dot(intercept.Vt)/muKerbin + 2/intercept.Rt.length());
  out.e = 1-_rl/out.a;                     //eccentricity
  out.p = out.a*(1-out.e*out.e);

  out.aop = atan2(intercept.Rt[1],intercept.Rt[0])-intercept.theta;

  out.time = intercept.toa-intercept.tfl;

  return out;
}
*/

//This function is very likely known to work for outbound
//Returns the position and velocity conditions at the Munar patch point
//Assumes an outbound trajectory
MunIntercept getIntercept(double rl, vmml::vector<3,double> Rt, double toa, double angle){
  vmml::vector<3,double> Rtm;
  vmml::vector<3,double> Vt;
  vmml::vector<3,double> Vtm;

  MunIntercept intercept;
  Orbit traj;

  //Step 1:
  //Choose an r_T and a t_F, then find the trajectory
  Rtm = Rt - munPosition(toa);

  double rt = Rt.length();                                       //arbitrary rt

  traj = getOrbit(rl,Rt,angle);

  //solved for all necessary orbital parameters

  //Step 2:
  //Find R_tm and V_tm - already have Rtm

  //This might be adaptable if angle -> 2*PI-angle [true assuming we want hit a periapsis, which we do]
  //Comes from IttMaMoA pg. 445
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
Orbit findOrbitOutbound(double periapsis, double toa, double angle){
  MunIntercept I1;
  vmml::vector<3,double> Ra;

  //Step 3:
  //Iterate over angles to guess a good starting intercept angle
  const double RAguessingThreshold = 20000.0*5.0/3.0*1000.0 * munRadius/1737000.0 *0.5;//times .5 for extra good guess

  double phi;
  for(phi = 0; phi < 1.57; phi += .1){
    //Intercept point relative to Mun
    vmml::vector<3,double> RtGuess (-munSOI*sin( munAngle(toa)+phi ), munSOI*cos( munAngle(toa)+phi ), 0 );

    RtGuess = RtGuess + munPosition(toa);


    I1 = getIntercept(periapsis, RtGuess, toa, angle);

    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
    printf("Phi: %f  R_a: %f\n", phi, Ra.length());
    if(Ra.length()<RAguessingThreshold){
      break;
    }
  }

  //TODO: This might be wayyy to many iterations
  for(int x=0; x<14; x++){
    //printf("%d: Ra = %f\n",x,Ra.length());
    I1.Rt = munPosition(toa) - munSOI * I1.Vtm/I1.Vtm.length();
    I1 = getIntercept(periapsis, I1.Rt, toa, angle);
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

void printVector(vmml::vector<3,double> v){
  printf("x: %f   y: %f\n", v[0], v[1]);
}

//X_y, Xy : vector
//x_y, xy : magnitude of same vector


int main(){
  //Set up
  printf("\n\nMun Planner v0.1, (June 17, 2016)\n----------------------------\n");
  setParamters();

  //Estimate for burnout-Mun angle
  double thetaFL = 3.0;
  double thetaFR = 2.5;
  Orbit O1;
  Orbit O2;

  //Step 7: Vary t_fl (and repeat step 6) so that r_m matches desired value
  for(int c = 0; c< 10; c++){
    O1 = findOrbit(desired_rl, desired_ta, thetaFL);

    //Step 4: Get ToF through Mun SOI
    //Need Vtm and Rtm
    double tSOI = getMunSOItime(O1.intercept.Vtm.length());
    printf("Mun SOI time: %f\n", tSOI);

    //Step 6: Vary t_fr so that Munar entry and exit velocities match
    for(int x = 0; x< 10; x++){

      //Step 5: find return orbit
      //findOrbit will need some modifications to accomodate a return trajectory
      O2 = findOrbit(desired_rr, desired_ta+tSOI, t_fr);
      //printf("VTM outbound: %f   VTM inbound: %f\n", O1.intercept.Vtm.length(),O2.intercept.Vtm.length());

      if(abs(O1.intercept.Vtm.length()-O2.intercept.Vtm.length())<20){
        break;
      }

      t_fr += 1200;
    }
    for(int x = 0; x < 4; x++){
      double deriv = findOrbit(desired_rr, desired_ta+tSOI, t_fr+1).intercept.Vtm.length()-O2.intercept.Vtm.length();

      t_fr = t_fr + (O1.intercept.Vtm.length() - O2.intercept.Vtm.length())/deriv;
      O2 = findOrbit(desired_rr, desired_ta+tSOI, t_fr);

      //printf("VTM outbound: %f   VTM inbound: %f\n", O1.intercept.Vtm.length(),O2.intercept.Vtm.length());
    }

    printf("Calculated R_m: %f\n\n", getMunRm(O1.intercept, O2.intercept));
  }



  printf("\n-------Results:--------\n");

  printf("a: %f\n", O1.a);
  printf("e: %f\n", O1.e);
  printf("AoP: %f\n", O1.aop+2*3.141592653);
  printf("ToP: %f\n", O1.time);
  printf("v: %f\n",sqrt(muKerbin * (2/desired_rl-1/O1.a)));
  printf("\nrun tothemun(%f, %f, %f).\n\n",O1.a,O1.aop+2*3.141592653,O1.time);

  printf("\n");
  return 0;
}
