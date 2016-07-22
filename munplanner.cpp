#include <iostream>
#include <stdint.h>
#include <cmath>
#include "stdio.h"
#include "stdlib.h"
#include <string>

#include <vmmlib/vector.hpp>
#include <vmmlib/math.hpp>

// g++ munplanner.cpp -I ./ -o munplanner -std=c++0x

using namespace std;

/*
Things needed:
[x] Mun's position as function of time
[x] vector math
[x] time of flight orbit solver
[x] position and velocity of orbit at a time - hodograph, vis viva?
[ ] six parameters
*/

double desired_rm;      //mun periapsis
double desired_ta = 0;      //time of arrival
double desired_il = 0;
double desired_ir = 0;
double desired_rl;      //Kerbin periapsis outbound
double desired_rr;      //Kerbin periapsis inbound

const double kerbinRadius = 600000; //meters
const double munRadius =    200000;
const double munAltitude = 12000000;
const double muKerbin = 3.5316e12; //m^3 s^-2
const double muMun = 6.5138398e10;
const double munSOI = 2429559.1;  //meters
const vmml::vector<3,double> i (0,0,1);           //inclination vector - just z for 0 inclination

void setParamters(){
  int toa[5];

  //equivalences to seconds for years, days, hours, minutes, and seconds
  int toSeconds[5] = {9203400, 6*3600,3600,60,1};

  int p = 0;

  printf("Time of arrival (y:d:h:m:s): ");
  scanf("%d:%d:%d:%d:%d", &toa[0], &toa[1], &toa[2], &toa[3], &toa[4]);

  printf("Initial Kerbin periapsis (m): ");
  scanf("%d",&p);
  desired_rl = (double)p+kerbinRadius;

  printf("Kerbin return periapsis (m): ");
  scanf("%d",&p);
  desired_rr = (double)p+kerbinRadius;

  printf("Mun periapsis (m): ");
  scanf("%d", &p);
  desired_rm = (double)p+munRadius;

  toa[0] = toa[0]-1;
  toa[1] = toa[1]-1;
  for(int i = 0; i < 5; i++){
    desired_ta += toSeconds[i] * toa[i];
  }
  printf("\n----------------------------\n");
}

double parabolicTime(double rl, double rt){     //return time of flight for parabolic trajectory
  return 1.0/3.0 * sqrt(2*(rt-rl)/muKerbin) * (rt+2*rl);
}

double semimajor(double rl, double rt, double theta){     //returns a, the semimajor axis
  double num = rl*(rl-rt*cos(theta));
  double denom = 2*rl-rt*(1+cos(theta));
  return num/denom;
}

double tfl(double rl, double rt, double theta){
  const double a = semimajor(rl,rt,theta);
  const double c = sqrt(rl*rl+rt*rt-2*rl*rt*cos(theta));
  const double s = .5*(rl+rt+c);
  const double alpha = 2*asin(sqrt(s/(2*a)));
  const double beta = 2*asin(sqrt( (s-c)/(2*a) ));

  return sqrt(a*a*a/muKerbin) * ( (alpha-sin(alpha)) - (beta-sin(beta)) );
}

double delS(double rl, double rt, double theta){
  const double c = sqrt(rl*rl+rt*rt-2*rl*rt*cos(theta));
  return rl*rt*sin(theta)/(2*c);
}

double delAOverA(double rl, double rt, double theta){
  double a = semimajor(rl, rt, theta);
  double c2 = rl*rl+rt*rt-2*rl*rt*cos(theta);

  double num = -2*rt*(a-rl)*sin(theta);
  double denom = rl*rl+c2-rt*rt;

  return num/denom;
}


double delTfl(double rl, double rt, double theta){
  const double a = semimajor(rl,rt,theta);
  const double c = sqrt(rl*rl+rt*rt-2*rl*rt*cos(theta));
  const double s = .5*(rl+rt+c);
  const double alpha = 2*asin(sqrt(s/(2*a)));
  const double beta = 2*asin(sqrt( (s-c)/(2*a) ));

  double dtfl = 3*tfl(rl,rt,theta)*delAOverA(rl,rt,theta)/2;
  dtfl += sqrt(a/muKerbin)*( tan(alpha/2)*(delS(rl,rt,theta) - s*delAOverA(rl,rt,theta)) +
                             tan(beta/2) *(delS(rl,rt,theta) + (s-c)*delAOverA(rl,rt,theta)));
  return dtfl;
}

double munAngle(double time){
  return 1.7 + (542.5/munAltitude)*time;
}

vmml::vector < 3, double > munPosition(double time){
  double angle = munAngle(time);
  vmml::vector< 3, double > p( munAltitude*cos(angle), munAltitude*sin(angle), 0 );
  return p;
}

vmml::vector <3, double> munVelocity(double time){
  double angle = munAngle(time);
  vmml::vector<3, double> p(-542.5*sin(angle), 542.5*cos(angle), 0 );
  return p;
}

struct MunIntercept{
  vmml::vector<3,double> Rt;
  vmml::vector<3,double> Rtm;
  vmml::vector<3,double> Vt;
  vmml::vector<3,double> Vtm;
  double tfl;
  double toa;
  double theta;
};


//2D Orbit type
struct Orbit{
  double a;       //semimajor axis
  double e;       //eccentricity

  double aop;     //angle of periapsis
  double time;    //time of periapsis

  double p;       //parameter

  MunIntercept intercept;
};

Orbit getOrbit(double rl, vmml::vector<3,double> Rt, double theta){
  Orbit out;

  out.a = semimajor(rl,Rt.length(), theta);
  out.e = 1-rl/out.a;                     //eccentricity
  out.p = out.a*(1-out.e*out.e);

  out.aop = atan2(Rt[1],Rt[0])-theta;

  return out;
}

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


MunIntercept getIntercept(double rl, vmml::vector<3,double> Rt, double toa, double tf){
  vmml::vector<3,double> Rtm;
  vmml::vector<3,double> Vt;
  vmml::vector<3,double> Vtm;

  MunIntercept intercept;
  Orbit traj;

  //Step 1:
  //Choose an r_T and a t_F, then find the trajectory
  Rtm = Rt - munPosition(toa);

  double rt = Rt.length();                                       //arbitrary rt
  //Now passing tf as function parameter
  //double tf = 1.2*parabolicTime(rl, rt);                        //tf must be greater than parabolic time
  double theta = 3.1;                                           //guess

  for(int x = 0; x< 10; x++){
    //printf("%f\n",theta);
    theta = theta + (tf-tfl(rl,rt,theta))/delTfl(rl,rt,theta);
  }

  traj = getOrbit(rl,Rt,theta);

  //solved for all necessary orbital parameters

  //Step 2:
  //Find R_tm and V_tm - already have Rtm

  Vt = 1.0/rt * (sqrt(muKerbin/traj.p)*traj.e*sin(theta)*Rt + sqrt(muKerbin*traj.p)/rt * i.cross(Rt) );
  Vtm = Vt-munVelocity(desired_ta);

  intercept.Rt = Rt;
  intercept.Rtm = Rtm;
  intercept.Vt = Vt;
  intercept.Vtm = Vtm;
  intercept.toa = toa;
  //is this tfl the same as tf? test
  intercept.tfl = tfl(rl,rt,theta);
  intercept.theta = theta;

  return intercept;
}

Orbit findOrbit(double periapsis, double toa, double tof){
  MunIntercept I1;
  vmml::vector<3,double> Ra;

  //Step 3:
  //Iterate over angles to guess a good starting intercept angle
  const double RAguessingThreshold = 20000.0*5.0/3.0*1000.0 * munRadius/1737000.0 *0.5;//times .5 for extra good guess
  //printf("Ra threshold: %f\n",RAguessingThreshold);

  double phi;
  for(phi = 0; phi < 1.7; phi += .1){
    vmml::vector<3,double> RtGuess (-munSOI*sin(munAngle(toa)+phi), munSOI*cos(munAngle(toa)+phi), 0 );
    RtGuess = RtGuess + munPosition(toa);
    I1 = getIntercept(periapsis,RtGuess,toa,tof);

    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
    //printf("\t-R_a: %f\n", Ra.length());
    if(Ra.length()<RAguessingThreshold){
      break;
    }
  }

  for(int x=0; x<14; x++){
    //printf("%d: Ra = %f\n",x,Ra.length());
    I1.Rt = munPosition(toa) - munSOI * I1.Vtm/I1.Vtm.length();
    I1 = getIntercept(periapsis,I1.Rt,toa,tof);
    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
  }

  Orbit O1 = getOrbit(desired_rl,I1);
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

//X_y, Xy : vector
//x_y, xy : magnitude of same vector

int main(){
  //Set up
  printf("\n\nMun Planner v0.1, (June 17, 2016)\n----------------------------\n");
  setParamters();

  double t_fl = 1.2*parabolicTime(desired_rl, munAltitude+munSOI);
  double t_fr = 1.2*parabolicTime(desired_rr, munAltitude+munSOI);
  Orbit O1;
  Orbit O2;

  //Estimate for tfl - just has to be greater that parabolic time
  //Keeping tfl in main scope is important because it must be varied later
  //Might want to keep Rtm and other variable out here as well
  O1 = findOrbit(desired_rl, desired_ta, t_fl);   //tof variable is unused - expose later!

  //Step 4: Get ToF through Mun SOI
  //Need Vtm and Rtm
  double tSOI = 10000;//getMunSOItime([speed wrt mun]);

  //Step 6: Vary t_fr so that Munar entry and exit velocities match
  for(int x = 0; x< 10; x++){

    //Step 5: find return orbit
    //findOrbit will need some modifications to accomodate a return trajectory
    O2 = findOrbit(desired_rr, desired_ta+tSOI, t_fr);
    printf("VTM outbound: %f   VTM inbound: %f\n", O1.intercept.Vtm.length(),O2.intercept.Vtm.length());

    if(abs(O1.intercept.Vtm.length()-O2.intercept.Vtm.length())<20){
      break;
    }

    t_fr += 1200;
  }
  for(int x = 0; x < 10; x++){
    double deriv = findOrbit(desired_rr, desired_ta+tSOI, t_fr+1).intercept.Vtm.length()-O2.intercept.Vtm.length();

    t_fr = t_fr + (O1.intercept.Vtm.length() - O2.intercept.Vtm.length())/deriv;
    O2 = findOrbit(desired_rr, desired_ta+tSOI, t_fr);

    printf("VTM outbound: %f   VTM inbound: %f\n", O1.intercept.Vtm.length(),O2.intercept.Vtm.length());
  }

  //Step 7: Vary t_fl (and repeat step 6) so that r_m matches desired value


  printf("\n-------Results:--------\n");

  printf("a: %f\n", O1.a);
  printf("e: %f\n", O1.e);
  printf("AoP: %f\n", O1.aop+2*3.141592653);
  printf("ToP: %f\n", O1.time);
  printf("v: %f\n",sqrt(muKerbin * (2/desired_rl-1/O1.a)));
  printf("\nrun tothemun(%f, %f, %f).\n\n",O1.a,O1.aop+2*3.141592653,O1.time);

  printf("a: %f\n", O2.a);
  printf("e: %f\n", O2.e);
  printf("AoP: %f\n", O2.aop+2*3.141592653);
  printf("ToP: %f\n", O2.time);
  printf("v: %f\n",sqrt(muKerbin * (2/desired_rl-1/O2.a)));

  printf("\n");
  return 0;
}
