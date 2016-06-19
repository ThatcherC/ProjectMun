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
[ ] Mun's position as function of time
[ ] vector math
[ ] time of flight orbit solver
[ ] position and velocity of orbit at a time - hodograph, vis viva?
[ ] six parameters
*/

double rm;      //mun periapsis
double ta = 0;      //time of arrival
double il = 0;
double ir = 0;
double rl;      //Kerbin periapsis outbound
double rr;      //Kerbin periapsis inbound

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
  rl = (double)p+kerbinRadius;
  printf("Mun periapsis (m): ");
  scanf("%d", &p);
  rm = (double)p+munRadius;

  for(int i = 0; i < 5; i++){
    ta += toSeconds[i] * toa[i];
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
/*
  printf("-\ta : %f\n",a);
  printf("-\tc : %f\n",c);
  printf("-\ts : %f\n",s);
  printf("-\talpha : %f\n",alpha);
  printf("-\tbeta : %f\n",beta);
*/
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
};


MunIntercept getIntercept(double rl, vmml::vector<3,double> Rt, double toa){
  vmml::vector<3,double> Rtm;
  vmml::vector<3,double> Vt;
  vmml::vector<3,double> Vtm;

  double a;
  double e;
  double p;

  MunIntercept intercept;

  //Step 1:
  //Choose an r_T and a t_F, then find the trajectory
  Rtm = Rt - munPosition(toa);

  double rt = Rt.length();                                       //arbitrary rt
  double tf = 1.2*parabolicTime(rl, rt);                        //tf must be greater than parabolic time
  double theta = 3.1;                                           //guess

  for(int x = 0; x< 10; x++){
    //printf("%f\n",theta);
    theta = theta + (tf-tfl(rl,rt,theta))/delTfl(rl,rt,theta);
  }
  a = semimajor(rl,rt, theta);
  e = 1-rl/a;                     //eccentricity
  //solved for all necessary orbital parameters

  //Step 2:
  //Find R_tm and V_tm - already have Rtm
  p = a*(1-e*e);

  Vt = 1.0/rt * (sqrt(muKerbin/p)*e*sin(theta)*Rt + sqrt(muKerbin*p)/rt * i.cross(Rt) );
  Vtm = Vt-munVelocity(ta);

  intercept.Rt = Rt;
  intercept.Rtm = Rtm;
  intercept.Vt = Vt;
  intercept.Vtm = Vtm;

  return intercept;
}


//X_y, Xy : vector
//x_y, xy : magnitude of same vector

int main(){
  //Set up
  printf("\n\nMun Planner v0.1, (June 17, 2016)\n----------------------------\n");
  setParamters();

  MunIntercept I1;
  vmml::vector<3,double> Ra;

  //Step 3:
  //Iterate over angles to guess a good starting intercept angle
  const double RAguessingThreshold = 20000.0*5.0/3.0*1000.0 * munRadius/1737000.0 *0.5;//time .8 for extra good guess

  double phi;
  for(phi = 0; phi < 1.7; phi += .1){
    printf("Phi: %f\n",phi);

    vmml::vector<3,double> RtGuess (-munSOI*sin(munAngle(ta)+phi), munSOI*cos(munAngle(ta)+phi), 0 );
    RtGuess = RtGuess + munPosition(ta);
    I1 = getIntercept(rl,RtGuess,ta);

    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
    //printf("\t-R_a: %f\n", Ra.length());
    if(Ra.length()<RAguessingThreshold){
      break;
    }
  }

  for(int x=0; x<6; x++){
    printf("%d: Ra = %f\n",x,Ra.length());
    I1.Rt = munPosition(ta) - munSOI * I1.Vtm/I1.Vtm.length();
    I1 = getIntercept(rl,I1.Rt,ta);
    Ra = I1.Rtm - I1.Rtm.dot(I1.Vtm)/I1.Vtm.dot(I1.Vtm) * I1.Vtm;
  }

  printf("\n");
  return 0;
}
