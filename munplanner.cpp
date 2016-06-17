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

const int kerbinRadius = 600000; //meters
const int munRadius =    200000;
const double muKerbin = 3.5316e12; //m^3 s^-2
const double muMun = 6.5138398e10;
const double munSOI = 2429559.1;  //meters


void setParamters(){
  int toa[5];

  //equivalences to seconds for years, days, hours, minutes, and seconds
  int toSeconds[5] = {9203400, 6*3600,3600,60,1};

  int p = 0;

  printf("Time of arrival (y:d:h:m:s): ");
  scanf("%d:%d:%d:%d:%d", &toa[0], &toa[1], &toa[2], &toa[3], &toa[4]);
  printf("Initial Kerbin periapsis (m): ");
  scanf("%d",&p);
  rl = (double)p;
  printf("Mun periapsis (m): ");
  scanf("%d", &p);
  rm = (double)p;

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
  double denom = 2*rl-rt*(1-cos(theta));
  return num/denom;
}

double tfl(double rl, double rt, double theta){
  double a = semimajor(rl,rt,theta);
  double c = sqrt(rl*rl+rt*rt-2*rl*rt*cos(theta));
  double s = .5*(rl+rt+c);
  double alpha = 2*asin(sqrt(s/(2*a)));
  double beta = 2*asin(sqrt( (s-c)/(2*a) ));

  return sqrt(a*a*a/muKerbin) * ( (alpha-sin(alpha)) - (beta-sin(beta)) );
}


int main(){
  printf("\n\nMun Planner v0.1, (June 17, 2016)\n----------------------------\n");

  setParamters();

  printf("\n");
  return 0;
}
