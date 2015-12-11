#include <iostream>
#include <stdint.h>
#include <cmath>
#include "stdio.h"
#include "stdlib.h"
#include <string>

// g++ integrator.cpp -o integrator -std=c++0x

using namespace std;

double long energy;
double long radius;
double long velocity;

long double function(double r){
  return 1/sqrt(2*(-278500+3.5316e12/r - 3135L*3135L*680000L*680000L/(2*r*r)));
}

// integrator lower upper steps
int main(int argc, char **argv){
  if(argc != 7){
    printf("\nInvalid number of arguments! Expected arguments:\n");
    printf("Energy (J), radius (m), velocity (m/s), lower bound, upper bound, steps\n\n");
    return 0;
  }
  printf("%i",argc);
  energy = atof(argv[1]);
  radius = atof(argv[2]);
  velocity = atof(argv[3]);
  long double lower = atof(argv[4]);
  long double upper = atof(argv[5]);
  float steps = atof(argv[6]);

  long double dr = (upper-lower)/steps;
  printf("\n");
  printf("Lower Bound: %Le\n",lower);
  printf("Upper bound: %Le\n",upper);
  printf("dr:          %Le\n",dr);
  printf("-------------------------\n");
  long double total = 0;
  long double prev = function(lower);
  long double next;

  for(int i = 1; i < steps+1; i++){
    next = function(lower+i*dr);
    total += (prev+next)*dr/2;
    prev = next;
  }
  printf("%Le s\n",total);
  printf("%ih %im %is\n\n",(int)total/3600,((int)total%3600)/60,(int)total%60);
  return 0;
}
