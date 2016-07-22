

const double kerbinRadius = 600000; //meters
const double munRadius =    200000;
const double munAltitude = 12000000;
const double muKerbin = 3.5316e12; //m^3 s^-2
const double muMun = 6.5138398e10;
const double munSOI = 2429559.1;  //meters
const vmml::vector<3,double> i (0,0,1);           //inclination vector - just z for 0 inclination


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
