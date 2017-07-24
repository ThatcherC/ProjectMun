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
