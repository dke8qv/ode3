///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d = 0.075;              // meters (correct)
  pars.b = 1.6e-4 * pars.d;    // correct formula
  pars.c = 0.25 * pars.d * pars.d; // correct formula

  void *p_par = (void*) &pars;

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double vPitch = 0;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here
  



  // 2D ODE system for baseball flight with drag
  // State: y[0]=x, y[1]=vx, y[2]=z, y[3]=vz


  auto f_x = [](double t, const vector<double>& y, void* p)->double {
    return y[1];
  };
  auto f_z = [](double t, const vector<double>& y, void* p)->double {
    return y[3];
  };
  auto f_vx = [](double t, const vector<double>& y, void* p)->double {
    Params* pars = (Params*)p;
    double vx = y[1], vz = y[3];
    double v = sqrt(vx*vx + vz*vz);
    if (v < 1e-12) return 0.0;
    double F = pars->b * v + pars->c * v * v;   // magnitude of drag force
    return -(F/pars->m) * (vx / v);
  };
  auto f_vz = [](double t, const vector<double>& y, void* p)->double {
    Params* pars = (Params*)p;
    double vx = y[1], vz = y[3];
    double v = sqrt(vx*vx + vz*vz);
    double drag = 0.0;
    if (v >= 1e-12) {
      double F = pars->b * v + pars->c * v * v;
      drag = (F/pars->m) * (vz / v);
    }
    return -drag - pars->g;
  };

  vector<pfunc_t> f(4);
  f[0] = f_x;
  f[1] = f_vx;
  f[2] = f_z;
  f[3] = f_vz;


  // Function to compute z(xend) for a given v0

  auto z_at_plate = [&](double v0)->double {

    vector<double> y(4);
    double th = theta0 * M_PI / 180.0;
    y[0] = 0.0;
    y[1] = v0 * cos(th);
    y[2] = z0;
    y[3] = v0 * sin(th);

    double t = 0.0;
    double tmax = 2.0;
    int nsteps = 6000;
    double h = tmax / nsteps;

    double x_prev = y[0], z_prev = y[2];

    for (int i = 0; i < nsteps; i++) {
      x_prev = y[0];
      z_prev = y[2];

      y = RK4StepN(f, y, t, h, p_par);
      t += h;

      if (y[2] < 0) return -1e9;  // hits ground
      if (y[0] >= xend) {
        double x_now = y[0], z_now = y[2];
        double frac = (xend - x_prev) / (x_now - x_prev);
        return z_prev + frac * (z_now - z_prev);
      }
    }
    return -1e9;
  };


// Root finding for vPitch: want z(xend)=0.9 m

  double target_z = 0.9;
  double v_lo = 5, v_hi = 80;

  double e_lo = z_at_plate(v_lo) - target_z;
  double e_hi = z_at_plate(v_hi) - target_z;

  // expand bracket if needed
  for (int i = 0; i < 20 && e_lo * e_hi > 0; i++) {
    v_hi += 20;
    e_hi = z_at_plate(v_hi) - target_z;
  }

  // bisection iterations
  for (int it = 0; it < 50; it++) {
    double v_mid = 0.5 * (v_lo + v_hi);
    double e_mid = z_at_plate(v_mid) - target_z;
    if (e_lo * e_mid <= 0) {
      v_hi = v_mid;
      e_hi = e_mid;
    } else {
      v_lo = v_mid;
      e_lo = e_mid;
    }
  }

  vPitch = 0.5 * (v_lo + v_hi);




  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

