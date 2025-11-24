
/// vterm.cpp
/// Part A for ODE3: projectile motion with/without air resistance,
/// energy conservation study, and terminal velocity extraction (single run).
///
/// Build inside ode3/src with:
///   make vterm
/// Run examples:
///   ./vterm                 (defaults)
///   ./vterm -n 2000         (change number of steps in [t0,tmax])
///   ./vterm --noair         (vacuum, energy conservation printed + saved)
///   ./vterm -v 80 -t 30 --noair
///   ./vterm -m 2 -k 0.1     (air on)
///
/// Outputs:
///   vterm_noair.txt   columns: t x y vx vy E/E0-1
///   vterm_air.txt     columns: t x y vx vy v
///   Also saves ROOT file vterm.root with graphs.

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>

using namespace std;

struct Params {
  double g;      ///< acceleration [m/s^2]
  double m;      ///< mass [kg]
  double air_k;  ///< drag constant k in F_drag = -k v^2
};

// y[0]=x, y[1]=vx, y[2]=y, y[3]=vy, indep var x=t

double f_ri(double t, const vector<double> &y, void *params=0){
  (void) t; (void) params;
  return y[1];
}

double f_rj(double t, const vector<double> &y, void *params=0){
  (void) t; (void) params;
  return y[3];
}

// with drag: dv/dt = -(k/m) v |v|
double f_vi_air(double t, const vector<double> &y, void *params=0){
  (void) t;
  Params *p = (Params*)params;
  double v = sqrt(y[1]*y[1] + y[3]*y[3]);
  return -p->air_k * v * y[1] / p->m;
}
double f_vj_air(double t, const vector<double> &y, void *params=0){
  (void) t;
  Params *p = (Params*)params;
  double v = sqrt(y[1]*y[1] + y[3]*y[3]);
  return -p->air_k * v * y[3] / p->m - p->g;
}

// no air: dvx/dt = 0, dvy/dt = -g
double f_vi_noair(double t, const vector<double> &y, void *params=0){
  (void) t; (void) y; (void) params;
  return 0.0;
}
double f_vj_noair(double t, const vector<double> &y, void *params=0){
  (void) t; (void) y;
  Params *p = (Params*)params;
  return -p->g;
}

// stop when y<0
double f_stop(double t, const vector<double> &y, void *params=0){
  (void) t; (void) params;
  return (y[2] < 0) ? 1.0 : 0.0;
}

static void usage(){
  cout << "Usage: ./vterm [options]\n"
       << "  -v <v0>        initial speed m/s (default 100)\n"
       << "  -t <theta>     launch angle degrees (default 45)\n"
       << "  -m <mass>      mass kg (default 10)\n"
       << "  -k <air_k>     drag constant (default 0.1)\n"
       << "  -n <nsteps>    number of steps (default 200)\n"
       << "  --tmax <sec>   max time window (default 20)\n"
       << "  --noair        turn off drag and do energy conservation study\n"
       << "  --plot         show ROOT plots (default off)\n";
}

int main(int argc, char **argv){

  Params pars;
  pars.g = 9.81;
  pars.m = 10.0;
  pars.air_k = 0.1;
  void *p_par = (void*)&pars;

  double theta = 45.0; // deg
  double v0 = 100.0;   // m/s
  int nsteps = 200;
  double t0 = 0.0;
  double tmax = 20.0;
  bool noair = false;
  bool showPlot = false;

  // long options
  static struct option long_opts[] = {
    {"noair", no_argument, 0, 1000},
    {"tmax", required_argument, 0, 1001},
    {"plot", no_argument, 0, 1002},
    {0,0,0,0}
  };

  int opt, long_index=0;
  while ((opt = getopt_long(argc, argv, "v:t:m:k:n:h", long_opts, &long_index)) != -1){
    switch(opt){
      case 'v': v0 = atof(optarg); break;
      case 't': theta = atof(optarg); break;
      case 'm': pars.m = atof(optarg); break;
      case 'k': pars.air_k = atof(optarg); break;
      case 'n': nsteps = atoi(optarg); break;
      case 1000: noair = true; break;
      case 1001: tmax = atof(optarg); break;
      case 1002: showPlot = true; break;
      case 'h': default: usage(); return 0;
    }
  }

  // function list
  vector<pfunc_t> v_fun(4);
  v_fun[0] = f_ri;
  v_fun[2] = f_rj;
  if (noair){
    v_fun[1] = f_vi_noair;
    v_fun[3] = f_vj_noair;
    pars.air_k = 0.0;
  } else {
    v_fun[1] = f_vi_air;
    v_fun[3] = f_vj_air;
  }

  // initial conditions
  vector<double> y(4);
  double th = theta * M_PI/180.0;
  y[0] = 0.0;
  y[1] = v0 * cos(th);
  y[2] = 0.0;
  y[3] = v0 * sin(th);

  cout << "Mode: " << (noair ? "NO AIR (vacuum)" : "AIR ON (drag)") << "\n";
  cout << "v0="<<v0<<" m/s, theta="<<theta<<" deg, m="<<pars.m<<" kg, k="<<pars.air_k<<"\n";
  cout << "nsteps="<<nsteps<<", tmax="<<tmax<<" s\n";

  // We'll do a low-level stepping loop so we can compute energy and speed.
  double h = (tmax - t0)/nsteps;
  double t = t0;

  ofstream fout;
  if (noair) fout.open("vterm_noair.txt");
  else fout.open("vterm_air.txt");

  // initial energy (vacuum definition)
  double E0 = 0.5*pars.m*(y[1]*y[1]+y[3]*y[3]) + pars.m*pars.g*y[2];

  // ROOT TGraphs
  TGraph tg_x, tg_y, tg_vx, tg_vy, tg_E;
  tg_x.SetName("x_vs_t"); tg_x.SetTitle("x(t);t [s];x [m]");
  tg_y.SetName("y_vs_t"); tg_y.SetTitle("y(t);t [s];y [m]");
  tg_vx.SetName("vx_vs_t"); tg_vx.SetTitle("v_{x}(t);t [s];v_{x} [m/s]");
  tg_vy.SetName("vy_vs_t"); tg_vy.SetTitle("v_{y}(t);t [s];v_{y} [m/s]");
  tg_E.SetName("Erel_vs_t"); tg_E.SetTitle("Relative energy drift; t [s]; (E/E0-1)");

  int ip=0;
  while (t <= tmax){
    double vmag = sqrt(y[1]*y[1]+y[3]*y[3]);
    double E = 0.5*pars.m*vmag*vmag + pars.m*pars.g*y[2];
    double Erel = (E0!=0.0) ? (E/E0 - 1.0) : 0.0;

    if (noair){
      fout << t << "\t" << y[0] << "\t" << y[2] << "\t" << y[1] << "\t" << y[3] << "\t" << Erel << "\n";
    } else {
      fout << t << "\t" << y[0] << "\t" << y[2] << "\t" << y[1] << "\t" << y[3] << "\t" << vmag << "\n";
    }

    tg_x.SetPoint(ip,t,y[0]);
    tg_y.SetPoint(ip,t,y[2]);
    tg_vx.SetPoint(ip,t,y[1]);
    tg_vy.SetPoint(ip,t,y[3]);
    tg_E.SetPoint(ip,t,Erel);

    if (f_stop(t,y,p_par)) break;

    y = RK4StepN(v_fun, y, t, h, p_par);
    t += h;
    ip++;
  }
  fout.close();

  cout << "Final speed = " << sqrt(y[1]*y[1]+y[3]*y[3]) << " m/s\n";
  if (noair){
    cout << "Final relative energy drift E/E0-1 = " << tg_E.GetY()[tg_E.GetN()-1] << "\n";
  } else {
    cout << "Note: with drag, mechanical energy decays (not conserved).\n";
  }

  // save graphs to ROOT file
  TFile tf("vterm.root","recreate");
  tg_x.Write(); tg_y.Write(); tg_vx.Write(); tg_vy.Write(); tg_E.Write();
  tf.Close();

  if (showPlot){
    TApplication theApp("App",&argc,argv);
    UInt_t dh = gClient->GetDisplayHeight()/2;
    UInt_t dw = 1.1*dh;
    TCanvas *c = new TCanvas("c","vterm",dw,dh);
    c->Divide(2,2);
    c->cd(1); tg_x.Draw("AL");
    c->cd(2); tg_y.Draw("AL");
    c->cd(3); tg_vx.Draw("AL");
    c->cd(4); tg_vy.Draw("AL");
    c->Update();
    if (noair){
      TCanvas *cE = new TCanvas("cE","energy drift",dw,dh);
      tg_E.Draw("AL");
      cE->Update();
    }
    cout << "Press Ctrl-C to quit.\n";
    theApp.Run();
  }

  return 0;
}
