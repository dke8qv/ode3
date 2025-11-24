///
/// Starter template for second baseball problem
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
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"


using namespace std;

//---------GLOBALS FOR ODE------------------------------------
struct Params {
    double g;
    double Bmag;
    double wx, wy, wz;
};

static Params pars_global;  // global so ODE funcs can see it

double F_of_v(double v){
    return 0.0039 + 0.0058/(1.0 + exp((v - 35.0)/5.0));
}

double speed_global(const vector<double>& y){
    return sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
}

// x'
double fx(double t,const vector<double>& y, void* p){
    return y[1];
}

// y'
double fy(double t,const vector<double>& y, void* p){
    return y[3];
}

// z'
double fz(double t,const vector<double>& y, void* p){
    return y[5];
}

// vx'
double fvx(double t,const vector<double>& y, void* p){
    double vx=y[1], vy=y[3], vz=y[5];
    double v = speed_global(y);

    if(v < 1e-12) return 0.0;

    double Fv = F_of_v(v);
    double aDx = -Fv*v*vx;
    double aMx = pars_global.Bmag*(pars_global.wy*vz - pars_global.wz*vy);

    return aDx + aMx;
}

// vy'
double fvy(double t,const vector<double>& y, void* p){
    double vx=y[1], vy=y[3], vz=y[5];
    double v = speed_global(y);

    if(v < 1e-12) return -pars_global.g;

    double Fv = F_of_v(v);
    double aDy = -Fv*v*vy;
    double aMy = pars_global.Bmag*(pars_global.wz*vx - pars_global.wx*vz);

    return aDy + aMy;
}

// vz'
double fvz(double t,const vector<double>& y, void* p){
    double vx=y[1], vy=y[3], vz=y[5];
    double v = speed_global(y);

    double Fv = F_of_v(v);
    double aDz = -Fv*v*vz;
    double aMz = pars_global.Bmag*(pars_global.wx*vy - pars_global.wy*vx);

    return aDz + aMz - pars_global.g;
}



int main(int argc, char **argv){

  // we have 6 initial conditions for this problem
  // y[0] = y[2] = y[4] = 0;  // init x,y,z
  // y[1] = v0*cos(theta0);   // vx  "x is line towards the plate
  // y[3] = 0;                // vy  "y" is measured as left/right divergence from line to plate
  // y[5] = v0*sin(theta0);   // vz  "z" is vertival measure
  vector<double> y0(6);

  bool showPlot=false;
  // pitches
  // slider ip=0
  // curve ip=1
  // screwball ip=2
  // fast ip=3
  int ip=1;    // default pitch
  int c;
  while ((c = getopt (argc, argv, "p:n")) != -1)
    switch (c) {
    case 'p':
      ip = atoi(optarg);
      showPlot=true;
      break;
    case 'n':
      showPlot=false;
      break;
    }

  TString title;
  if (ip==0){
    cout << "Setting up initial conditions for slider" << endl;
    //SetupSlider(y0);
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
    //SetupCurve(y0);
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
    //SetupScrewball(y0);
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
    //SetupFastball(y0);
  }

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  double xend=60;   // feet
  double yend=0;    // tbd
  double zend=0;    // tbd
  double vxend=0;
  double vyend=0;
  double vzend=0;

  // write code here

  // Physical constants
  const double ft_to_m = 0.3048;
  const double m_to_ft = 1.0/ft_to_m;

  pars_global.g = 9.81;
  pars_global.Bmag = 4.1e-4;   // Fitzpatrick eq 3.43


  // Select pitch parameters
  double v0;
  double theta = 1.0*M_PI/180.0;

  double omega_rpm = 1800.0;
  double omega = omega_rpm*(2*M_PI/60.0);

  double phi_deg;
  if(ip==0)      phi_deg = 0.0;     // slider
  else if(ip==1) phi_deg = 45.0;    // curveball
  else if(ip==2) phi_deg = 135.0;   // screwball
  else           phi_deg = 225.0;   // fastball

  if(ip==3) v0 = 42.47;   // fastball ~95 mph
  else      v0 = 38.01;   // slider, curve, screwball ~85 mph

  double phi = phi_deg*M_PI/180.0;

  pars_global.wx = 0.0;
  pars_global.wy = omega*sin(phi);
  pars_global.wz = omega*cos(phi);

 
  // Initial state
  y0.assign(6,0.0);
  y0[1] = v0*cos(theta);
  y0[5] = v0*sin(theta);


  // Build ODE function array (true function pointers)
  vector<pfunc_t> f(6);
  f[0] = fx;
  f[1] = fvx;
  f[2] = fy;
  f[3] = fvy;
  f[4] = fz;
  f[5] = fvz;


// Create graphs for plotting
TGraph *gyx = new TGraph();  // z vs x
TGraph *gzx = new TGraph();  // y vs x

vector<double> y = y0;

double t=0.0;
double dt=0.0001;
double x_target_m = 60.6*ft_to_m;

double x_prev=y[0], y_prev2=y[2], z_prev=y[4];


// Integration 
for(int i=0; i<500000; i++){

    // convert to feet for plotting
    double x_ft = y[0]*m_to_ft;
    double y_ft = y[2]*m_to_ft;
    double z_ft = y[4]*m_to_ft;

    gyx->SetPoint(gyx->GetN(), x_ft, z_ft);  // z(x)
    gzx->SetPoint(gzx->GetN(), x_ft, y_ft);  // y(x)

    if(y[0] >= x_target_m){
        double frac = (x_target_m - x_prev)/(y[0] - x_prev);

        xend = 60.6;
        yend = (y_prev2 + frac*(y[2]-y_prev2))*m_to_ft;
        zend = (z_prev  + frac*(y[4]-z_prev ))*m_to_ft;

        vxend = y[1]*m_to_ft;
        vyend = y[3]*m_to_ft;
        vzend = y[5]*m_to_ft;
        break;
    }

    x_prev = y[0];
    y_prev2 = y[2];
    z_prev = y[4];

    y = RK4StepN(f,y,t,dt,nullptr);
    t += dt;
}


TString canname;
if(ip==0) canname="Slider";
if(ip==1) canname="Curveball";
if(ip==2) canname="Screwball";
if(ip==3) canname="Fastball";

// Create canvas 
TCanvas *c1 = new TCanvas("c1", canname, 900, 600);

// Style adjustments
gStyle->SetLineWidth(2);

// Solid line for vertical displacement z(x)
gyx->SetLineColor(kBlue);
gyx->SetLineWidth(3);
gyx->SetLineStyle(1);   // solid
gyx->SetTitle(canname + "; x (ft); displacement (ft)");

// Dashed line for horizontal displacement y(x)
gzx->SetLineColor(kRed);
gzx->SetLineWidth(3);
gzx->SetLineStyle(2);   // dashed

// Draw both on same axes
gyx->Draw("AL");
gzx->Draw("L SAME");

// Add legend
auto legend = new TLegend(0.15, 0.75, 0.45, 0.88);
legend->AddEntry(gyx, "Vertical break z(x)", "l");
legend->AddEntry(gzx, "Horizontal break y(x)", "l");
legend->Draw();

// Update canvas
c1->Update();


TString pdfname = canname + ".pdf";
c1->SaveAs(pdfname);





  // to compare to the plots in Fitzpatrick, output your results in **feet**
  // do not change these lines
  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,x) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("********************************\n");

  // plot the trajectory.  See Fitzpatrick for plot details
  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

