// vterm_mass.cpp
// Part A
// Produces:
//   vt_vs_mass.txt          (ASCII table)
//   vt_mass.root            (ROOT graphs)
//   vt_vs_mass.png          (exported ROOT canvas)


#include "RKn.hpp"
#include "TROOT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGClient.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

struct Params {
    double g;
    double m;
    double air_k;
};

double f_ri(double t, const vector<double> &y, void *params=0){
    (void)t; (void)params;
    return y[1];
}

double f_rj(double t, const vector<double> &y, void *params=0){
    (void)t; (void)params;
    return y[3];
}

double f_vi_air(double t, const vector<double> &y, void *params=0){
    (void)t;
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3]);
    return -p->air_k * v * y[1] / p->m;
}

double f_vj_air(double t, const vector<double> &y, void *params=0){
    (void)t;
    Params *p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3]);
    return -p->air_k * v * y[3] / p->m - p->g;
}

// Numerical terminal velocity estimator.
// Runs until projectile is descending, then checks |dv/dt| < tol.
double estimate_vt(double m, double k, double g){
    Params pars;
    pars.m = m;
    pars.air_k = k;
    pars.g = g;
    void *p_par = (void*)&pars;

    vector<pfunc_t> v_fun(4);
    v_fun[0] = f_ri;
    v_fun[1] = f_vi_air;
    v_fun[2] = f_rj;
    v_fun[3] = f_vj_air;

    vector<double> y(4);
    y[0] = 0;
    y[1] = 0;      // Start falling straight down
    y[2] = 100;    // Start 100 m above ground
    y[3] = 0;

    double t = 0;
    double h = 0.002;
    double tmax = 200;

    double vprev = 0;
    int stable = 0;
    const int need = 2000;
    const double tol = 1e-5;

    while (t < tmax){
        y = RK4StepN(v_fun, y, t, h, p_par);
        t += h;

        double v = sqrt(y[1]*y[1] + y[3]*y[3]);
        double dvdt = fabs((v - vprev) / h);
        vprev = v;

        if (dvdt < tol) stable++;
        else stable = 0;

        if (stable > need) return v;

        if (y[2] < 0) return v; // hit ground
    }
    return vprev;
}

int main(){
    const double g = 9.81;
    const double k = 0.1;

    ofstream fout("vt_vs_mass.txt");
    fout << "# mass(kg)\tvt_num\tvt_analytic\n";

    TGraph tg_num;
    TGraph tg_an;

    tg_num.SetTitle("Terminal Velocity vs Mass;Mass m [kg];vt [m/s]");
    tg_an.SetTitle("Terminal Velocity vs Mass;Mass m [kg];vt [m/s]");

    int N = 60;
    double mmin = 1e-3, mmax = 10.0;

    int ip = 0;
    for (int i=0; i<N; i++){
        double frac = (double)i/(N-1);
        double m = mmin * pow(mmax/mmin, frac);

        double vt_num = estimate_vt(m, k, g);
        double vt_an = sqrt(m * g / k);

        fout << m << "\t" << vt_num << "\t" << vt_an << "\n";

        tg_num.SetPoint(ip, m, vt_num);
        tg_an.SetPoint(ip, m, vt_an);

        cout << "m = " << m << " kg   vt_num = " << vt_num
             << "   vt_an = " << vt_an << endl;
        ip++;
    }

    fout.close();

TCanvas *c = new TCanvas("c","vt vs mass",900,700);

    tg_num.SetMarkerStyle(20);
    tg_num.SetMarkerColor(kBlue);
    tg_num.Draw("AP");

    tg_an.SetLineColor(kRed);
    tg_an.SetLineWidth(2);
    tg_an.Draw("L SAME");

    c->SaveAs("vt_vs_mass.png");

    // save ROOT file
    TFile f("vt_mass.root","RECREATE");
    tg_num.Write("vt_num");
    tg_an.Write("vt_an");
    f.Close();

    cout << "\nWrote vt_vs_mass.png and vt_mass.root\n";
    return 0;
}

