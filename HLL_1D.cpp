//coding: utf-8
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;


/////////////////////////declare functions/////////////////////////////////////
double hll_update(double rho[], double u[], double p[], double m[], double E[]);
double speed(double rho[], double u[], double p[], int i);
double decide_dt(double rho[], double u[], double p[]);


///////////////////////////define valiable/////////////////////////////////////
const int nx1 = 100;//-----------------mesh number of space
const int NGHOST = 1;//----------------ghost cells of edge
const int is = NGHOST;//---------------start point of calculation
const int ie = is + nx1 -1;//----------end point of calculation
const int ncells1 = nx1 + 2*NGHOST;//--mesh number include ghost cells
const double Lx = 1.0;//---------------length of tube
const double dx = Lx/nx1;//------------length of a cell
double x[ncells1];//-------------------x-coordinate
double rho[ncells1];//-----------------density
double u[ncells1];//-------------------speed of gas
double p[ncells1];//-------------------pressure
double gamma = 1.4;//------------------effective heat capacity raito
double m[ncells1];//-------------------momentum
double E[ncells1];//-------------------Energy
double rho_new[ncells1], m_new[ncells1], E_new[ncells1];
double t=0, dt, time_out, time_lim=3.0, dt_out=0.1;
double sl,sr;//------------------------max and min eigenvalue


/////////////////////////////////main part/////////////////////////////////////
int main(void)
{
    ofstream ofs1("rho.csv");
    ofstream ofs2("u.csv");
    ofstream ofs3("p.csv");

    //-----------initial condition------------
    //if x<0.5 ... (rho,u,p)=(1.0, 0.0, 1.0)  if x>0.5...(rho,u,p)=(0.125, 0.0, 0.1)
    for (int i = 0 ; i < ncells1 ; i++) {
        x[i] = dx*i;
        u[i] = 0.0;
        if (i <= (ie-is+1)/2)
        {
          rho[i] = 1.0;
          p[i] = 1.0;
        }
        else
        {
          rho[i] = 0.125;
          p[i] = 0.1;
        }
        m[i] = rho[i]*u[i];
        E[i] = p[i]/(gamma-1) + 0.5*u[i]*u[i]*rho[i];
    }

    //------------output initial state---------
    for(int i = is; i <= ie; i++)
    {
      ofs1  <<  rho[i];
      ofs2  <<  u[i];
      ofs3  <<  p[i];
      if (i<ie){
        ofs1 << ",";
        ofs2 << ",";
        ofs3 << ",";
      }
    }
      ofs1  << "\n";
      ofs2  << "\n";
      ofs3  << "\n";
      time_out = dt_out;


    //-------Solve Riemann Problems----------
    while ( t <= time_lim ) {

        //decide dt
        dt = decide_dt(rho,u,p);

        //update u_new, m_new, E_new
        hll_update(rho,u,p,m,E);

        //deal with initial and end point
        rho_new[0] = rho_new[1];
        m_new[0] = m_new[1];
        E_new[0] = E_new[1];
        rho_new[ie+1] = rho_new[ie];
        m_new[ie+1] = m_new[ie];
        E_new[ie+1] = E_new[ie];

        //update rho, m, E, u, p
        for (int i = 0; i < ncells1; i++) {
            rho[i] = rho_new[i];
            m[i] = m_new[i];
            E[i] = E_new[i];
            u[i] = m[i]/rho[i];
            p[i] = (gamma-1)*(E[i] - 0.5*m[i]*m[i]/rho[i]);
        }

        //update time and output
        if (t >= time_out)
        {
          for(int i = is; i <= ie; i++)
          {
            ofs1  <<  rho[i];
            ofs2  <<  u[i];
            ofs3  <<  p[i];
            if (i<ie){
              ofs1 << ",";
              ofs2 << ",";
              ofs3 << ",";
            }
          }
            ofs1  << "\n";
            ofs2  << "\n";
            ofs3  << "\n";
            time_out += dt_out;
          }
        t = t + dt;
    }
}


//////////////////////update U=(rho,m,E) using HLL scheme///////////////////////
double hll_update(double rho[], double u[], double p[], double m[], double E[])
{
    //f_rho, f_u, f_p are elements of flux F_{i} (i=1,2,3...)
    double f_rho[ncells1];
    double f_u[ncells1];
    double f_p[ncells1];
    //fc_rho, fc_u, fc_p are elements of flux F_{i} (i=1/2, 3/2, 5/2...)
    double fc_rho[ncells1];
    double fc_u[ncells1];
    double fc_p[ncells1];
    //define FL, FR, F_HLL
    double FL_rho, FL_u, FL_p;
    double FR_rho, FR_u, FR_p;
    double F_HLL_rho, F_HLL_u, F_HLL_p;

    //-----------input---------------
    for (int i = 0; i < ncells1; i++) {
        f_rho[i] = m[i];
        f_u[i] = rho[i]*u[i]*u[i] + p[i];
        f_p[i] = u[i]*(E[i] + p[i]);
        fc_rho[i] = 0.0;
        fc_u[i]= 0.0;
        fc_p[i] = 0.0;
    }

    //----------calculate fc---------
    for (int i = 0; i <= ie; i++) {
          speed(rho, u, p, i);//---------------update sl,sr
          //decide FL, FR, F_HLL
          FL_rho = f_rho[i];
          FL_u = f_u[i];
          FL_p = f_p[i];
          FR_rho = f_rho[i+1];
          FR_u = f_u[i+1];
          FR_p = f_p[i+1];
          F_HLL_rho = (sr*FL_rho - sl*FR_rho + sl*sr*(rho[i+1]-rho[i]))/(sr-sl);
          F_HLL_u = (sr*FL_u - sl*FR_u + sl*sr*(m[i+1]-m[i]))/(sr-sl);
          F_HLL_p = (sr*FL_p - sl*FR_p + sl*sr*(E[i+1]-E[i]))/(sr-sl);

          //update fc
          if (sl>0)
          {
            fc_rho[i] = FL_rho;
            fc_u[i] = FL_u;
            fc_p[i] = FL_p;
          }
          else if(sl<=0 && sr>=0)
          {
            fc_rho[i] = F_HLL_rho;
            fc_u[i] = F_HLL_u;
            fc_p[i] = F_HLL_p;
          }
          else
          {
            fc_rho[i] = FR_rho;
            fc_u[i] = FR_u;
            fc_p[i] = FR_p;
          }
    }

    //--------output---------------
    for ( int i = is; i <= ie; i++) {
        rho_new[i] = rho[i] + dt/dx*(fc_rho[i] - fc_rho[i-1]);
        m_new[i] = m[i] + dt/dx*(fc_u[i] - fc_u[i-1]);
        E_new[i] = E[i] + dt/dx*(fc_u[i] - fc_u[i-1]);
    }

}

////////////////////////////decide max speed S/////////////////////////////////
double speed(double rho[], double u[], double p[], int i)
{
    double rho_l = rho[i];
    double rho_r = rho[i+1];
    double ul = u[i];
    double ur = u[i+1];
    double pl = p[i];
    double pr = p[i+1];
    double c_l = sqrt(gamma*pl/rho_l);
    double c_r = sqrt(gamma*pr/rho_r);
    sl = ul - c_l;
    sr = ur + c_r;
}

////////////////////////////////decide dt//////////////////////////////////////
double decide_dt(double rho[], double u[], double p[])
{
    double dt = 10000, dt_i;
    for (int i = is; i <= ie; i++)
    {
        dt_i = dx/abs(u[i] + sqrt(gamma*p[i]/rho[i]));
        dt = min(dt_i,dt);
    }
    //std::cout<<dt << "\n";
    return dt;
}
