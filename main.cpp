//coding: utf-8
//created by Kengo Shibata(Osaka University)
// >> g++ -O3 -std=c++11 main.cpp 
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>      // sqrt()
#include <csignal>    // ISO C/C++ signal() and sigset_t, sigemptyset() POSIX C extensions
#include <cstdint>    // int64_t
#include <cstdio>     // sscanf()
#include <cstdlib>    // strtol
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <exception>  // exception
#include <iomanip>    // setprecision()
#include <fstream>    // ofstream
#include <limits>     // max_digits10
#include <new>        // bad_alloc
#include <string>     // string
#include <sstream>    // ostringstream
#include "athena_arrays.hpp"

//----------グローバル変数を定義
namespace{
    int nx1;//-------------------mesh number of x-space 
    int NGHOST;//----------------ghost cells of the edge
    int is;//--------------------start point of calculation
    int ie;//--------------------end point of calculation
    int ncells1;//---------------mesh number include ghost cells
    double xmin, xmax, xc;
    double dx;//-----------------length of a horizonal cell
    double dy;//-----------------length of a vertical cell
    double gamma = 1.666667;//---effective heat capacity raito
    double tlim, nlim, dt_out;
    double den0_L, p0_L, vx0_L, den0_R, p0_R, vx0_R;
    double sr,sl;
    enum{DEN,MOM,EN};//-----------set as DEN=0, MOM=1, EN=2
    enum{RHO,P,VX};
}

//-------------関数を宣言
std::vector<std::string> split(std::string str, std::string separator);
//void  LoadInputFile();
double decide_dt(AthenaArray<double> &prim);
void initial_condition(AthenaArray<double> &u);
void primitive_valiable(AthenaArray<double> &prim, AthenaArray<double> &u);
void output(AthenaArray<double> &prim, int out_n, double t, int cycle);
void decide_eigenvalues(AthenaArray<double> &prim, int i);
void hll_update(AthenaArray<double> &prim, AthenaArray<double> &u, double dt);
void boundary_condition(AthenaArray<double> &u);


//-------------メインループ
int main(void)
{
    // fileの読み込み
    std::ifstream ifs("input.txt");

    if (ifs.fail()) {
        std::cerr << "Failed to open the input file" << std::endl;
        return -1;
    }
    std::string str;

    //デフォルト値をセット
    xmin=-1.0;xmax=1.0;xc=0.0;nx1=100;den0_L=1.0;p0_L=1.0;vx0_L=0;den0_R=1.0;p0_R=1.0;vx0_R=0;tlim=10;nlim=100;dt_out=1.0;

    while (std::getline(ifs, str)) {
        //全文ifsの一行ずつを収めたstrを=の前後で分割
        std::vector<std::string> ary = split(str,"=");
        // remove spaces
        for (unsigned int i = 0; i < ary.size(); i++) {
            std::string temp = ary[i];
            temp.erase(remove(temp.begin(), temp.end(),' '), temp.end());
            ary[i] = temp;
        }

        if (ary[0]=="xmin"){
            xmin = std::stod(ary[1]);
        } else if (ary[0]=="xmax"){
            xmax = std::stod(ary[1]); 
        } else if (ary[0]=="xc"){
            xc = std::stod(ary[1]); 
        }else if (ary[0]=="nx1"){
            nx1 = std::stoi(ary[1]);
        } else if (ary[0]=="den0_L"){
            den0_L = std::stod(ary[1]);
        } else if (ary[0]=="p0_L"){
            p0_L = std::stod(ary[1]);
        } else if (ary[0]=="vx0_L"){
            vx0_L = std::stod(ary[1]);
        } else if (ary[0]=="den0_R"){
            den0_R = std::stod(ary[1]);
        } else if (ary[0]=="p0_R"){
            p0_R = std::stod(ary[1]);
        } else if (ary[0]=="vx0_R"){
            vx0_R = std::stod(ary[1]);
        }else if (ary[0]=="tlim"){
            tlim = std::stod(ary[1]);
        } else if (ary[0]=="nlim"){
            nlim = std::stod(ary[1]);
        } else if (ary[0]=="dt_out"){
            dt_out = std::stod(ary[1]);
        } else {
            continue;
        }

    }
    ifs.close(); // close the file

    //set parameters
    NGHOST = 1;
    is = NGHOST;
    ie = is + nx1 -1;
    ncells1 = nx1 + 2*NGHOST;
    dx = (xmax-xmin)/nx1;

    //----------define local valiable-----------------
    double t = 0.0;
    double dt;
    double output_time = 0.0;
    int    cycle = 0;
    int    out_n = 0;

    //------set initial condition------------------
    //LoadInputFile();
    AthenaArray<double> u,prim;
    u.NewAthenaArray(3, ncells1);
    prim.NewAthenaArray(3, ncells1);
    initial_condition(u);
    primitive_valiable(prim, u);


    //------output initial condition--------------
    out_n = 0;
    //output_time += dt_out;
    cycle = 0;
    output(prim, out_n, t, cycle);


    //---------time integration--------------------
    while (t < tlim && cycle < nlim)
    {
        dt = decide_dt(prim);
        hll_update(prim, u, dt);
        boundary_condition(u);
        primitive_valiable(prim, u);
        
        //-------------output---------------------
        if (t > output_time)
        {
            out_n += 1;
            output_time += dt_out;
            output(prim, out_n, t, cycle);
            
        }

        //-------update t and cycle----------------
        cycle += 1;
        t += dt;
        std::cout << "cycle = " << cycle << "     time = " << t << "     dt = " << dt << "\n";

    }
}


//----------初期条件をセット（ここは自由に変えていいよ！）
void initial_condition(AthenaArray<double> &u)
{
    for (int i = 0; i < ncells1; i++)
    {
        if (xmin+dx*i < xc){
            u(DEN,i) = den0_L;
            u(MOM,i) = den0_L * vx0_L;
            u(EN,i) = p0_L/(gamma-1) + 0.5 * den0_L*pow(vx0_L,2);
        }else{
            u(DEN,i) = den0_R;
            u(MOM,i) = den0_R * vx0_R;
            u(EN,i) = p0_R/(gamma-1) + 0.5 * den0_R*pow(vx0_R,2);
        }
    }
}

//----------保存量uからPrimitive Valiable（基本変数）を計算
void primitive_valiable(AthenaArray<double> &prim, AthenaArray<double> &u)
{
    for (int i = 0; i < ncells1; i++)
    {
        prim(RHO,i) = u(DEN,i);
        prim(VX,i)  = u(MOM,i)/prim(RHO,i);
        prim(P,i)   = (u(EN,i) - 0.5*prim(RHO,i)*pow(prim(VX,i),2))*(gamma-1);
    }
}

//-----------データ出力
void output(AthenaArray<double> &prim, int out_n, double t, int cycle)
{
    FILE *pfile;
    std::string fname;
    char tlabel[5];

    std::snprintf(tlabel, sizeof(tlabel), "%04d", out_n);
    fname.assign("result_");
    fname.append(tlabel);
    fname.append(".tab");
    pfile = fopen(fname.c_str(), "w");
    
    std::fprintf(pfile, "#at time = %4lf, cycle = %04d\n", t, cycle);
    std::fprintf(pfile, "#i    #x       #rho      #vx      #p\n");

    for(int i = is; i <= ie; i++)
    {
        double rho = prim(RHO,i);
        double vx =  prim(VX,i);
        double p  =  prim(P,i);
        std::fprintf(pfile, "%04d %04.4lf %04.4lf %04.4lf %04.4lf\n",i, xmin+i*dx, rho, vx, p); 
    }
    fclose(pfile);

}

//---------ヤコビアン固有値（速度）の最大・最小値であるsr, slを計算（uに音速を足し引きするだけ）
void decide_eigenvalues(AthenaArray<double> &prim, int i)
{
    double rho_l = prim(RHO,i);
    double rho_r = prim(RHO,i+1);
    double ul    = prim(VX,i);
    double ur    = prim(VX,i+1);
    double pl    = prim(P,i);
    double pr    = prim(P,i+1);
    double c_l   = sqrt(gamma*pl/rho_l);
    double c_r   = sqrt(gamma*pr/rho_r);
    sl = std::min(ul - c_l, ur - c_r);
    sr = std::max(ul + c_l ,ur + c_r);
}

//-----------dt を決定する（「dx/音速」の最小値に0.5をかけたくらい）
double decide_dt(AthenaArray<double> &prim)
{
    double cfl = 0.5;
    double dt = 1e+10, dt_i;
    for (int i = is; i <= ie; i++)
    {
        dt_i = dx/abs(prim(VX,i) + sqrt(gamma*prim(P,i)/prim(RHO,i)));
        dt = std::min(dt_i,dt);
    }
    return dt*cfl;
}

//-----------HLLの心臓部
void hll_update(AthenaArray<double> &prim, AthenaArray<double> &u, double dt)
{
    //define fluxes
    //Fcはセル中心でのフラックス（添え字が整数のやつ）、Ffはセル境界のフラックス（添え字が半整数のやつ）
    AthenaArray<double> Fc, Ff, FL, FR;
    Fc.NewAthenaArray(3, ncells1);
    Ff.NewAthenaArray(3,ncells1);
    FL.NewAthenaArray(1,3);
    FR.NewAthenaArray(1,3);
  
    //Fcを計算
    for (int i = 0; i < ncells1; i++)
    {
        Fc(0,i) = u(MOM,i);
        Fc(1,i) = u(MOM,i)*prim(VX,i) + prim(P,i);
        Fc(2,i) = prim(VX,i)*(u(EN,i) + prim(P,i));
    }

    //FcをもとにFR、FLを決定し、リーマン問題を近似的に解き、Ffを決定
    for (int i = is-NGHOST; i <= ie; i++) {
        //FcからFL,FRを決定する。いわゆるreconstructionというやつ。今は簡単のため一次精度にしてある
        FL(0) = Fc(0,i);
        FL(1) = Fc(1,i);
        FL(2) = Fc(2,i);
        FR(0) = Fc(0,i+1);
        FR(1) = Fc(1,i+1);
        FR(2) = Fc(2,i+1);
        
        double F_HLL_rho = (sr*FL(0) - sl*FR(0) + sl*sr*(prim(RHO,i+1)-prim(RHO,i)))/(sr-sl);
        double F_HLL_u   = (sr*FL(1) - sl*FR(1) + sl*sr*(u(MOM,i+1)-u(MOM,i)))/(sr-sl);
        double F_HLL_p   = (sr*FL(2) - sl*FR(2) + sl*sr*(u(EN,i+1)-u(EN,i)))/(sr-sl);

        //HLL法でFf（F_{i+1/2}, F_{i-1/2}に対応）を決定
        decide_eigenvalues(prim, i);//---ここで、場合分けに使うパラメータsr,slを求める
        if (sl>0)
        {
        Ff(0,i) = FL(0);
        Ff(1,i) = FL(1);
        Ff(2,i) = FL(2);
        }
        else if(sl<=0 && sr>=0)
        {
        Ff(0,i) = F_HLL_rho;
        Ff(1,i) = F_HLL_u;
        Ff(2,i) = F_HLL_p;
        }
        else
        {
        Ff(0,i) = FR(0);
        Ff(1,i) = FR(1);
        Ff(2,i) = FR(2);
        }
    }

    //求めたFfから保存量uをアップデート
    for ( int i = is; i <= ie; i++) {
        u(0,i) -= dt/dx*(Ff(0,i) - Ff(0,i-1));
        u(1,i) -= dt/dx*(Ff(1,i) - Ff(1,i-1));
        u(2,i) -= dt/dx*(Ff(2,i) - Ff(2,i-1));
    }

}

//--------------境界条件
void boundary_condition(AthenaArray<double> &u)
{
    //自由端境界
    for (int i =0 ; i< NGHOST; i++){
        u(0,is-i) = u(0,is-i+1);
        u(1,is-i) = u(1,is-i+1);
        u(2,is-i) = u(2,is-i+1);
        u(0,ie+i) = u(0,ie+i-1);
        u(1,ie+i) = u(1,ie+i-1);
        u(2,ie+i) = u(2,ie+i-1);
    }
}

//----------（データ読み込み用）
std::vector<std::string> split(std::string str, std::string separator) {
    if (separator == "") return {str};
    std::vector<std::string> result;
    std::string tstr = str + separator;
    unsigned long l = tstr.length(), sl = separator.length();
    std::string::size_type pos = 0, prev = 0;
    
    for (;pos < l && (pos = tstr.find(separator, pos)) != std::string::npos; prev = (pos += sl)) {
        result.emplace_back(tstr, prev, pos - prev);
    }
    
    return result;
}

//----------（データ読み込み用）
//int LoadInputFile()
//void LoadInputFile()
//{
//    std::ifstream ifs("iput.txt");
//
//    if (ifs.fail()) {
//        std::cerr << "Failed to open the input file" << std::endl;
//        //return -1;
//    }
//    std::string str;
//
//    //デフォルト値をセット
//    xmin=-1.0;xmax=1.0;xc=0.0;nx1=100;den0_L=1.0;p0_L=1.0;vx0_L=0;den0_R=1.0;p0_R=1.0;vx0_R=0;tlim=10;nlim=100;dt_out=1.0;
//
//    while (std::getline(ifs, str)) {
//        //全文ifsの一行ずつを収めたstrを=の前後で分割
//        std::vector<std::string> ary = split(str,"=");
//        // remove spaces
//        for (unsigned int i = 0; i < ary.size(); i++) {
//            std::string temp = ary[i];
//            temp.erase(remove(temp.begin(), temp.end(),' '), temp.end());
//            ary[i] = temp;
//        }
//
//        if (ary[0]=="xmin"){
//            xmin = std::stod(ary[1]);
//        } else if (ary[0]=="xmax"){
//            xmax = std::stod(ary[1]); 
//        } else if (ary[0]=="xc"){
//            xc = std::stod(ary[1]); 
//        }else if (ary[0]=="nx1"){
//            nx1 = std::stoi(ary[1]);
//        } else if (ary[0]=="den0_L"){
//            den0_L = std::stod(ary[1]);
//        } else if (ary[0]=="p0_L"){
//            p0_L = std::stod(ary[1]);
//        } else if (ary[0]=="vx0_L"){
//            vx0_L = std::stod(ary[1]);
//        } else if (ary[0]=="den0_R"){
//            den0_R = std::stod(ary[1]);
//        } else if (ary[0]=="p0_R"){
//            p0_R = std::stod(ary[1]);
//        } else if (ary[0]=="vx0_R"){
//            vx0_R = std::stod(ary[1]);
//        }else if (ary[0]=="tlim"){
//            tlim = std::stod(ary[1]);
//        } else if (ary[0]=="nlim"){
//            nlim = std::stod(ary[1]);
//        } else if (ary[0]=="dt_out"){
//            dt_out = std::stod(ary[1]);
//        } else {
//            continue;
//        }
//
//    }
//    ifs.close(); // close the file
//
    //set parameters
//    NGHOST = 1;
//    is = NGHOST;
//    ie = is + nx1 -1;
//    ncells1 = nx1 + 2*NGHOST;
//    dx = (xmax-xmin)/nx1;
//}


