#include<iostream>
#include<fstream>
#include<cmath>
#include"Vector.h"
#include"Random64.h"

using namespace std;

const double g=100, K=1e4;
const double Kx=1,Ky=1,Kz=1;               /*Spring's Constant in each axis*/
const double Lx=100,Ly=100,Lz=100;         /*Lenght for parametrize the sides of the box*/
const int Nx=3,Ny=3,Nz=3;                  /*Number of springs in each axis*/
const int N2D=Nx*Ny ,N3D=Nx*Ny*Nz;         /*Total number of springs*/
const double Gamma=50, Kcundall=10, MU=0.4;/*Contact Force*/

const double Xi=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Chi=-0.06626458266981849;

class Body;
class Collider;
class Oscillator;
/*--------------------------------------------------CLASS BODY------------------------------------------------------------------------*/
class Body{
private:
  double Rx,Ry,Rz,Vx,Vy,Vz,Fx,Fy,Fz,m,r,xBreak,yBreak,zBreak; 
public:
  void Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
	     double m0,double r0,double x0Break,double y0Break,double z0Break);
  /*Move the Position and Velocity in each component*/
  void Move_x(double dt, double Constante);
  void Move_y(double dt, double Constante);
  void Move_z(double dt, double Constante);
  void Move_Vx(double dt, double Constante);
  void Move_Vy(double dt, double Constante);
  void Move_Vz(double dt, double Constante);
  /*Inline Functions*/
  double Getx(void){return Rx;};
  double Gety(void){return Ry;};
  double Getz(void){return Rz;};
  /*Add/Erases Forces and Torsion*/
  void EraseForce(void);
  void AddForce_x(double F0x);
  void AddForce_y(double F0y);
  void AddForce_z(double F0z);
  /*Draw the Trajectory of the Particle*/
  void Draw(void);
  void Draw3D(void);
  friend class Collider;
  friend class Oscillator;
};
/*-------------------------------------------------BODY'S FUNCTIONS-------------------------------------------------------------------*/
void Body::Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
		 double m0,double r0,double x0Break,double y0Break,double z0Break)
{Rx=Rx0; Ry=Ry0; Rz=Rz0; Vx=Vx0; Vy=Vy0; Vz=Vz0; m=m0; r=r0; xBreak=x0Break; yBreak=y0Break; zBreak=z0Break;}
void Body::Move_x(double dt, double Constante){Rx+=Vx*(Constante*dt);}
void Body::Move_y(double dt, double Constante){Ry+=Vy*(Constante*dt);}
void Body::Move_z(double dt, double Constante){Rz+=Vz*(Constante*dt);}
void Body::Move_Vx(double dt, double Constante){Vx+=Fx*(Constante*dt/m);}
void Body::Move_Vy(double dt, double Constante){Vy+=Fy*(Constante*dt/m);}
void Body::Move_Vz(double dt, double Constante){Vz+=Fz*(Constante*dt/m);}
void Body::EraseForce(void){Fx=0;Fy=0;Fz=0;}
void Body::AddForce_x(double F0x){Fx+=F0x;}
void Body::AddForce_y(double F0y){Fy+=F0y;}
void Body::AddForce_z(double F0z){Fz+=F0z;}
void Body::Draw(void){cout<<", "<<xBreak+Rx<<"+"<<r<<"*cos(t),"<<yBreak+Ry<<"+"<<r<<"*sin(t)";}
void Body::Draw3D(void){cout<<", "<<xBreak+Rx<<"+"<<r<<"*sin(v)*cos(u),"<<yBreak+Ry<<"+"<<r<<"*sin(v)*sin(u),"<<zBreak+Rz<<"+"<<r<<"*cos(v)";}
/*--------------------------------------CLASS COLLIDER-------------------------------------------------------------------------------*/
class Collider{
private:
  
public:
  void AllForces(Body* Spring);
  void InteractionForce(Body & Spring1, Body & Spring2);
};
/*--------------------------------------COLLIDER'S FUNCTIONS------------------------------------------------------------------------*/
void Collider::AllForces(Body* Spring){
  int i,j,k;
  /*Erase all Forces and  Torsions*/
  for(i=0;i<N3D;i++){Spring[i].EraseForce();}
  /*Add Gravitational Force due by the Earth's Gravitational Field*/
  for(i=0;i<N3D;i++){Spring[i].AddForce_z(-Spring[i].m*g);}
  /*Add forces on x for extremum springs in x axis*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){for(k=0;k<Nz;k++){
	Spring[0+j*Ny*Nz+k*Nz].AddForce_x(-Kx*Spring[0+j*Ny*Nz+k*Nz].Rx);
	Spring[Nx-1+j*Ny*Nz+k*Nz].AddForce_x(-Kx*Spring[Nx-1+j*Ny*Nz+k*Nz].Rx);}}}
  /*Add forces on y for extremum springs in y axis*/
  for(i=0;i<Nx;i++){for(k=0;k<Nz;k++){
      Spring[i*Nx*Nz+k].AddForce_y(-Ky*Spring[i*Nx*Nz+k].Ry);
      Spring[Nx*(Ny-1)+i*Nx*Ny+k].AddForce_y(-Ky*Spring[Nx*(Ny-1)+i*Nx*Ny+k].Ry);}}
  /*Add forces on z for extremum springs in z axis*/
  for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
      Spring[0+i+j*Ny].AddForce_z(-Kz*Spring[0+i+j*Ny].Rz);
      Spring[Nx*Ny*(Nz-1)+i+j*Ny].AddForce_z(-Kz*Spring[Nx*Ny*(Nz-1)+i+j*Ny].Rz);}}
  /*Local Interaction*/
  for(i=0;i<N2D;i++){for(k=0;k<Nz;k++){
      if(k!=Nz-1){
	if(i%Nx==Nx-1){
	  Spring[i+k*Nx*Ny].AddForce_x(0);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx+k*Nx*Ny]);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx*Ny+k*Nx*Ny]);}/*Last y-lines*/
	else if((N2D-i)<=Ny){
	  Spring[i+k*Nx*Ny].AddForce_y(0);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+1+k*Nx*Ny]);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx*Ny+k*Nx*Ny]);}/*Last x-lines*/
	else if(i==N2D-1){
	  Spring[i+k*Nx*Ny].AddForce_x(0);
	  Spring[i+k*Nx*Ny].AddForce_y(0);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx*Ny+k*Nx*Ny]);}/*Last Spring xy Corners*/
	else{
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+1+k*Nx*Ny]);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx+k*Nx*Ny]);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx*Ny+k*Nx*Ny]);}}
      /*Last Z plane of Springs*/
      else if(k==Nz-1){
	if(i%Nx==Nx-1){
	  Spring[i+k*Nx*Ny].AddForce_x(0);
	  Spring[i+k*Nx*Ny].AddForce_z(0);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx+k*Nx*Ny]);}/*Last y-line*/
	else if((N2D-i)<=Ny){
	  Spring[i+k*Nx*Ny].AddForce_y(0);
	  Spring[i+k*Nx*Ny].AddForce_z(0);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+1+k*Nx*Ny]);}/*Last x-line*/
	else if(i==N2D-1){
	  Spring[i+k*Nx*Ny].AddForce_x(0);
	  Spring[i+k*Nx*Ny].AddForce_y(0);
	  Spring[i+k*Nx*Ny].AddForce_z(0);}/*Last Spring*/
	else{
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+1+k*Nx*Ny]);
	  InteractionForce(Spring[i+k*Nx*Ny],Spring[i+Nx+k*Nx*Ny]);
	  Spring[i+k*Nx*Ny].AddForce_z(0);}}}}}

void Collider::InteractionForce(Body & Spring1, Body & Spring2){
  double drx=Spring2.Rx-Spring1.Rx;
  double dry=Spring2.Ry-Spring1.Ry;
  double drz=Spring2.Rz-Spring1.Rz;
  double F1x,F1y,F1z;
  F1x=Kx*drx;F1y=Ky*dry;F1z=Kz*drz;
  Spring1.AddForce_x(F1x);
  Spring1.AddForce_y(F1y);
  Spring1.AddForce_z(F1z);
  Spring2.AddForce_x(-F1x);
  Spring2.AddForce_y(-F1y);
  Spring2.AddForce_z(-F1z);
}
/*-------------------------------------CLASS OSCILLATOR--------------------------------------------------*/
class Oscillator{
private:
  
public:
  void DrawInteractionSpring(Body & Spring1,Body & Spring2);
};
/*-------------------------------------OSCILLATOR'S FUNCTIONS--------------------------------------------*/
void Oscillator::DrawInteractionSpring(Body & Spring1,Body & Spring2){//Falta cambiar
  vector3D D,X; double d,a;
  D.cargue(Spring2.Rx-Spring1.Rx,Spring2.Ry-Spring1.Ry,Spring2.Rz-Spring1.Rz); d=norma(D);
  X.cargue(1,0,0);
  a=(D*X)/d; /*Angle between D and x-axis*/
  /*Cicloid Parametrization
  x=r*(t-sin(t))
  y=r*(1-cos(t))
  Active Rotation Matrix [cos -sin; sin cos]*/
  cout<<", "<<Spring1.Rx<<"+"<<d<<"*(t-sin(4*t))*"<<cos(a)<<"-"<<d<<"*(1-cos(4*t))*"<<sin(a)
      <<"," <<Spring1.Ry<<"+"<<d<<"*(t-sin(4*t))*"<<sin(a)<<"+"<<d<<"*(1-cos(4*t))*"<<cos(a);
}
/*-------------------------------------GLOBAL'S FUNCTIONS------------------------------------------------*/
/*-------------------------------------3D----------------------------------------------------------------*/
void StartAnimation3D(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'MySpring3D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xlabel 'x'"<<endl;
  cout<<"set ylabel 'y'"<<endl;
  cout<<"set zlabel 'z'"<<endl;
  cout<<"set xrange [-150:350]"<<endl;
  cout<<"set yrange [-200:500]"<<endl;
  cout<<"set zrange [-400:160]"<<endl; 
  //cout<<"set autoscale"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set urange [0:2*pi]"<<endl;
  cout<<"set vrange [0:pi]"<<endl;
  cout<<"set isosamples 12,12"<<endl;
  //cout<<"set view 90, 90, 0.5, 3"<<endl;  
}
void StartSquare3D(void){cout<<"splot 0,0,0 ";}
void StartFixedMasses3D(double x0Break,double y0Break,double z0Break,double r0){
  int i,j;
  //for(i=0;i<Nx;i++){for(j=0;j<Nz;j++){
  //    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*Front Plane*/
  //for(i=0;i<Nx;i++){for(j=0;j<Nz;j++){
  //    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*BackPlane*/
  //for(i=0;i<Ny;i++){for(j=0;j<Nz;j++){
  //    cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*Left Plane*/
  //for(i=0;i<Ny;i++){for(j=0;j<Nz;j++){
  //    cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*RigPlane*/
  //for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
  //    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*j<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}}/*Down Plane*/
  //for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
  //    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*j<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}}/*TopPlane*/
  /*Lines Corners*/
  //for(i=0;i<Nx;i++){
  //  cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Front*/
  //for(i=0;i<Nx;i++){
  //  cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Back*/
  //for(i=0;i<Ny;i++){
  //  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Left*/
  //for(i=0;i<Ny;i++){
  //  cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Right*/
  for(i=0;i<Nx;i++){
    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Front*/
  for(i=0;i<Nx;i++){
    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Back*/
  for(i=0;i<Ny;i++){
    cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Left*/
  for(i=0;i<Ny;i++){
    cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Right*/
  //for(i=0;i<Nz;i++){
  //  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Front Left*/
  //for(i=0;i<Nz;i++){
  //  cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Front Right*/
  //for(i=0;i<Nz;i++){
  //  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Back Left*/
  //for(i=0;i<Nz;i++){
  //  cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Back Right*/
  /*Corners*/
  //cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower Left Front*/
  //cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower Right Front*/
  //cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower Left Back*/
  //cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower RightBack*/
  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top Left Front*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top Right Front*/
  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top Left Back*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top RightBack*/
}

void FinishSquare(void){cout<<endl;}
/*-------------------------------------MAIN PROGRAM------------------------------------------------------*/

int main(void)
{
  int i,j,k,M=Ny*(Nx-1)+Nx*(Ny-1);
  double t,dt=1e-1;          /*Times*/
  double tdrawings,Ndrawings;/*Accountants*/
  Body Spring[N3D];          /*Oscillators Masses*/
  Oscillator SpringK[M];     /*Oscillators Springs*/
  Collider Hooke;            /*Collider*/
  
  double m0=1,r0=5,Rx0=Lx/8,Ry0=Lx/10,Rz0=Lx/5,Vx0=0,Vy0=0,Vz0=0; /*Initial Conditions*/
  double x0Break=Lx, y0Break=3*Ly/2, z0Break=Lz/2;                /*Distance of Separation in each axis*/
  double omegax=sqrt(Kx/m0),omegay=sqrt(Ky/m0),omegaz=sqrt(Kz/m0);/*One of Frequencies in each axis*/
  double Tx=2*M_PI/omegax, tmax=1*Tx;                             /*Time Step and Tmax*/

  /*Start Springs*/
  for(k=0;k<Nz;k++){for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
	Spring[i+Nx*j+N2D*k].Start(Rx0,Ry0,Rz0,Vx0,Vy0,Vz0,m0,r0,i*x0Break,j*y0Break,k*z0Break);}}}
  /*Start Animation*/
  StartAnimation3D();  Ndrawings=100;
  /*Integration*/
  for(t=tdrawings=0;t<tmax;t+=dt,tdrawings+=dt)
    {
      
      if(tdrawings>tmax/Ndrawings){
	StartSquare3D();
	StartFixedMasses3D(x0Break,y0Break,z0Break,r0);
	for(i=0;i<N3D;i++){Spring[i].Draw3D();}
	FinishSquare();
	tdrawings=0;}
      /*
      for(i=0;i<N3D;i++){cout<<t<<"\t"<<i<<"\t"<<Spring[i].Getx()<<"\t"<<Spring[i].Gety()<<"\t"<<Spring[i].Getz()<<endl;}
      */
      //Move with Omelyan PEFRL
      for(i=0;i<N3D;i++){Spring[i].Move_x(dt,Xi);}
      for(i=0;i<N3D;i++){Spring[i].Move_y(dt,Xi);}
      for(i=0;i<N3D;i++){Spring[i].Move_z(dt,Xi);}
      Hooke.AllForces(Spring);
      for(i=0;i<N3D;i++){Spring[i].Move_Vx(dt,(1-2*Lambda)/2);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vy(dt,(1-2*Lambda)/2);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vz(dt,(1-2*Lambda)/2);}
      for(i=0;i<N3D;i++){Spring[i].Move_x(dt,Chi);}
      for(i=0;i<N3D;i++){Spring[i].Move_y(dt,Chi);}
      for(i=0;i<N3D;i++){Spring[i].Move_z(dt,Chi);}
      Hooke.AllForces(Spring);
      for(i=0;i<N3D;i++){Spring[i].Move_Vx(dt,Lambda);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vy(dt,Lambda);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vz(dt,Lambda);}
      for(i=0;i<N3D;i++){Spring[i].Move_x(dt,1-2*(Chi+Xi));}
      for(i=0;i<N3D;i++){Spring[i].Move_y(dt,1-2*(Chi+Xi));}
      for(i=0;i<N3D;i++){Spring[i].Move_z(dt,1-2*(Chi+Xi));}
      Hooke.AllForces(Spring);
      for(i=0;i<N3D;i++){Spring[i].Move_Vx(dt,Lambda);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vy(dt,Lambda);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vz(dt,Lambda);}
      for(i=0;i<N3D;i++){Spring[i].Move_x(dt,Chi);}
      for(i=0;i<N3D;i++){Spring[i].Move_y(dt,Chi);}
      for(i=0;i<N3D;i++){Spring[i].Move_z(dt,Chi);}
      Hooke.AllForces(Spring);
      for(i=0;i<N3D;i++){Spring[i].Move_Vx(dt,(1-2*Lambda)/2);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vy(dt,(1-2*Lambda)/2);}
      for(i=0;i<N3D;i++){Spring[i].Move_Vz(dt,(1-2*Lambda)/2);}
      for(i=0;i<N3D;i++){Spring[i].Move_x(dt,Xi);}
      for(i=0;i<N3D;i++){Spring[i].Move_y(dt,Xi);}
      for(i=0;i<N3D;i++){Spring[i].Move_z(dt,Xi);}
    }
  return 0;
}
