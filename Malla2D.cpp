#include<iostream>
#include<fstream>
#include<cmath>
#include"Vector.h"
#include"Random64.h"
using namespace std;

const double g=1;
const double Kx=1,Ky=1,Kz=1;
const double Lx=100,Ly=100,Lz=100;
const int Nx=3,Ny=3,Nz=1;
const int N2D=Nx*Ny,N3D=Nx*Ny*Nz;

const double Xi=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Chi=-0.06626458266981849;

class Body;
class Collider;
class Oscillator;
/*--------------------------------------------------CLASS BODY------------------------------------------------------------------------*/
class Body{
private:
  double Rx,Ry,Rz,Vx,Vy,Vz,Fx,Fy,Fz,m,r,xcorrido,ycorrido,zcorrido; 
public:
  void Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
	     double m0,double r0,double x0corrido,double y0corrido,double z0corrido);
  void Move_x(double dt, double Constante);
  void Move_y(double dt, double Constante);
  void Move_z(double dt, double Constante);
  void Move_Vx(double dt, double Constante);
  void Move_Vy(double dt, double Constante);
  void Move_Vz(double dt, double Constante);
  double Getx(void){return Rx;};
  double Gety(void){return Ry;};
  double Getz(void){return Rz;};
  void EraseForce(void);
  void AddForce_x(double F0x);
  void AddForce_y(double F0y);
  void AddForce_z(double F0z);
  void Draw(void);
  friend class Collider;
  friend class Oscillator;
};
/*-------------------------------------------------BODY'S FUNCTIONS-------------------------------------------------------------------*/
void Body::Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,double m0,double r0,double x0corrido,double y0corrido,double z0corrido)
{Rx=Rx0; Ry=Ry0; Rz=Rz0; Vx=Vx0; Vy=Vy0; Vz=Vz0; m=m0; r=r0; xcorrido=x0corrido; ycorrido=y0corrido; zcorrido=z0corrido;}
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
void Body::Draw(void){cout<<", "<<xcorrido+Rx<<"+"<<r<<"*cos(t),"<<ycorrido+Ry<<"+"<<r<<"*sin(t)";}
/*--------------------------------------CLASS COLLIDER-------------------------------------------------------------------------------*/
class Collider{
private:
  
public:
  void AllForces(Body* Spring);
  void InteractionForce(Body & Spring1, Body & Spring2);
};
/*--------------------------------------COLLIDER'S FUNCTIONS------------------------------------------------------------------------*/
void Collider::AllForces(Body* Spring){
  int i,j;
  for(i=0;i<N2D;i++){Spring[i].EraseForce();}
  //for(i=0;i<N2D;i++){Spring[i].AddForce_y(-Spring[i].m*g);}
  //Add forces on x for extremum springs in x axis
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
      Spring[0+j*Ny].AddForce_x(-Kx*Spring[0+j*Ny].Rx); Spring[Nx-1+j*Ny].AddForce_x(-Kx*Spring[Nx-1+j*Ny].Rx);}}
  //Add forces on y for extremum springs in y axis
  for(i=0;i<Nx;i++){Spring[i].AddForce_y(-Ky*Spring[i].Ry); Spring[Nx*(Ny-1)+i].AddForce_y(-Ky*Spring[Nx*(Ny-1)+i].Ry);}
  //Interaction local
  for(i=0;i<N2D;i++){
    if(i%Nx==Nx-1){Spring[i].AddForce_x(0);InteractionForce(Spring[i],Spring[i+Nx]);}      //Last y-line
    else if((N2D-i)<=Ny){Spring[i].AddForce_y(0);InteractionForce(Spring[i],Spring[i+1]);} //Last x-line
    else if(i==N2D-1){Spring[i].AddForce_x(0);Spring[i].AddForce_y(0);}                    //Last Spring
    else{InteractionForce(Spring[i],Spring[i+1]);InteractionForce(Spring[i],Spring[i+Nx]);}}
}

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
//class Oscillator{
//private:

//public:
//void HookeLaw(Body & Spring1,double K);
//};
/*-------------------------------------OSCILLATOR'S FUNCTIONS--------------------------------------------*/

/*-------------------------------------GLOBAL'S FUNCTIONS------------------------------------------------*/
void StartAnimation(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'MiResorte2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  //cout<<"set xrange [-20:400]"<<endl;
  //cout<<"set yrange [-20:20]"<<endl;
  cout<<"set autoscale"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:2*pi]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void StartSquare(void){cout<<"plot 0,0 ";}
void FinishSquare(void){cout<<endl;}
/*-------------------------------------MAIN PROGRAM------------------------------------------------------*/

int main(void)
{
  double t,dt=1e-1;
  double tdrawings,Ndrawings;
  Body Spring[N2D];
  Collider Hooke;
  int i,j;
  
  double m0=1,r0=5;
  double Rx0=Lx/8,Ry0=Lx/8,Rz0=0,Vx0=0,Vy0=0,Vz0=0;
  double x0corrido=50, y0corrido=50, z0corrido=0;
  double omegax=sqrt(Kx/m0);
  double Tx=2*M_PI/omegax, tmax=5*Tx;
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
      Spring[i+Nx*j].Start( Rx0, Ry0, Rz0, Vx0, Vy0, Vz0, m0, r0,i*x0corrido,j*y0corrido,z0corrido);}}
  
  StartAnimation();  Ndrawings=1000;
  
  for(t=tdrawings=0;t<tmax;t+=dt,tdrawings+=dt)
    {
      
      if(tdrawings>tmax/Ndrawings){
	StartSquare();
	for(i=0;i<N2D;i++){Spring[i].Draw();}
	FinishSquare();
	tdrawings=0;}
      /*
      cout<<t<<"\t"<<Spring[0].Getx()<<endl;
      */
      //Movese con Omelyan PEFRL
      for(i=0;i<N2D;i++){Spring[i].Move_x(dt,Xi);}
      for(i=0;i<N2D;i++){Spring[i].Move_y(dt,Xi);}
      for(i=0;i<N2D;i++){Spring[i].Move_z(dt,Xi);}
      Hooke.AllForces(Spring);
      for(i=0;i<N2D;i++){Spring[i].Move_Vx(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vy(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vz(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D;i++){Spring[i].Move_x(dt,Chi);}
      for(i=0;i<N2D;i++){Spring[i].Move_y(dt,Chi);}
      for(i=0;i<N2D;i++){Spring[i].Move_z(dt,Chi);}
      Hooke.AllForces(Spring);
      for(i=0;i<N2D;i++){Spring[i].Move_Vx(dt,Lambda);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vy(dt,Lambda);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vz(dt,Lambda);}
      for(i=0;i<N2D;i++){Spring[i].Move_x(dt,1-2*(Chi+Xi));}
      for(i=0;i<N2D;i++){Spring[i].Move_y(dt,1-2*(Chi+Xi));}
      for(i=0;i<N2D;i++){Spring[i].Move_z(dt,1-2*(Chi+Xi));}
      Hooke.AllForces(Spring);
      for(i=0;i<N2D;i++){Spring[i].Move_Vx(dt,Lambda);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vy(dt,Lambda);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vz(dt,Lambda);}
      for(i=0;i<N2D;i++){Spring[i].Move_x(dt,Chi);}
      for(i=0;i<N2D;i++){Spring[i].Move_y(dt,Chi);}
      for(i=0;i<N2D;i++){Spring[i].Move_z(dt,Chi);}
      Hooke.AllForces(Spring);
      for(i=0;i<N2D;i++){Spring[i].Move_Vx(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vy(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D;i++){Spring[i].Move_Vz(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D;i++){Spring[i].Move_x(dt,Xi);}
      for(i=0;i<N2D;i++){Spring[i].Move_y(dt,Xi);}
      for(i=0;i<N2D;i++){Spring[i].Move_z(dt,Xi);}
    }
  return 0;
}
