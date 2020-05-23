#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_poly.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>

#define id 11
#define id1 6
#define M 1
#define N (M+1)*id
#define nhost 3
#define niter 100000
#define nt 158
#define atl 23
#define npart 10000

class sys_para{
public:
    double ep,a,a0,c,d,tpi,w,NN,T,K,x1,x2,bx,by,xm,ym,xmp,ymp,xmp1,ymp1,ri,rt2,rmax,ang,Prr1,Prr2,ri2,rb,rb2,SI,NH10,NHic;
double PIH1,PIH2,PIH3,PIS;// birth rates from host 1 to 3 and sand fly
double NS,NH,N2,N3,Ntot,NNav,Nr1,Nr2,Nr3,NSr;//
double muH1,muH2,muH3,muS;//death rates form host 1 to 3 and sand fly
double sig1,rH,d1,gam1;//sig1 is from E1 to I1,rp is rate of reporting the infectious,d1 is disease induced mortality,gam1 is the rate of tranfer from treated class to PKDL
double sig2,sig3;//sig2 is rate of transfer from E2 to I2,sig3 is rate of transfer from E3 to I3
double e1,tau1,taup1,gam2;//e1 is rate of recovery directly from latent state, tau1 from infectious state, taup1 from reported state,gam2 from PKDL state of host 1
double del1,del2,del3;//rate of loss of immunity from host 1 to 3
double lamH1,lamH1j,lamH2,lamH2j,lamH3,lamH3j,lamS,lamSj;//Force of infection for hosts 1 to 3 and sand fly;
double Q1,Q2,Q3,Q4,Q5;//proportion of host of populations 1 to 3;
double r11,r12,r13,r14,r21,r22,r23,r24;//
double rho1,rho2,th1,th2;
double tau2,tau3,d2,d3;
double R02,R0,R0d;
double C1,C2,C3;
double b,bet1,bet2,bet3,betS;
double br1,br20,br21,br22,br23,br24,br3,br40,br41,br42,br5,br60,br61,br62;
double sigv;
double F2,F1,A,B,C,lH1,lH2,lV1,lV2,NH1,NH2,Ih1,Ih2,Sh2,Sh1,dlH2,f,g,eta,Sv1,Sv2,Svv1,Svv2;
double Sh10,Ih10,Sh20,Ih20,Sh30,Ih30,Sv0,Iv0;
double Eh1,Eh2,Rh1,Rh2,Ev1,Ev2,Iv1,Iv2;
double QQ1,QQ2,QQ3,QQ4,QQ5,R0T,PIST,p1,p2,rK,q1,q2,Rc0,Ic,yesc;
double NH1j,N2j,N3j;

double llh1[id1+1],llh2[id1+1],llh3[id1+1],llv[id1+1],EEv[id1+1],IIv[id1+1],Imllh1[id1+1],Imllh2[id1+1],Imllh3[id1+1],NNH1[id1+1],NNH2[id1+1],NNH3[id1+1],IIh1[id1+1],IIh2[id1+1],IIh3[id1+1];
double ImNNH1[id1+1],ImNNH2[id1+1],ImNNH3[id1+1],ImIIh1[id1+1],ImIIh2[id1+1],ImIIh3[id1+1],SSh1[id1+1],SSh2[id1+1],SSh3[id1+1],ImSSh1[id1+1],ImSSh2[id1+1],ImSSh3[id1+1],ddlH2;
double EEh1[id1+1],RRh1[id1+1],SSv[id1+1];
double m0,m[nhost+1],FF2[nhost+1],FF1[nhost+1],N0[nhost+1],N0r[nhost+1],Nr[nhost+1],alph[nhost+1],Qv,NN0av,NS0r,NN0rav;
double ra,rc,re;
double Q01,Q02,Q03,Q04,R0M;
double ep1,ep2;
double M1,M2,M3,G1,G2,G3;
double Bv,L0,L1,L2,O0,O1,O2,O3,LHSf,D0,D1,D2,D3,D4,D5,D6,MG1,MG2,MG3,MG4,MG5,MG6;
double A0,A1,A2,A3,A4,A5,A6;
double Norm[id1+1][2];

};
void gs(double *y,double *sum);
void rootmulti(sys_para &sp);
double DeterMmultihost(sys_para &sp,double A);

void rootmultiend(sys_para &sp);
double DeterMmultihostend(sys_para &sp,double A);

void bubbleSort(double norm[][3],int nn);

void bubbleSortpoly(double norm[][2]);//this sorts the polynomial roots
void highpowers(sys_para &sp);

int
func (double t, const double y0[], double dy[],
void *params)
{
//double mu = *(double *)params;
int i,j;
sys_para sp=*(sys_para*)params;
  
sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
sp.NNav=sp.alph[1]*sp.NH+sp.alph[2]*sp.N2+sp.alph[3]*sp.N3;

sp.lamS=sp.b*sp.betS*(sp.alph[1]*sp.C1*(y0[3]+sp.p1*y0[2])+sp.alph[2]*sp.C2*y0[6]+sp.alph[3]*sp.C3*y0[8])/sp.NNav;

sp.lamH1=sp.b*sp.alph[1]*sp.bet1*y0[11]/sp.NNav;
sp.lamH2=sp.b*sp.alph[2]*sp.bet2*y0[11]/sp.NNav;
sp.lamH3=sp.b*sp.alph[3]*sp.bet3*y0[11]/sp.NNav;
//----------Human host section with SE(I-P)(PKDL)RS
  dy[1]=sp.PIH1-sp.lamH1*y0[1]+sp.del1*y0[4]-sp.muH1*y0[1];//susceptible
  dy[2]=sp.lamH1*y0[1]-(sp.sig1+sp.muH1)*y0[2];//Latent
  dy[3]=sp.f*sp.sig1*y0[2]-sp.Q2*y0[3];//infectious
  dy[4]=sp.sig1*(1-sp.f)*y0[2]+sp.tau1*y0[3]-sp.Q3*y0[4];//Recovered

//Second host
  dy[5]=sp.PIH2-sp.lamH2*y0[5]+sp.del2*y0[6]-sp.muH2*y0[5];
  dy[6]=sp.lamH2*y0[5]-sp.Q4*y0[6];
//Third host
  dy[7]=sp.PIH3-sp.lamH3*y0[7]+sp.del3*y0[8]-sp.muH3*y0[7];
  dy[8]=sp.lamH3*y0[7]-sp.Q5*y0[8];
//---Sand fly Section----//
  dy[9]=sp.PIS-sp.lamS*y0[9]-sp.muS*y0[9];//Susceptible
  dy[10]=sp.lamS*y0[9]-sp.Qv*y0[10];//Latent
  dy[11]=sp.sigv*y0[10]-sp.muS*y0[11];//Infectious

  for(int itr=1;itr<=M;itr++){
    j=itr*id;
sp.NH1j=y0[1+j]+y0[2+j]+y0[3+j]+y0[4+j];
sp.N2j=y0[5+j]+y0[6+j];
sp.N3j=y0[7+j]+y0[8+j];

sp.lamH1j=sp.b*sp.bet1*sp.alph[1]*y0[11+j]/sp.NNav-sp.b*sp.bet1*sp.alph[1]*y0[11]*(sp.alph[1]*sp.NH1j+sp.alph[2]*sp.N2j+sp.alph[3]*sp.N3j)/sp.NNav;
sp.lamH2j=sp.b*sp.bet2*sp.alph[2]*y0[11+j]/sp.NNav-sp.b*sp.bet2*sp.alph[2]*y0[11]*(sp.alph[1]*sp.NH1j+sp.alph[2]*sp.N2j+sp.alph[3]*sp.N3j)/sp.NNav;
sp.lamH3j=sp.b*sp.bet3*sp.alph[3]*y0[11+j]/sp.NNav-sp.b*sp.bet3*sp.alph[3]*y0[11]*(sp.alph[1]*sp.NH1j+sp.alph[2]*sp.N2j+sp.alph[3]*sp.N3j)/sp.NNav;

sp.lamSj=sp.b*sp.betS*(sp.alph[1]*sp.C1*(y0[3+j]+sp.p1*y0[2+j])+sp.alph[2]*sp.C2*y0[6+j]+sp.alph[3]*sp.C3*y0[8+j])/sp.NNav-sp.b*sp.betS*(sp.alph[1]*sp.C1*(y0[3]+sp.p1*y0[2])+sp.alph[2]*sp.C2*y0[6]+sp.alph[3]*sp.C3*y0[8])*(sp.alph[1]*sp.NH1j+sp.alph[2]*sp.N2j+sp.alph[3]*sp.N3j)/sp.NNav;

  dy[1+j]=-sp.lamH1j*y0[1]-sp.lamH1*y0[1+j]+sp.del1*y0[4+j]-sp.muH1*y0[1+j];//susceptible
  dy[2+j]=sp.lamH1j*y0[1]+sp.lamH1*y0[1+j]-sp.Q1*y0[2+j];//Latent
  dy[3+j]=sp.f*sp.sig1*y0[2+j]-sp.Q2*y0[3+j];//infectious
  dy[4+j]=sp.sig1*(1-sp.f)*y0[2+j]+sp.tau1*y0[3+j]-sp.Q3*y0[4+j];//Recovered
//Second host
  dy[5+j]=-sp.lamH2j*y0[5]-sp.lamH2*y0[5+j]+sp.del2*y0[6+j]-sp.muH2*y0[5+j];
  dy[6+j]=sp.lamH2j*y0[5]+sp.lamH2*y0[5+j]-sp.Q4*y0[6+j];
//Third host
  dy[7+j]=-sp.lamH3j*y0[7]-sp.lamH3*y0[7+j]+sp.del3*y0[8+j]-sp.muH3*y0[7+j];
  dy[8+j]=sp.lamH3j*y0[7]+sp.lamH3*y0[7+j]-sp.Q5*y0[8+j];
//---Sand fly Section----//
  dy[9+j]=-sp.lamSj*y0[9]-sp.lamS*y0[9+j]-sp.muS*y0[9+j];//Susceptible
  dy[10+j]=sp.lamSj*y0[9]+sp.lamS*y0[9+j]-sp.Qv*y0[10+j];//Latent
  dy[11+j]=sp.sigv*y0[10+j]-sp.muS*y0[11+j];//Infectious

  }


 
return GSL_SUCCESS;
}
int
jac (double t, const double y0[], double *dfdy,
double dfdt[], void *params)
{
/*
//double mu = *(double *)params;
sys_para sp=*(sys_para*)params;
gsl_matrix_view dfdy_mat
= gsl_matrix_view_array (dfdy, 2, 2);
gsl_matrix * m = &dfdy_mat.matrix;

gsl_matrix_set (m, 0, 0, -sp.a*y[1]-sp.b);
gsl_matrix_set (m, 0, 1, -sp.a*y[0]+sp.c);

gsl_matrix_set (m, 1, 0, sp.a*y[1]);
gsl_matrix_set (m, 1, 1, sp.a*y[0]-sp.b-sp.c);

dfdt[0] = 0.0;
dfdt[1] = 0.0;*/
return NULL;//GSL_SUCCESS;
}
int
main (void)
{
FILE *fp=fopen("SImultihost_twodeadendhosts.dat","w");
int i,j,nt1,atl1,it=0;;
const gsl_odeiv_step_type * T
= gsl_odeiv_step_rkck;
gsl_odeiv_step * s
= gsl_odeiv_step_alloc (T, N);
gsl_odeiv_control * c
= gsl_odeiv_control_y_new (1e-5, 0.0);
gsl_odeiv_evolve * e
= gsl_odeiv_evolve_alloc (N);
//double mu = 10;
sys_para sp;

 sp.NN=1;
sp.tpi=8*atan(1.0); 




gsl_odeiv_system sys = {func, jac, N, &sp};
int l;
double t = 0.0, t1 = 100.0;
double h = 1e-4;
double y0[N+1] = {0,0.1, 0.2, 0.3, 0.1 },sum[M+1];
double d_ep,epmax,z_min,z_max,zz_max,zz[niter+1][id],yy[atl+1][3],Pr[nt+1][3],Pr2[nt+1][3],Pr3[nt+1][3],s1,s2,ti;
double cf[214],z[214],*x1,*x0,nrm;;

//gsl_vector *norm=gsl_vector_alloc(id1+1);
gsl_permutation *p=gsl_permutation_alloc(id1+1);


sp.muH1=1/(68.1*365);//human death rate//https://www.disabled-world.com/calculators-charts/in-lifespan.php
sp.muS=1/14.0;//Sandfly death rate//Shyam Sundar paper and http://howmed.net/community-medicine/sandfly-characteristics-life-cycle-and-control-measures/

//Brith rates
//Population of Muzaffarpur District in Bihar in 2011 was 4801062. Therefore with the life span of 68.1 years in bihar rural we get birth rate of 193.151167703 births per day
//Out of the total Muzaffarpur population for 2011 census, 9.86 percent lives in urban regions of district. As per 2011 census, 90.14 % population of Muzaffarpur districts lives in rural areas of villages.
//https://www.census2011.co.in/census/district/68-muzaffarpur.html
sp.PIH1=193.151167703;//humans birth rate

sp.PIS=100;//sand flies birth rate

sp.sig1=0.006;//intrinsic incubation period (i.e. for humans) 
sp.sigv=1/5.0;//extrinsic incubation period(i.e. for sandflies) Shyam Sundar

sp.d1=0.009;//0.000280;//0.009;////disease induced death rate in humans//

//recovery from I,P stages//
sp.gam1=0.001;//only 12.3 % percent of reporting happend in the indian state of bihar in goverment hospitals and the rest  happened at private clinics and NGOs
          //https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-3156.2006.01647.x
//We are taking 1 % reporting of infectious cases per day i.e. sp.gam1=0.01

sp.f=0.14;//Approximately one in seven aysmptomatics proceeds to have VL

sp.tau1=0.00476;//From infectious stage
sp.rH=1/30.0;//recovery from treatment (one month treatment)



//Transmission rates
sp.bet1=0.74;//Sandflies to humans
sp.del1=0.000913;//Rate of loss of immunity humans
sp.C1=1;
sp.alph[1]=1;
sp.betS=0.158;//Host to sandflies when no host competence is involved
sp.p1=0.01;//weight of infection in sand flies coming from latent humans

//Sandflies biting rate
sp.b=0.25;// vector biting rate

//host 2 parameters (dogs)
sp.bet2=0.74;
sp.muH2=0.000167;
sp.PIH2=sp.muH2*500000;
sp.alph[2]=1;
sp.C2=1;
sp.d2=0;
sp.del2=sp.del1;



//host 3 parameters
sp.bet3=0.74;
sp.muH3=0.01;
sp.PIH3=60000*sp.muH3;
sp.alph[3]=1;
sp.C3=1;
sp.d3=0;
sp.del3=sp.del1;

sp.ep1=1.0000;
sp.ep2=1;

sp.bet2=sp.ep1*sp.bet2;
sp.bet3=sp.ep2*sp.bet3;

sp.R0=0.35;//critical R0M= 0.982670 when R0= 1.09729 for dead end hosts C2=C3=0 and alph[1]=alph[2]=alph[3]=1
               //ciritcal R0M=0.765256 when R0=0.267493 for tow competeive hosts C2=C3=1 and crosses R0M=1 when R0=0.35

//declaring the attractor points here
//for(sp.R0=0.2;sp.R0<=1.13;sp.R0+=0.00001)
{
sp.Q01=sp.sig1+sp.muH1;
sp.Q02=sp.d1+sp.tau1+sp.muH1;
sp.Q03=sp.del1+sp.muH1;

sp.Q1=sp.sig1+sp.muH1;
sp.Q2=sp.d1+sp.tau1+sp.muH1;
sp.Q3=sp.del1+sp.muH1;

sp.Q4=sp.d2+sp.del2+sp.muH2;
sp.Q5=sp.d3+sp.del3+sp.muH2;

sp.Qv=sp.sigv+sp.muS;

sp.FF1[1]=((1-sp.f)*sp.sig1*sp.Q2+sp.tau1*sp.f*sp.sig1)/(sp.Q1*sp.Q2*sp.Q3);
sp.FF1[2]=1/sp.Q4;
sp.FF1[3]=1/sp.Q5;

sp.FF2[1]=1/sp.Q1+sp.f*sp.sig1/(sp.Q1*sp.Q2)+sp.FF1[1];
sp.FF2[2]=sp.FF1[2];
sp.FF2[3]=sp.FF1[3];

sp.m0=(sp.f*sp.sig1+sp.p1*sp.Q02)/(sp.Q01*sp.Q02);
sp.m[1]=sp.alph[1]*sp.C1*(sp.f*sp.sig1+sp.p1*sp.Q2)/(sp.Q1*sp.Q2);
sp.m[2]=sp.alph[2]*sp.C2/sp.Q4;
sp.m[3]=sp.alph[3]*sp.C3/sp.Q5;

sp.ra=sp.alph[1]*(1-sp.del1*sp.FF1[1])/(sp.alph[1]*sp.muH1);
sp.rc=sp.alph[2]*(1-sp.del2*sp.FF1[2])/(sp.alph[1]*sp.muH2);
sp.re=sp.alph[3]*(1-sp.del3*sp.FF1[3])/(sp.alph[1]*sp.muH3);

sp.N0[1]=sp.PIH1/sp.muH1;
sp.N0[2]=sp.PIH2/sp.muH2;
sp.N0[3]=sp.PIH3/sp.muH3;

sp.NN0av=sp.alph[1]*sp.N0[1]+sp.alph[2]*sp.N0[2]+sp.alph[3]*sp.N0[3];

sp.N0r[1]=sp.N0[1]/sp.NN0av;
sp.N0r[2]=sp.N0[2]/sp.NN0av;
sp.N0r[3]=sp.N0[3]/sp.NN0av;

sp.F1=sp.del1*((1-sp.f)*sp.sig1*sp.Q02+sp.f*sp.sig1*sp.tau1)/(sp.Q01*sp.Q02*sp.Q03);
sp.F2=(1/sp.Q01)+(sp.f*sp.sig1/(sp.Q01*sp.Q02))+(sp.F1/(sp.del1));
sp.PIS=(sp.R0*sp.R0*sp.muS*sp.muS*sp.PIH1*sp.Q01*sp.Q02*sp.Qv)/(sp.b*sp.b*sp.betS*sp.bet1*sp.sigv*(sp.f*sp.sig1+sp.p1*sp.Q02)*sp.muH1);//
sp.NS=sp.PIS/sp.muS;

sp.A=sp.F2*sp.F2+(sp.b*sp.betS*(sp.f*sp.sig1+sp.p1*sp.Q02)*sp.F2)/(sp.muS*sp.Q01*sp.Q02);
sp.B=2*sp.F2+(sp.b*sp.betS*(sp.f*sp.sig1+sp.p1*sp.Q02))/(sp.muS*sp.Q01*sp.Q02)-((sp.R0*sp.R0)*(1-sp.F1))/sp.muH1;
sp.C=1-sp.R0*sp.R0;

sp.lH1=(-sp.B+sqrt(sp.B*sp.B-4*sp.A*sp.C))/(2*sp.A);
sp.lH2=(-sp.B-sqrt(sp.B*sp.B-4*sp.A*sp.C))/(2*sp.A);
sp.Sh1=sp.PIH1/(sp.lH1*(1-sp.F1)+sp.muH1);
sp.Sh2=sp.PIH1/(sp.lH2*(1-sp.F1)+sp.muH1);
sp.NH1=(1+sp.lH1*sp.F2)*sp.Sh1;
sp.NH2=(1+sp.lH2*sp.F2)*sp.Sh2;
sp.Ih1=sp.f*sp.sig1*sp.lH1*sp.Sh1/(sp.Q01*sp.Q02);

sp.Ih2=sp.f*sp.sig1*sp.lH2*sp.Sh2/(sp.Q01*sp.Q02);


sp.NS0r=sp.NS/sp.NN0av;

sp.NN0rav=sp.alph[1]+sp.alph[2]*sp.N0[2]/sp.N0[1]+sp.alph[3]*sp.N0[3]/sp.N0[1];
sp.Bv=sp.b*sp.betS*sp.sigv*sp.NS0r/(sp.muS*sp.Qv);
//sp.R0M=sqrt(((sp.R0*sp.R0*sp.m[1])/(sp.alph[1]*sp.m0*sp.NN0rav))+sp.Bv*(sp.b*sp.alph[2]*sp.bet2*sp.m[2]*sp.N0r[2]+sp.b*sp.alph[3]*sp.bet3*sp.N0r[3]*sp.m[3]));

sp.R0M=sqrt(sp.Bv*(sp.b*sp.bet1*sp.alph[1]*sp.m[1]*sp.N0r[1]+sp.b*sp.bet2*sp.alph[2]*sp.m[2]*sp.N0r[2]+sp.b*sp.bet3*sp.alph[3]*sp.m[3]*sp.N0r[3]));
highpowers(sp);
//printf("%lf %lf %lf %lf %lf %lf\n",sp.A,sp.A2,sp.B,sp.A1,sp.C,sp.A0);
//Declaring the quintic polynomial
cf[id1]=sp.A6;//sextic term
cf[id1-1]=sp.A5;//quintic term
cf[id1-2]=sp.A4;//quartic term
cf[id1-3]=sp.A3;//cubic term
cf[id1-4]=sp.A2;//quadratic term
cf[id1-5]=sp.A1;//linear term
cf[id1-6]=sp.A0;//constant term
gsl_poly_complex_workspace *w=gsl_poly_complex_workspace_alloc(id1+1);

gsl_poly_complex_solve(cf,id1+1,w,z);
gsl_poly_complex_workspace_free(w);

for(l=0;l<id1;l++)
{nrm=z[2*l];
sp.Norm[l][0]=nrm;
sp.Norm[l][1]=l;
}

bubbleSortpoly(sp.Norm);

//Disease free equilibrium defined here

sp.Sh10=sp.PIH1/sp.muH1;
sp.Sh20=sp.PIH2/sp.muH2;
sp.Sh30=sp.PIH3/sp.muH3;
sp.Ih10=0;sp.Ih20=0;sp.Ih30=0;

sp.Sv0=sp.PIS/sp.muS;
sp.Iv0=0;
//------------------------------------------------//

//Endemic equilibria calculated here//
int sadle=0;//if sadle=1 then a sadle node bifurcation occurs-and is tested in the following loop
if(z[2*int(sp.Norm[0][1])+1]==0 || z[2*int(sp.Norm[1][1])+1]==0)
{
for(l=0;l<id1;l++)
{
int  ll=sp.Norm[l][1];
sp.llh1[l+1]=z[2*ll];
sp.llh2[l+1]=sp.alph[2]*z[2*ll]/sp.alph[1];
sp.llh3[l+1]=sp.alph[3]*z[2*ll]/sp.alph[1];


sp.SSh1[l+1]=sp.PIH1/(sp.llh1[l+1]*(1-sp.del1*sp.FF1[1])+sp.muH1);
sp.SSh2[l+1]=sp.PIH2/(sp.llh2[l+1]*(1-sp.del2*sp.FF1[2])+sp.muH2);
sp.SSh3[l+1]=sp.PIH3/(sp.llh3[l+1]*(1-sp.del3*sp.FF1[3])+sp.muH3);

sp.IIh1[l+1]=sp.f*sp.sig1*sp.llh1[l+1]*sp.SSh1[l+1]/(sp.Q1*sp.Q2);
sp.IIh2[l+1]=sp.llh2[l+1]*sp.SSh2[l+1]/(sp.Q4);
sp.IIh3[l+1]=sp.llh3[l+1]*sp.SSh3[l+1]/(sp.Q5);

sp.NNH1[l+1]=(1+sp.llh1[l+1]*sp.FF2[1])*sp.SSh1[l+1];
sp.NNH2[l+1]=(1+sp.llh2[l+1]*sp.FF2[2])*sp.SSh2[l+1];
sp.NNH3[l+1]=(1+sp.llh3[l+1]*sp.FF2[3])*sp.SSh3[l+1];

sp.llv[l+1]=sp.b*sp.betS*(sp.m[1]*sp.llh1[l+1]*sp.SSh1[l+1]+sp.m[2]*sp.llh2[l+1]*sp.SSh2[l+1]+sp.m[3]*sp.llh3[l+1]*sp.SSh3[l+1])/(sp.alph[1]*sp.NNH1[l+1]+sp.alph[2]*sp.NNH2[l+1]+sp.alph[3]*sp.NNH3[l+1]);
sp.SSv[l+1]=sp.PIS/(sp.llv[l+1]+sp.muS);
sp.EEv[l+1]=sp.llv[l+1]*sp.SSv[l+1]/sp.Qv;
sp.IIv[l+1]=sp.sigv*sp.llv[l+1]*sp.SSv[l+1]/(sp.Qv*sp.muS);

sp.EEh1[l+1]=sp.llh1[l+1]*sp.SSh1[l+1]/sp.Q1;
sp.RRh1[l+1]=sp.FF1[1]*sp.llh1[l+1]*sp.SSh1[l+1];

if(z[2*ll+1]==0)
{
sadle=1;//saddle node is born if the loop reaches this line
}
}

}

//sadle node condition and calculation of endemic equlibrium ends here

printf("done with EE %d %lf %lf\n",sadle,sp.R0,sp.R0M);

if(sadle==1)//condition over saddle node
{

//The basin boundary calculations for first attractor require to tell the centroid obtained from mean square distplacement of the attractor

 y0[1]=sp.Sh10;y0[2]=0; y0[3]=0; y0[4]=0;
 y0[5]=sp.Sh20;y0[6]=0;
 y0[7]=sp.Sh30;y0[8]=0;
 y0[9]=sp.Sv0;y0[10]=0;y0[11]=0;



sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
sp.NNav=sp.alph[1]*sp.NH+sp.alph[2]*sp.N2+sp.alph[3]*sp.N3;
    for(i=1;i<=M;i++)
       {
        sum[i]=0.0;
       }
  t=0;t1=800000;
it=0;sp.xmp=0;sp.ymp=0;
while (t < t1)
{
sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
sp.NNav=sp.alph[1]*sp.NH+sp.alph[2]*sp.N2+sp.alph[3]*sp.N3;
int status = gsl_odeiv_evolve_apply (e, c, s,
&sys,
&t, t1,
&h, y0);
gs(y0,sum);
if (status != GSL_SUCCESS)
break;
//printf ("%.5e %.5e %.5e  %.5e %.5e\n", t/365., log10(y[0]), log10(y[1]),log10(y[2]),log10(y[3]));
if(t>700000)
{
it++;
}
}
sp.xmp=fabs(y0[3]/sp.NH);//sqrt(sp.xmp/it);
//printf("%lf\n",sp.xmp);
sp.ymp=sqrt(sp.ymp/it);


//attractor declaration ends here



it=1;
//for(sp.ang=sp.tpi*0.5;sp.ang>=sp.tpi*0.25;sp.ang-=0.01)
for(sp.by=354342.216522;sp.by>=1;sp.by/=1.3)//sp.Ih1+1000000
{
//sp.by=sp.Ih1-sp.rt2;
sp.rt2=sp.by-sp.Ih10;
rootmulti(sp);
//sp.rt2=1;//*drand48();
sp.bx=sp.Sh10-sp.ri;//sp.ri*cos(sp.ang)+sp.Sh0;

 sp.rb=sqrt(pow(sp.Sh10-sp.bx,2)+pow(sp.Ih10-sp.by,2)+pow(0,2));
Pr[it][0]=sp.rb;
Pr[it][1]=sp.bx;
Pr[it][2]=sp.by;
it++;

//printf("%lf %lf %lf %d %lf %lf\n",sp.rt2,sp.ri,sp.rb,it,sp.bx,sp.by);

}
nt1=it-1;
if(it==1)
nt1=1;
bubbleSort(Pr,nt1);
sp.Prr1=Pr[nt1][0];//shortest distance to basin from first attractor
//printf("%lf\n",sp.Prr1);



//calculating shortest distance to limit cycle/second attractor given that we now know the basin boundary of the first attractor

//The first step is to calculate points on the second attractor
y0[1]=sp.SSh1[1];y0[2]=sp.EEh1[1]; y0[3]=sp.IIh1[1]; y0[4]=sp.RRh1[1];
y0[5]=sp.SSh2[1];y0[6]=sp.IIh2[1];
y0[7]=sp.SSh3[1];y0[8]=sp.IIh3[1];
y0[9]=sp.SSv[1];y0[10]=sp.EEv[1];y0[11]=sp.IIv[1];



sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
sp.NNav=sp.alph[1]*sp.NH+sp.alph[2]*sp.N2+sp.alph[3]*sp.N3;
for(i=1+id;i<=N;i++)
y0[i]=drand48();

t=0;t1=800000;

it=0;int it1=1;
while (t < t1)
{
sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
sp.NNav=sp.alph[1]*sp.NH+sp.alph[2]*sp.N2+sp.alph[3]*sp.N3;
int status = gsl_odeiv_evolve_apply (e, c, s,
&sys,
&t, t1,
&h, y0);
gs(y0,sum);
if (status != GSL_SUCCESS)
break;
//printf ("%.5e %.5e %.5e  %.5e %.5e\n", t/365., log10(y[0]), log10(y[1]),log10(y[2]),log10(y[3]));
if(t>700000)
{

//it++;
//if(it%10==0 && it1<atl){
//yy[it1][1]=y0[1];
//yy[it1][2]=y0[3];
//it1++;
sp.xmp1=fabs(100*y0[3]/sp.NH);

}
}//while loop for integration ends here

y0[1]=sp.SSh1[1];y0[2]=sp.EEh1[1]; y0[3]=sp.IIh1[1]; y0[4]=sp.RRh1[1];
y0[5]=sp.SSh2[1];y0[6]=sp.IIh2[1];
y0[7]=sp.SSh3[1];y0[8]=sp.IIh3[1];
y0[9]=sp.SSv[1];y0[10]=sp.EEv[1];y0[11]=sp.IIv[1];


it=1;int yes=1,yes1=1;
for(sp.by=sp.IIh1[1]+1000000;sp.by>=1;sp.by/=1.3)
{
//sp.by=sp.Ih1-sp.rt2;
sp.rt2=sp.IIh1[1]-sp.by;
if(sp.by>0 )
{
rootmultiend(sp);
sp.bx=sp.SSh1[1]+sp.ri;//sp.ri*cos(sp.ang)+sp.Sh0;
sp.NHic=sp.bx+sp.EEh1[1]+sp.by+sp.RRh1[1];
if(sp.NHic<=sp.Sh10)
  {
if(sp.ri==0)
  {
   yes++;
  }
   sp.rb2=sqrt(pow(sp.SSh1[1]-sp.bx,2)+pow(sp.IIh1[1]-sp.by,2));
   Pr2[it][0]=sp.rb2;
   Pr2[it][1]=sp.bx;
   Pr2[it][2]=sp.by;

   it++;
//printf("%lf %lf %lf %lf\n",sp.rt2,sp.ri,sp.bx,sp.by);
  }//if condition over testing the population is within sp.PIH1/sp.muH1
}//if condition testing the infectious population is in bilogical meaningful range of greater than zero

}//for loop over rt2 ends here

nt1=it-1;
if(it==1)
nt1=1;

//bubbleSort(Pr,nt1);

//sp.Prr1=Pr[nt1][0];//shortest distance to basin from first attractor

bubbleSort(Pr2,nt1);
sp.Prr2=Pr2[nt1][0];//shortest distance to basin from second attractor

//printf("%lf\n",sp.Prr1);
//printf("%lf\n",sp.Prr2);

//}//this bracket is felxibele and closes the loop over sp.R0



//third step is to calculate local Lyapuno exponents for trajectories staring from the basin boundary and going towards respective attractors and summ all of them for each attractor

t=0;t1=npart;s1=0;s2=0;
for(int ii=1;ii<=nt1;ii++)
{
//first attractor LE
y0[2]=0; y0[4]=0;
y0[5]=sp.Sh20;y0[6]=0;
y0[7]=sp.Sh30;y0[8]=0;
y0[9]=sp.Sv0;y0[10]=0;y0[11]=0;

sp.bx=Pr[ii][1]+Pr[ii][1]/10;
sp.by=Pr[ii][2]-Pr[ii][2]/10;
y0[1]=sp.bx;y0[3]=sp.by;
for(i=1+id;i<=N;i++)
y0[i]=drand48();
for(i=1;i<=M;i++)
     {
      sum[i]=0.0;
     }
for (int iter = 1; iter <= npart; iter++)
{
ti=iter*t1/npart;
while (t < ti)
{
int status = 
gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y0);
//if (h < hmin) { h = hmin; } 


if (status != GSL_SUCCESS)
break;
}
gs(y0,sum);
}
for(i=1;i<=M;i++)
{
sum[i]=sum[i]/t1;
}
s1+=-sum[1];
//printf("%lf\n",y0[3]*100/sp.NH);
//for second attractor
y0[2]=sp.EEh1[1];  y0[4]=sp.RRh1[1];
y0[5]=sp.SSh2[1];y0[6]=sp.IIh2[1];
y0[7]=sp.SSh3[1];y0[8]=sp.IIh3[1];
y0[9]=sp.SSv[1];y0[10]=sp.EEv[1];y0[11]=sp.IIv[1];
sp.bx=Pr2[ii][1]-Pr2[ii][1]/10.0;
sp.by=Pr2[ii][2]+Pr2[ii][2]/10.0;
y0[1]=sp.bx;y0[3]=sp.by;
for(i=1+id;i<=N;i++)
y0[i]=drand48();
for(i=1;i<=M;i++)
     {
      sum[i]=0.0;
     }
for (int iter = 1; iter <= npart; iter++)
{
ti=iter*t1/npart;
while (t < ti)
{
int status = 
gsl_odeiv_evolve_apply (e, c, s, &sys, &t, ti, &h, y0);
//if (h < hmin) { h = hmin; } 


if (status != GSL_SUCCESS)
break;
}
gs(y0,sum);
}
for(i=1;i<=M;i++)
{
sum[i]=sum[i]/t1;
}
s2+=-sum[1];
//printf("%lf\n",y0[3]*100/sp.NH);
}

if(yes==it)
  {
   sp.Prr2=1;
   sp.Prr1=0;
  }
//printf("%lf %lf %lf %lf %d %d\n",sp.Prr1,sp.Prr2,s1,s2,yes,it);
sp.SI=sp.Prr1*s1/(sp.Prr2*s2+sp.Prr1*s1);
fprintf(fp,"%lf %lf %lf %lf %lf\n",sp.R0,sp.R0M,sp.NS,sp.SI,sp.Prr2*s2/(sp.Prr2*s2+sp.Prr1*s1));

}//if condition over sadle node condition
}



gsl_odeiv_evolve_free (e);
gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);
fclose(fp);
}


double DeterMmultihost(sys_para &sp,double A)
{
int l,i,it;
double In,nrm,Inn;


const gsl_odeiv_step_type * T
= gsl_odeiv_step_rkck;
gsl_odeiv_step * s
= gsl_odeiv_step_alloc (T, N);
gsl_odeiv_control * c
= gsl_odeiv_control_y_new (1e-5, 0.0);
gsl_odeiv_evolve * e
= gsl_odeiv_evolve_alloc (N);

gsl_odeiv_system sys = {func, jac, N, &sp};
double t = 0.0, t1 = 100.0;
double h = 1e-4;
double y0[N+1] = {0,0.1, 0.2, 0.3, 0.1 },sum[M+1];
double d_ep,epmax,z_min,z_max,zz_max,zz[niter+1][id],yy[1001][3];

sp.ri=A;

y0[2]=0;  y0[4]=0;
y0[5]=sp.Sh20;y0[6]=0;
y0[7]=sp.Sh30;y0[8]=0;
 y0[9]=sp.Sv0;y0[10]=0;y0[11]=0;

//sp.bx=sp.ri*cos(sp.ang)+sp.Sh0;
//sp.by=sp.rt2*sin(sp.ang)+sp.Ih0;
sp.bx=sp.Sh10-sp.ri;//sp.ri*cos(sp.ang)+sp.Sh0;
sp.by=sp.Ih10+sp.rt2;//sp.rt2*sin(sp.ang)+sp.Ih0;
y0[1]=sp.bx;
y0[3]=sp.by;
for(i=1+id;i<=N;i++)
y0[i]=drand48();

t=0;t1=800000;

for(i=1;i<=M;i++)
     {
      sum[i]=0.0;
     }
it=0;sp.xm=0;
while (t < t1)
{
sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
int status = gsl_odeiv_evolve_apply (e, c, s,
&sys,
&t, t1,
&h, y0);
gs(y0,sum);
if (status != GSL_SUCCESS)
break;
if(t>700000)
{
//sp.xm+=pow((y0[3]/sp.NH,sp.tpi),2);
it++;
}
}
sp.xm=fabs(y0[3]/sp.NH);//sqrt(sp.xm/it);

gsl_odeiv_evolve_free (e);
gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);

if(fabs(sp.xm-sp.xmp)<=0.001)
In=1;
if(fabs(sp.xm-sp.xmp)>=0.001)
In=-1;


return (In);
}

double DeterMmultihostend(sys_para &sp,double A)
{
int l,i,it;
double In,nrm,Inn;


const gsl_odeiv_step_type * T
= gsl_odeiv_step_rkck;
gsl_odeiv_step * s
= gsl_odeiv_step_alloc (T, N);
gsl_odeiv_control * c
= gsl_odeiv_control_y_new (1e-5, 0.0);
gsl_odeiv_evolve * e
= gsl_odeiv_evolve_alloc (N);

gsl_odeiv_system sys = {func, jac, N, &sp};
double t = 0.0, t1 = 100.0;
double h = 1e-4;
double y0[N+1] = {0,0.1, 0.2, 0.3, 0.1 },sum[M+1];
double d_ep,epmax,z_min,z_max,zz_max,zz[niter+1][id],yy[1001][3];
sp.ri=A;

y0[1]=sp.SSh1[1];y0[2]=sp.EEh1[1]; y0[3]=sp.IIh1[1]; y0[4]=sp.RRh1[1];

y0[5]=sp.SSh2[1];y0[6]=sp.IIh2[1];
y0[7]=sp.SSh3[1];y0[8]=sp.IIh3[1];
     y0[9]=sp.SSv[1];y0[10]=sp.EEv[1];y0[11]=sp.IIv[1];

//sp.bx=sp.ri*cos(sp.ang)+sp.Sh0;
//sp.by=sp.rt2*sin(sp.ang)+sp.Ih0;
sp.bx=sp.SSh1[1]+sp.ri;//sp.ri*cos(sp.ang)+sp.Sh0;
sp.by=sp.IIh1[1]-sp.rt2;//sp.rt2*sin(sp.ang)+sp.Ih0;
y0[1]=sp.bx;
y0[3]=sp.by;
for(i=1+id;i<=N;i++)
y0[i]=drand48();

t=0;t1=990000;

for(i=1;i<=M;i++)
     {
      sum[i]=0.0;
     }
it=0;sp.xm=0;
while (t < t1)
{

int status = gsl_odeiv_evolve_apply (e, c, s,
&sys,
&t, t1,
&h, y0);
gs(y0,sum);
if (status != GSL_SUCCESS)
break;
if(t>890000)
{
//sp.xm+=pow((y0[3]/sp.NH,sp.tpi),2);
it++;
}
}
sp.NS=y0[9]+y0[10]+y0[11];
sp.NH=y0[1]+y0[2]+y0[3]+y0[4];
sp.N2=y0[5]+y0[6];
sp.N3=y0[7]+y0[8];
//printf("%lf\n",sp.NH);

sp.xm=fabs(100*y0[3]/sp.NH);//sqrt(sp.xm/it);

gsl_odeiv_evolve_free (e);
gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);

if(fabs(sp.xm-sp.xmp1)<=0.001)
In=1;
if(fabs(sp.xm-sp.xmp1)>=0.001)
In=-1;


return (In);
}

void rootmulti(sys_para &sp)
{
double yl,yr,ym;
double acc=1e-9;
double Llam,Rlam;

Llam=1000;
Rlam=2000000;

if(DeterMmultihost(sp,Llam)*DeterMmultihost(sp,Rlam)<0)
{

if(DeterMmultihost(sp,Llam)<0)
{
yl=Llam;
yr=Rlam;
}
if(DeterMmultihost(sp,Llam)>0)
{
yl=Rlam;
yr=Llam;
}
do
{
  ym=(yl+yr)/2.;
if(DeterMmultihost(sp,ym)<0)
{
yl=ym;
yr=yr;
}
if(DeterMmultihost(sp,ym)>0)
{
yr=ym;
yl=yl;
}
 }
 while(fabs((yr-yl)/(yr+yl))>=acc);
}
sp.ri=ym;
}

void rootmultiend(sys_para &sp)
{
double yl,yr,ym;
double acc=1e-9;
double Llam,Rlam;

Llam=1000;
Rlam=6000000;

if(DeterMmultihostend(sp,Llam)*DeterMmultihostend(sp,Rlam)<0)
{

if(DeterMmultihostend(sp,Llam)<0)
{
yl=Llam;
yr=Rlam;
}
if(DeterMmultihostend(sp,Llam)>0)
{
yl=Rlam;
yr=Llam;
}
do
{
  ym=(yl+yr)/2.;
if(DeterMmultihostend(sp,ym)<0)
{
yl=ym;
yr=yr;
}
if(DeterMmultihostend(sp,ym)>0)
{
yr=ym;
yl=yl;
}
 }
 while(fabs((yr-yl)/(yr+yl))>=acc);
}
sp.ri=ym;
}

void highpowers(sys_para &sp)
{
int i,j;


sp.Bv=sp.b*sp.betS*sp.sigv*sp.NS0r/(sp.muS*sp.Qv);
//RHS components Nu=O*L
sp.O0=1;
sp.O1=sp.ra+sp.rc+sp.re;
sp.O2=sp.ra*sp.rc+sp.ra*sp.re+sp.rc*sp.re;
sp.O3=sp.ra*sp.rc*sp.re;

sp.L0=sp.R0M*sp.R0M+sp.b*sp.Bv*(sp.m[1]*sp.N0r[1]*sp.alph[1]*(sp.bet2+sp.bet3)+sp.m[2]*sp.N0r[2]*sp.alph[2]*(sp.bet1+sp.bet3)+sp.m[3]*sp.N0r[3]*sp.alph[3]*(sp.bet1+sp.bet2));

sp.L1=sp.b*sp.Bv*(sp.m[1]*sp.N0r[1]*sp.alph[1]*(sp.rc+sp.re)+sp.m[2]*sp.N0r[2]*sp.alph[2]*(sp.ra+sp.re)+sp.m[3]*sp.N0r[3]*sp.alph[3]*(sp.ra+sp.rc))*(sp.bet1+sp.bet2+sp.bet3);

sp.L2=sp.b*sp.Bv*(sp.m[1]*sp.N0r[1]*sp.alph[1]*(sp.rc*sp.re)+sp.m[2]*sp.N0r[2]*sp.alph[2]*(sp.ra*sp.re)+sp.m[3]*sp.N0r[3]*sp.alph[3]*(sp.ra*sp.rc))*(sp.bet1+sp.bet2+sp.bet3);

//LHS Components D=(1+F)^2 + M*(1+F)
sp.M1=sp.b*sp.betS*(sp.m[1]*sp.N0r[1]*sp.alph[1]+sp.m[2]*sp.N0r[2]*sp.alph[2]+sp.m[3]*sp.N0r[3]*sp.alph[3])/(sp.muS*sp.alph[1]);

sp.M2=sp.b*sp.betS*(sp.m[1]*sp.N0r[1]*sp.alph[1]*(sp.rc+sp.re)+sp.m[2]*sp.N0r[2]*sp.alph[2]*(sp.ra+sp.re)+sp.m[3]*sp.N0r[3]*sp.alph[3]*(sp.ra+sp.rc))/(sp.muS*sp.alph[1]);

sp.M3=sp.b*sp.betS*(sp.m[1]*sp.N0r[1]*sp.alph[1]*(sp.rc*sp.re)+sp.m[2]*sp.N0r[2]*sp.alph[2]*(sp.ra*sp.re)+sp.m[3]*sp.N0r[3]*sp.alph[3]*(sp.ra*sp.rc))/(sp.muS*sp.alph[1]);

sp.G1=(1/sp.alph[1])*(sp.alph[1]*sp.N0r[1]*sp.FF2[1]*sp.alph[1]+sp.alph[2]*sp.N0r[2]*sp.FF2[2]*sp.alph[2]+sp.alph[3]*sp.N0r[3]*sp.FF2[3]*sp.alph[3]) + sp.alph[1]*sp.N0r[1]*(sp.rc+sp.re)+sp.alph[2]*sp.N0r[2]*(sp.ra+sp.re)+sp.alph[3]*sp.N0r[3]*(sp.ra+sp.rc);


sp.G2=sp.alph[1]*sp.N0r[1]*sp.rc*sp.re+sp.alph[2]*sp.N0r[2]*sp.ra*sp.re+sp.alph[3]*sp.N0r[3]*sp.ra*sp.rc+(sp.alph[1]*sp.alph[1]*sp.N0r[1]*sp.FF2[1]*(sp.rc+sp.re)+sp.alph[2]*sp.alph[2]*sp.N0r[2]*sp.FF2[2]*(sp.ra+sp.re)+sp.alph[3]*sp.alph[3]*sp.N0r[3]*sp.FF2[3]*(sp.ra+sp.rc))/sp.alph[1];

sp.G3=(sp.alph[1]*sp.alph[1]*sp.N0r[1]*sp.FF2[1]*(sp.rc*sp.re)+sp.alph[2]*sp.alph[2]*sp.N0r[2]*sp.FF2[2]*(sp.ra*sp.re)+sp.alph[3]*sp.alph[3]*sp.N0r[3]*sp.FF2[3]*(sp.ra*sp.rc))/sp.alph[1];

sp.MG1=sp.M1;
sp.MG2=sp.M2+sp.M1*sp.G1;
sp.MG3=sp.M3+sp.M1*sp.G2+sp.M2*sp.G1;
sp.MG4=sp.M1*sp.G3+sp.M2*sp.G2+sp.M3*sp.G1;
sp.MG5=sp.M2*sp.G3+sp.M3*sp.G2;
sp.MG6=sp.M3*sp.G3;

//Deominator D=(1+F)^2 + M*(1+F)
//F = G1+G2+G3 
//M = M1+M2+M3
//F^2 = G1^2 + G2^2 + G3^2 + 2*(G1*G2 + G1*G3 + G2*G3); 
//M*(1+F)=MG1 + MG2 + MG3 + MG4 + MG5 + MG6
//This implies D=1+ 2*F + F^2 + M*(1+F)= 1+ 2*(G1+G2+G3) + (G1^2 + G2^2 + G3^2 + 2*[G1*G2 + G1*G3 + G2*G3] ) + MG1 + MG2 + MG3 + MG4 + MG5 + MG6 
//D = D0 + D1 + D2 + D3 + D4 + D5 + D6
//D0=1; D1=2*G1+MG1; D2= 2*G2+G1^2+MG2; D3=2*G3+ 2*G1*G2+MG3; D4=G2^2 + 2*G1*G3 + MG4; D5=2*G2*G3 + MG5; D6= G3^2 + MG6

//Numerator Nu=(O0+O1+O2+O3)*(L0+L1+L2)
//Nu=Nu0+Nu1+Nu2+Nu3+Nu4+Nu5
//Nu0=O0*L0; Nu1= O0*L1+ O1*L0; Nu2=O0*L2+O1*L1+O2*L0; Nu3=O1*L2+O2*L1+O3*L0; Nu4=O2*L2+O3*L1; Nu5= O3*L2
sp.D0=1;
sp.D1=2*sp.G1+sp.MG1;
sp.D2=2*sp.G2+sp.MG2+sp.G1*sp.G1;
sp.D3=2*sp.G3+2*sp.G1*sp.G2+sp.MG3;
sp.D4=sp.G2*sp.G2+2*sp.G1*sp.G3+sp.MG4;
sp.D5=2*sp.G2*sp.G3+sp.MG5;
sp.D6=sp.G3*sp.G3+sp.MG6;
sp.LHSf=(1+sp.ep1+sp.ep2);

sp.A0=sp.LHSf*sp.D0-(sp.O0*sp.L0);
sp.A1=sp.LHSf*sp.D1-(sp.O0*sp.L1+sp.O1*sp.L0);
sp.A2=sp.LHSf*sp.D2-(sp.O0*sp.L2+sp.O1*sp.L1+sp.O2*sp.L0);
sp.A3=sp.LHSf*sp.D3-(sp.O1*sp.L2+sp.O2*sp.L1+sp.O3*sp.L0);
sp.A4=sp.LHSf*sp.D4-(sp.O2*sp.L2+sp.O3*sp.L1);
sp.A5=sp.LHSf*sp.D5-(sp.O3*sp.L2);
sp.A6=sp.LHSf*sp.D6;
}

void bubbleSort(double norm[][3],int nn) 
    {
      bool swapped = true;
      int j = 0,sq,i,kk;
      const int q=101;
      double tmp,tmp1,tmp2;




      while (swapped) {
            swapped = false;
            j++;
            for (int i = 1; i <=nn - j; i++) {
                  if (norm[i][0] < norm[i + 1][0]) {
                        tmp = norm[i][0];
                        tmp1=norm[i][1];
                        tmp2=norm[i][2];
                        norm[i][0] = norm[i + 1][0];
                        norm[i][1]=norm[i+1][1];
                        norm[i][2]=norm[i+1][2];

                        norm[i + 1][0] = tmp;
                        norm[i+1][1]=tmp1;
                        norm[i+1][2]=tmp2;
                        swapped = true;
                  }
            }
      }
}


void bubbleSortpoly(double norm[][2]) {
      bool swapped = true;
      int j = 0,sq,i,kk;
      const int q=101;
      double tmp,tmp1;




      while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < id1 - j; i++) {
                  if (norm[i][0] < norm[i + 1][0]) {
                        tmp = norm[i][0];
                        tmp1=norm[i][1];
                        norm[i][0] = norm[i + 1][0];
                        norm[i][1]=norm[i+1][1];
                        norm[i + 1][0] = tmp;
                        norm[i+1][1]=tmp1;
                        swapped = true;
                  }
            }
      }
}


void gs(double *y,double *sum)
{
  double a[M+1];
  double bnorm;
  int ii,itt,it,iii;

  for(ii=1;ii<=M;ii++){
    a[ii]=0.0;
  }

  for(itt=1;itt<=M;itt++){
    bnorm=0.0;
    if(itt==1){
      for(it=1;it<=id;it++){
    bnorm=bnorm+pow(y[it+id],2.0);
      }
      for(it=1;it<=id;it++){
    y[it+id]=y[it+id]/sqrt(bnorm);
      }
      sum[itt]=sum[itt]+0.5*log(bnorm);
      bnorm=0.0;
    }
    else{
      for(iii=1;iii<=(itt-1);iii++){
    a[iii]=0.0;
    for(it=1;it<=id;it++){
      a[iii]=a[iii]+y[it+iii*id]*y[it+itt*id];
    }
    for(it=1;it<=id;it++){
      y[it+itt*id]=y[it+itt*id]-a[iii]*y[it+iii*id];
        }
      }
      for(it=1;it<=id;it++){
    bnorm=bnorm+pow(y[it+itt*id],2.0);
      }
      for(it=1;it<=id;it++){
    y[it+itt*id]=y[it+itt*id]/sqrt(bnorm);
      }
      sum[itt]=sum[itt]+0.5*log(bnorm);
      bnorm=0.0;
    }}}
