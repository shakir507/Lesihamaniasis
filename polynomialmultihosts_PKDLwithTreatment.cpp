//Life expectancy of goats is 15-18 years. Some put it at 10-12 years due to some problems. Life expectancy of buffaloes is 25-30 years (in general ) and those of cows is 8 to 10 years in india. Life span of street dogs is 13 to 16 years. Life span of mice 2 to 3 
//years.Chicken 2-3 years. We will consider goats and chicken in addition to humans. This means muH2=1/(10*365) and muH3=1/(12*365) for goats and cows respectively. Area of muzaffarpur Bihar is 93 SqKm. 
//Whereas goat density is 100-150 per SqKm (Reference: Goat reproduction scenario in Bihar, India--A. Dey, S. K. Barari, and B.P.S Yadav). This means there are 125*93=11625 goats on average in muzaffarpur. This comes out to be 11625/4801062=0.002421339. Population of chicken in Hazaribhag is 297007 in broiler and 2314 in desi chicken. We will assume twice this population of chickens for muzaffarpur bihar as the human population is about four times in muzaffarpur. The population of cows in muzaffarpur is 315656 (refence-Brief industrial profile of Muzaffarpur district, page 7). This means that PIH2=12000*muH2 and PIH3=600000*muH3. The singh et. al study observed that goats and to some extent cows are a potential reservoir for VL in India: 31 out of 867 goats and 1 in 161 cows were found to b e carrying Leishmania parasite. About five times more goats implied in testing. This indicates that goats have six times more chance of getting VL parasite.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_poly.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_sort_vector.h>

#define id 6
#define M 1
#define N id
#define nhost 3
#define niter 1000000
#define npart 2000000
#define period 365

class sys_para{
public:
double PIH1,PIH2,PIH3,PIS;// birth rates from host 1 to 3 and sand fly
double NS,N2,N3,Ntot,NNav,NSr;//
double muH1,muH2,muH3,muS;//death rates form host 1 to 3 and sand fly
double sig1,rH,rK,d1,gam1;//sig1 is from E1 to I1,rp is rate of reporting the infectious,d1 is disease induced mortality,gam1 is the rate of tranfer from treated class to PKDL
double sig2,sig3;//sig2 is rate of transfer from E2 to I2,sig3 is rate of transfer from E3 to I3
double e1,tau1,taup1,gam2;//e1 is rate of recovery directly from latent state, tau1 from infectious state, taup1 from reported state,gam2 from PKDL state of host 1
double del1,del2,del3;//rate of loss of immunity from host 1 to 3
double Q1,Q2,Q3,Q4,Q5,Q6,Q7;//proportion of host of populations 1 to 3;
double tau2,tau3,d2,d3;
double R02,R0,R0d;
double C1,C2,C3;
double b,bet1,bet2,bet3,betS;
double sigv;
double F2,F1,A,B,C,lH1,lH2,lV1,lV2,NH1,NH2,Ih1,Ih2,Sh2,Sh1,dlH2,f,g,eta;
double llh1[N+1],llh2[N+1],llh3[N+1],Imllh1[N+1],Imllh2[N+1],Imllh3[N+1],NNH1[N+1],NNH2[N+1],NNH3[N+1],IIh1[N+1],IIh2[N+1],IIh3[N+1];
double ImNNH1[N+1],ImNNH2[N+1],ImNNH3[N+1],ImIIh1[N+1],ImIIh2[N+1],ImIIh3[N+1],SSh1[N+1],SSh2[N+1],SSh3[N+1],ImSSh1[N+1],ImSSh2[N+1],ImSSh3[N+1],ddlH2;
double QQ1,QQ2,QQ3,QQ4,QQ5,R0T,PIST,p1,p2;
double m0,m[nhost+1],FF2[nhost+1],FF1[nhost+1],N0[nhost+1],NH[nhost+1],N0r[nhost+1],Nr[nhost+1],alph[nhost+1],Qv,NN0av,NS0r,NN0rav;
double ra,rc,re;
double Q01,Q02,Q03,Q04,R0M;
double ep1,ep2;
double M1,M2,M3,G1,G2,G3;
double Bv,L0,L1,L2,O0,O1,O2,O3,LHSf,D0,D1,D2,D3,D4,D5,D6,MG1,MG2,MG3,MG4,MG5,MG6;
double A0,A1,A2,A3,A4,A5,A6;
double Norm[N+1][2];
};

void highpowers(sys_para &sp);
void bubbleSort(double norm[][2]);
int
main (void)
{
FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
fp=fopen("polynomial_SEIRSnSEI.dat","w");
fp1=fopen("polynomial_multihost1Re_treatmentPKDL_equalpref.dat","w");
fp2=fopen("polynomial_multihost1Im_treatmentPKDL_equalpref.dat","w");
fp3=fopen("polynomial_multihost2Re_treatmentPKDL_equalpref.dat","w");
fp4=fopen("polynomial_multihost2Im_treatmentPKDL_equalpref.dat","w");
fp5=fopen("polynomial_multihost3Re_treatmentPKDL_equalpref.dat","w");
fp6=fopen("polynomial_multihost3Im_treatmentPKDL_equalpref.dat","w");
int i,j,l;

sys_para sp;
double t = 0.0, t1 = npart,ti=1500;
double h = 0.0031,hmin=1e-3;
double y0[N+1]= {0};
double d_ep,epmax,z_min,z_max,zz_max,nrm;
double cf[214],z[214],*x1,*x0;

//gsl_vector *norm=gsl_vector_alloc(N+1);
gsl_permutation *p=gsl_permutation_alloc(N+1);


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
sp.gam1=0.001;//0.001;//only 12.3 % percent of reporting happend in the indian state of bihar in goverment hospitals and the rest  happened at private clinics and NGOs
          //https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-3156.2006.01647.x
//We are taking 1 % reporting of infectious cases per day i.e. sp.gam1=0.01

sp.f=0.14;//Approximately one in seven aysmptomatics proceeds to have VL

sp.tau1=0.00476;//From infectious stage
sp.rH=1/30.0;//recovery from treatment (one month treatment)
sp.rK=1/547.5;//18 months treament period implies this much rate of recovery from PKDL
sp.g=0;//One in (10) 20 VL infected people proceed to PKDL
sp.eta=1/20.0;//one in 20 VL treated people proceed to PKDL

//Transmission rates
sp.bet1=0.74;//Sandflies to humans
sp.del1=0.000913;//Rate of loss of immunity humans
sp.C1=1;
sp.alph[1]=1;
sp.betS=0.158;//Host to sandflies when no host competence is involved
sp.p1=0.01;//weight of infection in sand flies coming from latent humans
sp.p2=.01;

//Sandflies biting rate
sp.b=0.25;// vector biting rate

//host 2 parameters (goats)
sp.bet2=0.74;
sp.muH2=1/(10.0*365);//life span of goats is 10 years
sp.PIH2=sp.muH2*12000;//check population of goats and their life span in the top most comments.
sp.alph[2]=1;//0.6
sp.C2=1;
sp.d2=0;
sp.del2=sp.del1;



//host 3 parameters (cows)
sp.bet3=0.74;
sp.muH3=1/(12.0*365);//life span of cows is 12 years on average
sp.PIH3=315656*sp.muH3;//Check population and life span of cows in the top most comments
sp.alph[3]=1;//0.1
sp.C3=1;
sp.d3=0;
sp.del3=sp.del1;

sp.ep1=1.0000;
sp.ep2=1;

sp.bet2=sp.ep1*sp.bet2;
sp.bet3=sp.ep2*sp.bet3;
sp.R0=1.4;
//-----------------Basic Reproduction number calculation---------------//
for(sp.R0=1.1;sp.R0>=0.2;sp.R0-=0.000001)
{

sp.Q01=sp.sig1+sp.muH1;
sp.Q02=sp.d1+sp.tau1+sp.muH1;
sp.Q03=sp.del1+sp.muH1;

sp.Q1=sp.sig1+sp.muH1;
sp.Q2=sp.d1+sp.gam1+sp.tau1+sp.muH1;
sp.Q3=sp.rH+sp.muH1;
sp.Q4=sp.rK+sp.muH1;
sp.Q5=sp.del1+sp.muH1;

sp.Q6=sp.d2+sp.del2+sp.muH2;
sp.Q7=sp.d3+sp.del3+sp.muH2;

sp.Qv=sp.sigv+sp.muS;

sp.FF1[1]=((1-sp.f)*sp.sig1*sp.Q2*sp.Q3*sp.Q4+(1-sp.g)*sp.tau1*sp.f*sp.sig1*sp.Q3*sp.Q4+(1-sp.eta)*sp.rH*sp.gam1*sp.f*sp.sig1*sp.Q4+sp.rK*sp.eta*sp.rH*sp.gam1*sp.f*sp.sig1+sp.rK*sp.g*sp.tau1*sp.f*sp.sig1*sp.Q3)/(sp.Q1*sp.Q2*sp.Q3*sp.Q4*sp.Q5);
sp.FF1[2]=1/sp.Q6;
sp.FF1[3]=1/sp.Q7;

sp.FF2[1]=1/sp.Q1+sp.f*sp.sig1/(sp.Q1*sp.Q2)+sp.gam1*sp.f*sp.sig1/(sp.Q1*sp.Q2*sp.Q3)+sp.g*sp.tau1*sp.f*sp.sig1/(sp.Q1*sp.Q2*sp.Q4)+sp.eta*sp.rH*sp.gam1*sp.f*sp.sig1/(sp.Q1*sp.Q2*sp.Q3*sp.Q4)+sp.FF1[1];
sp.FF2[2]=sp.FF1[2];
sp.FF2[3]=sp.FF1[3];

sp.m0=(sp.f*sp.sig1+sp.p1*sp.Q02)/(sp.Q01*sp.Q02);
sp.m[1]=sp.alph[1]*sp.C1*(sp.Q3*sp.Q4*(sp.f*sp.sig1+sp.p1*sp.Q2)+sp.p2*(sp.g*sp.tau1*sp.f*sp.sig1*sp.Q3+sp.eta*sp.rH*sp.gam1*sp.f*sp.sig1))/(sp.Q1*sp.Q2*sp.Q3*sp.Q4);
sp.m[2]=sp.alph[2]*sp.C2/sp.Q6;
sp.m[3]=sp.alph[3]*sp.C3/sp.Q7;

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
cf[N]=sp.A6;//sextic term
cf[N-1]=sp.A5;//quintic term
cf[N-2]=sp.A4;//quartic term
cf[N-3]=sp.A3;//cubic term
cf[N-4]=sp.A2;//quadratic term
cf[N-5]=sp.A1;//linear term
cf[N-6]=sp.A0;//constant term
gsl_poly_complex_workspace *w=gsl_poly_complex_workspace_alloc(N+1);

gsl_poly_complex_solve(cf,N+1,w,z);
gsl_poly_complex_workspace_free(w);





//-------------Solving the quintic polynomial-------------------------//

//if(z[1]==0 || z[3]==0 || z[5]==0 || z[7]==0 || z[9]==0 || z[11]==0)
{

for(l=0;l<N;l++)
{nrm=z[2*l];
sp.Norm[l][0]=nrm;
sp.Norm[l][1]=l;
//printf("%lf %lf\n",sp.Norm[l][0],z[2*l]);
}

bubbleSort(sp.Norm);

printf("%lf\n",sp.A0);

if(z[2*int(sp.Norm[0][1])+1]==0 || z[2*int(sp.Norm[1][1])+1]==0 && sp.R0M>0.7)
{
fprintf(fp1,"%lf %lf %lf ",sp.R0,sp.R0M,sp.PIS/sp.muS);
fprintf(fp2,"%lf %lf %lf ",sp.R0,sp.R0M,sp.PIS/sp.muS);
fprintf(fp3,"%lf %lf %lf ",sp.R0,sp.R0M,sp.PIS/sp.muS);
fprintf(fp4,"%lf %lf %lf ",sp.R0,sp.R0M,sp.PIS/sp.muS);
fprintf(fp5,"%lf %lf %lf ",sp.R0,sp.R0M,sp.PIS/sp.muS);
fprintf(fp6,"%lf %lf %lf ",sp.R0,sp.R0M,sp.PIS/sp.muS);
for(l=0;l<2;l++)
{
int  ll=sp.Norm[l][1];
//printf("%lf %lf\n",sp.Norm[l][0],z[2*ll]);
sp.llh1[l+1]=z[2*ll];
sp.llh2[l+1]=sp.alph[2]*z[2*ll]/sp.alph[1];
sp.llh3[l+1]=sp.alph[3]*z[2*ll]/sp.alph[1];

sp.SSh1[l+1]=sp.PIH1/(sp.llh1[l+1]*(1-sp.del1*sp.FF1[1])+sp.muH1);
sp.SSh2[l+1]=sp.PIH2/(sp.llh2[l+1]*(1-sp.del2*sp.FF1[2])+sp.muH2);
sp.SSh3[l+1]=sp.PIH3/(sp.llh3[l+1]*(1-sp.del3*sp.FF1[3])+sp.muH3);

sp.IIh1[l+1]=sp.f*sp.sig1*sp.llh1[l+1]*sp.SSh1[l+1]/(sp.Q1*sp.Q2);
sp.IIh2[l+1]=sp.llh2[l+1]*sp.SSh2[l+1]*sp.FF1[2];
sp.IIh3[l+1]=sp.llh3[l+1]*sp.SSh3[l+1]*sp.FF1[3];

sp.NNH1[l+1]=(1+sp.llh1[l+1]*sp.FF2[1])*sp.SSh1[l+1];
sp.NNH2[l+1]=(1+sp.llh2[l+1]*sp.FF2[2])*sp.SSh2[l+1];
sp.NNH3[l+1]=(1+sp.llh3[l+1]*sp.FF2[3])*sp.SSh3[l+1];


sp.Imllh1[l+1]=z[2*ll+1];
sp.Imllh2[l+1]=sp.alph[2]*z[2*ll+1]/sp.alph[1];
sp.Imllh3[l+1]=sp.alph[3]*z[2*ll+1]/sp.alph[1];

sp.ImSSh1[l+1]=sp.PIH1/(sp.Imllh1[l+1]*(1-sp.del1*sp.FF1[1])+sp.muH1);
sp.ImSSh2[l+1]=sp.PIH2/(sp.Imllh2[l+1]*(1-sp.del2*sp.FF1[2])+sp.muH2);
sp.ImSSh3[l+1]=sp.PIH3/(sp.Imllh3[l+1]*(1-sp.del3*sp.FF1[3])+sp.muH3);

sp.ImIIh1[l+1]=sp.f*sp.sig1*sp.Imllh1[l+1]*sp.ImSSh1[l+1]/(sp.Q1*sp.Q2);
sp.ImIIh2[l+1]=sp.Imllh2[l+1]*sp.ImSSh2[l+1]*sp.FF1[2];
sp.ImIIh3[l+1]=sp.Imllh3[l+1]*sp.ImSSh3[l+1]*sp.FF1[3];

sp.ImNNH1[l+1]=(1+sp.Imllh1[l+1]*sp.FF2[1])*sp.ImSSh1[l+1];
sp.ImNNH2[l+1]=(1+sp.Imllh2[l+1]*sp.FF2[2])*sp.ImSSh2[l+1];
sp.ImNNH3[l+1]=(1+sp.Imllh3[l+1]*sp.FF2[3])*sp.ImSSh3[l+1];
if(z[2*ll+1]==0 && isfinite(sp.IIh1[l+1]/sp.NNH1[l+1])!=0 && sp.R0M>0.7)
{
fprintf(fp1," %lf",100*sp.IIh1[l+1]/sp.NNH1[l+1]);
fprintf(fp2," %lf",sp.ImIIh1[l+1]/sp.ImNNH1[l+1]);
fprintf(fp3," %lf",100*sp.IIh2[l+1]/sp.NNH2[l+1]);
fprintf(fp4," %lf",sp.Imllh2[l+1]/sp.ImNNH2[l+1]);
fprintf(fp5," %lf",100*sp.IIh3[l+1]/sp.NNH3[l+1]);
fprintf(fp6," %lf",sp.Imllh3[l+1]/sp.ImNNH3[l+1]);
}
}
//nrm=sqrt(pow(z[2*l],2)+pow(z[(2*l)+1],2));
  //gsl_vector_set(norm,l,nrm);

fprintf(fp1,"\n");
fprintf(fp2,"\n");
fprintf(fp3,"\n");
fprintf(fp4,"\n");
fprintf(fp5,"\n");
fprintf(fp6,"\n");
}
}
//--------------------------------------------------------------------//

sp.lH1=(-sp.B+sqrt(sp.B*sp.B-4*sp.A*sp.C))/(2*sp.A);
sp.lH2=(-sp.B-sqrt(sp.B*sp.B-4*sp.A*sp.C))/(2*sp.A);
sp.Sh1=sp.PIH1/(sp.lH1*(1-sp.F1)+sp.muH1);
sp.Sh2=sp.PIH1/(sp.lH2*(1-sp.F1)+sp.muH1);
sp.NH1=(1+sp.lH1*sp.F2)*sp.Sh1;
sp.NH2=(1+sp.lH2*sp.F2)*sp.Sh2;
sp.Ih1=sp.f*sp.sig1*sp.lH1*sp.Sh1/(sp.Q01*sp.Q02);

sp.Ih2=sp.f*sp.sig1*sp.lH2*sp.Sh2/(sp.Q01*sp.Q02);


//printf("%lf %lf\n",sp.R0,sp.NS0r);
if(sp.B*sp.B-4*sp.A*sp.C>0)
{
sp.dlH2=(2*sp.R0)/(sp.B);
fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",sp.R0,sp.PIS/sp.muS,sp.lH1,sp.lH2,100*sp.Ih1/sp.NH1,100*sp.Ih2/sp.NH2,sp.dlH2);
}

//if(sp.B<0)
//break;
}



fclose(fp);
fclose(fp1);
fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp5);
fclose(fp6);
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



void bubbleSort(double norm[][2]) {
      bool swapped = true;
      int j = 0,sq,i,kk;
      const int q=101;
      double tmp,tmp1;




      while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < N - j; i++) {
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
