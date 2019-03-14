#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <deque>
#include <cmath>
#include <vector>


double A_north(double dy,double R,double a,double v)
{
	return -(1.0/(pow(dy,2)))+(R*a*v)/dy;
}

double A_south(double dy,double R,double a,double v)
{
	return -(1.0/(pow(dy,2)))-(R*a*v)/dy;
}

double B_coeff(double dx, double R, double a, double up, double uw)
{
	return (R*a*up*(up-uw)/dx) - ((pow(a,2))*(uw+up-2.0*up)/(pow(dx,2)));
}

double changing_coeff(double B,double dp)
{
	return B+dp;
};

double u_guess(double dp,double y)
{
	return (0.5)*dp*((pow(y,2))-y);
}

double contCheck(std::vector<double> uSolved)
{
	size_t size = uSolved.size();
	double sum = 0;
	for(int i = 0; i<size;++i)
		sum = sum+uSolved[i]*0.1;
	return sum;
}

double P_j(double An,double Ap, double As, double P)
{
	return An/(Ap-As*P);
}

double Q_j(double D, double As, double Ap, double Q, double P)
{
	return (D+As*Q)/(Ap-As*P);
}

double newU(double P, double Q, double u_new)
{
	return P*u_new + Q;
}

double updateV(double up, double uw)
{
	return up-uw;
}

void TDMA(double u[10], double v[10],double &dp, std::vector<double> &uSolution, double diff, double uw[10])
{
//	double ndp = dp;
	while(diff>0.001 || diff<-0.001)
	{
		double P_old = 0;
		double Q_old = 0;
		double D[10],P[10],Q[10];
		if(diff<0)
			dp = dp+0.01;
		if(diff>0)
			dp = dp-0.01;
		for(int i=0;i<9;++i)
		{
			double An = A_north(0.1,30,0.1,v[i]);
			double As = A_south(0.1,30,0.1,v[i]);
			D[i] = B_coeff(0.1,30,0.1,u[i],uw[i])+dp;
			P[i] = P_j(An,An+As,As,P_old);
			Q[i] = Q_j(D[i],As,As+An,Q_old,P_old);
			P_old = P[i];
			Q_old = Q[i];
		}
		uSolution.clear();
		uSolution.push_back(0);
		for(int i = 9;i>0;--i)
		{
			size_t uSize = uSolution.size();
			uSolution.insert(uSolution.begin(),P[i-1]*uSolution[0]+Q[i-1]);
		}
		diff = 1-contCheck(uSolution);
	}
//	for(int i = 0; i<10;++i)
//	{
//		uw[i] = uSolution[i];
//		u[i] = u_guess(dp,y);
//		v[i] = updateV(...);
//	}
}

int main()
{
	double utop,vtop, vend = 0;
	double uprev[10] = {1,1,1,1,1,1,1,1,1,1};
	double dp = -12;
	double P_old = 0;
	double Q_old = 0;
	double y_range[10], uGuess[10],vGuess[10],D[10],P[10],Q[10];
	std::vector<double> uSoln;
	uSoln.push_back(0);
	// Loop to define coefficients and fill out TDMA
	for(int i=0;i<10;++i)
	{
		y_range[i] = 0.1+(0.1)*i;
		uGuess[i] = u_guess(dp,y_range[i]);
		vGuess[i] = -0.01;
		if(i>3)
			vGuess[i] = 0.01;
		if(i==4)
			vGuess[i] = 0;
		if(i==9)
			vGuess[i] = vend;
		double An = A_north(0.1,30,0.1,vGuess[i]);
		double As = A_south(0.1,30,0.1,vGuess[i]);
		D[i] = B_coeff(0.1,30,0.1,uGuess[i],uprev[i])+dp;
		if(i==9)
			D[i] = 0;
		P[i] = P_j(An,An+As,As,P_old);
		Q[i] = Q_j(D[i],As,As+An,Q_old,P_old);
		if(i==9)
		{
			P[i] = 0;
			Q[i] = 0;
		}
		P_old = P[i];
		Q_old = Q[i];
		
		std::cout<<y_range[i]<<" "<<uGuess[i]<<" "<<vGuess[i]<<" "<<D[i]<<" "<<P[i]<<" "<<Q[i]<<std::endl;
	}
	for(int i = 9;i>0;--i)
	{
		size_t uSize = uSoln.size();
		uSoln.insert(uSoln.begin(),P[i-1]*uSoln[0]+Q[i-1]);
	}

	double continuity = contCheck(uSoln);
	std::cout<<continuity<<std::endl;
	std::cout<<"Old pressure:" <<dp<<std::endl;
	double diff = 0.9-continuity;
	std::cout<<"DIFF: "<<diff<<std::endl;
	
	TDMA(uGuess, vGuess,dp, uSoln, diff, uprev);
	std::cout<<"New pressure: "<<dp<<std::endl;
	for(int i = 0; i<uSoln.size();++i)
		std::cout<<uSoln[i]<<std::endl;

	return 0;
}