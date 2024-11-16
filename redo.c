#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RND drand48()
#define _RK

#define PI 3.14159

double closedform(double t, double omega, double phi, double sigma, double rho)
{
  double c1, c2, lambda1, lambda2;
  double ips = 1 + rho + sigma;
  double W;
  
  lambda1 = 0.5*( -ips + sqrt(ips*ips - 4*(sigma+rho)) );
  lambda2 = 0.5*( -ips - sqrt(ips*ips - 4*(sigma+rho)) );
  c1 = 1. / ( rho/(sigma+lambda1) - rho/(sigma+lambda2) );
  c2 = -c1;

  W = 0;
  W += 0.5*(c1/lambda1 * exp(lambda1*t) + c2/lambda2 * exp(lambda2*t));
  W += 0.5*( c1*exp(lambda1*t)/(lambda1*lambda1+omega*omega) * lambda1*sin(omega*t + phi) - omega*cos(omega*t+phi));
  W += 0.5*( c2*exp(lambda2*t)/(lambda2*lambda2+omega*omega) * lambda2*sin(omega*t + phi) - omega*cos(omega*t+phi));

  return W;
}

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// return max of two values (used for thresholding a zero)
double mymax(double x, double y)
{
  if(x > y) return x;
  return y;
}

// return min of two values (used for thresholding a zero)
double mymin(double x, double y)
{
  if(x < y) return x;
  return y;
}

// environmental function (nutrient available) non-negative sine wave, possible with noise
double envfn(double t, double omega, double phi, double eps)
{
  double core =  0.5+0.5*sin(omega*t + phi);
  if(eps > 0) {
    core += gsl_ran_gaussian(eps);
    return mymax(0, mymin(1, core));
  }
  return core;
}

// simulate an instance of the system
void Simulate(double k0, double kp, double ki, double kd, double beta, double omega, double phi, double eps, double *W, double *L, double *rdelta, double *tend, int output, char *fname)
{
  double stateI, stateA, stateW, stateL;  // state of system
  double dstateIdt, dstateAdt, dstateWdt, dstateLdt;   // d/dt
  double rho, sigma, kappa1, kappa2, delta;   // rate consts and associated quantities
  double t, dt = 0.01;
  FILE *fp;
  double env;
  double diff, prevdiff, ddiffdt, intdiffdt;  // used in PID calcs
 
  if(output)
    {
      fp = fopen(fname, "w");
      fprintf(fp, "t, env, stateI, stateA, stateW, stateL, sigma, rho, kappa1, kappa2, cfW, intdiffdt\n");
    }
      
  // initialise state
  stateI = 1; stateA = stateW = stateL = 0;
  //delta = 0;
  // fix a delta value for this parameterisation
#ifdef _IGJ
  delta = beta*(k0 + kp + ki + kd);
#endif
  // euler time simulation
  prevdiff = 0; intdiffdt = 0; 
  for(t = 0; t < 100 && (stateI + stateA) > 1e-6; t += dt)
    {
      // current environments
      env = envfn(t, omega, phi, eps);
      // environment statistics for PID
      diff = env-0.5;
      ddiffdt = (diff-prevdiff)/dt;
      intdiffdt = intdiffdt + diff*dt;

      // rate constants
      sigma = mymax(0, k0 + kp*diff + ki*intdiffdt + kd*ddiffdt);
      rho = mymax(0, k0 - kp*diff - ki*intdiffdt - kd*ddiffdt);
#ifdef _RK
      delta = mymax(0, beta*(kp*diff + ki*intdiffdt + kd*ddiffdt));
#endif
      kappa1 = mymax(0, env*(1.-delta));
      kappa2 = mymax(0, 1.-env*(1.-delta));
      //delta = delta + beta*(k0 + kp + ki + kd)*dt;

      // time derivatives
      dstateIdt = dt* ( -stateI*sigma + stateA*rho );
      dstateAdt = dt* (  stateI*sigma - stateA*rho - stateA*(kappa1+kappa2) );
      dstateWdt = dt* (  stateA*kappa1 );
      dstateLdt = dt* (  stateA*kappa2 );

      // update state and environment statistic
      stateI += dstateIdt;
      stateA += dstateAdt;
      stateW += dstateWdt;
      stateL += dstateLdt;
      prevdiff = diff;

      if(output)
	{
	  fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", t, env, stateI, stateA, stateW, stateL, sigma, rho, kappa1, kappa2, closedform(t, omega, phi, k0, k0), intdiffdt);
	}

    }
  *W = stateW; *L = stateL; *rdelta = delta; *tend = t;
}

int main(void)
{
  double beta, omega, phi;  // environment and cost
  double k0, kp, ki, kd;    // sensing parameters
  double tomega;
  double stateL, stateW, delta, tend;
  FILE *fp;
  double testomega = 1;
  double eps;
  int teps;
  
  //testomega = 2*PI;
  
  // collections of test cases outputting time series
  beta = 0; omega = 0; phi = 2.5; eps = 0; k0 = 1, kp = 0; ki = 0.0; kd = 0;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-0a.csv");
  beta = 0; omega = testomega; phi = 2.5; eps = 0; k0 = 1; kp = 0; ki = 0.0; kd = 0;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-0b.csv");
  beta = 0; omega = testomega; phi = 0; eps = 0; k0 = 1; kp = 0; ki = 0.0; kd = 0;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-0c.csv");
  beta = 0; omega = testomega; phi = 2.5; eps = 0; k0 = 0.25; kp = 0.75; ki = 0.75; kd = 1;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-1.csv");
  beta = 10; omega = testomega; phi = 2.5; eps = 0; k0 = 0.25; kp = 0.75; ki = 0.75; kd = 1;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-2.csv");

  beta = 0; omega = 0; phi = 4.0; eps = 0; k0 = 1; kp = 0.4; ki = 1; kd = 0.4;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-issue.csv");
  beta = 0; omega = testomega; phi = 2.5; eps = 0.5; k0 = 0.25; kp = 0.75; ki = 0.75; kd = 1;
  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 1, "example-noise.csv");

  // return 0;

#ifdef _RK
  fp = fopen("redo-out-rk.csv", "w");
#else
  fp = fopen("redo-out.csv", "w");
#endif
  fprintf(fp, "beta,epsilon,omega,phi,k0,kp,ki,kd,W,L,delta,tend\n");
  // loop through sensing cost
  for(beta = 0; beta <= 1; beta += 1)
    {
      for(teps = 0; teps < 3; teps++)
	{
	  switch(teps) {
	  case 0: eps = 0; break;
	  case 1: eps = 0.2; break;
	  case 2: eps = 0.5; break;
	  }
	  // loop through environmental frequency
	  //  for(tomega = 0.02; tomega < 5; tomega *= 2)
	  for(tomega = 0; tomega <= 2; tomega += 0.1)
	    {
	      omega = tomega;
	      //  if(tomega == 0.02) omega = 0; else omega = tomega;
	      // loop through environmental phase
	      for(phi = 0; phi < 2*PI; phi += 0.5)
		{
		  // loop through PID terms
		  for(k0 = 0; k0 <= 1; k0 += 0.2)
		    {
		      for(kp = 0; kp <= 1; kp += 0.2)
			{
			  for(ki = 0; ki <= 1; ki += 0.2)
			    {
			      for(kd = 0; kd <= 1; kd += 0.2)
				{
				  Simulate(k0, kp, ki, kd, beta, omega, phi, eps, &stateW, &stateL, &delta, &tend, 0, "tmp");
				  fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", beta, eps, tomega, phi, k0, kp, ki, kd, stateW, stateL, delta, tend);
				}
			    }
			}
		    }
		}
	      printf("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", beta, eps, tomega, phi, k0, kp, ki, kd, stateW, stateL, delta);
	    }
	}
    }
  fclose(fp);
  
  return 0;
}
						     
