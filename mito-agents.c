// code to simulate ATP concentration dynamics in a model 2D cell

// assume cell production is 10^9 ATP/s (https://www.molbiolcell.org/doi/10.1091/mbc.e14-09-1347)
// assume all OXPHOS: 100 mitos -> production (balanced) = 10^7 ATP/s/mito
// diffusion rate somewhere around (5?)*5*1.5e2 length^2 / s
// solution concentration ~ 10^5 ATP/length^2
// typical cellular ATP conc = 1-10mM (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8615055/)
// 1mM ATP = 10^-3 mol dm-3 = 10^-3 mol (10^15 um^3)-1 = 1e-3 * 6e23 / 1e15 = 6e5 molecules / um^3
// suggest picturing a 10um-thick cell reflected by the 2D profile we simulate: 100x100x10 um = 1e5, for around 6e10 molecules / cell
// update: cell dimensions are now parameters, so we can picture e.g. 50x50x10um

// Intracellular diffusion gradients of O2 and ATP https://journals.physiology.org/doi/epdf/10.1152/ajpcell.1986.250.5.C663

// lots of diffusion rates https://book.bionumbers.org/what-are-the-time-scales-for-diffusion-in-cells/
// generally cell is ~3 times slower than water
// fish muscle: ATP diffusion 2.5e-6 cm^2 / s = 2.5e2 um^2 / s (https://pubmed.ncbi.nlm.nih.gov/7547189/)
// water: 7e-6 cm^2 / s = 7e2 um^2 / s https://www.sciencedirect.com/science/article/abs/pii/0003986164902656
// so something like 2.5e2 um^2/s seems fine

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXNMITO 1000  // lazy, for memory allocation

int NMITO = 100;     // number of mitos in cell
int DEPTH = 10;      // depth of cell in um
int SUBDIV = 2;      // number of simulation cells in 1um3, controls grid granularity

#define RND drand48()

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

int myround(double x)
{
  int xi = (int)x;
  if(x-xi > 0.5) return xi+1;
  return xi;
}

int main(int argc, char *argv[])
{
  double *grid, *newgrid;
  double x[MAXNMITO], y[MAXNMITO];
  FILE *fp;
  char fstr[100];
  double loss, totalloss, totalgain;
  double h, l, r, u, d;
  int i, j;
  int m;
  double t;
  double dt = 0.00025;
  int lastt;
  double total;
  int expt;
  double totalATP, totalATP2, changeATP, minATP, maxATP, lastATP, meanATP, meanATP2, varATP;
  double kappa, delta;
  double vol;
  char masterfstr[100], snapfstr[100], mitosnapfstr[100];
  double mparam = 0;
  double dx, dy, gradx, grady, norm;
  int cellsize;
  int GRIDX, GRIDY;
  int CENTRE, CENTRESCALE;
  double D = 2.5e2;
  double speed;
  double gridx;

  // process command-line arguments
  // different experiments: (0) uniform mitos, uniform consumption; (1) uniform mitos, clustered consumption;
  //                        (2) clustered mitos, uniform consumption; (3) clustered mitos, clustered consumption
  //                        (4) mobile-1 mitos, uniform consumption; (5) mobile-1 mitos, clustered consumption;
  //                        (6) mobile-2 mitos, uniform consumption; (7) mobile-2 mitos, clustered consumption
  //                        (8) mobile-2 mitos, uniform consumption; (9) mobile-2 mitos, clustered consumption
  if(argc < 6) {
    printf("Please specify a cell size, cell depth, number of mitos, elements per um, and an experiment! Experiments 0-9\n");
    return 0;
  }
  SUBDIV = atoi(argv[4]);
  expt = atoi(argv[5]);
  if(expt >= 4) {
    if(argc < 7) {
      printf("This experiment needs a motion parameter!\n");
      return 0;
    }
    mparam = atof(argv[6]);
  }
  // physical properties of cell
  cellsize = atoi(argv[1]);
  DEPTH = atoi(argv[2]);
  NMITO = atoi(argv[3]);
  if(NMITO >= MAXNMITO) {
    printf("Too many mitochondria!\n");
    return 0;
  }

  // choose a timestep appropriate for our grid element scale (for numerical stability)
  if(SUBDIV <= 1) dt = 0.001;
  else if(SUBDIV <= 2) dt = 0.00025;
  else if(SUBDIV <= 4) dt = 0.00005;

  // helper geometric quantities
  gridx = 1./SUBDIV;
  GRIDX = GRIDY = cellsize*SUBDIV;
  CENTRE = GRIDX/2;
  
  // circle radius CENTRE/2 -> pi (grid/4)^2 vs grid^2 -> 1./(pi/16.) 
  CENTRESCALE = 1./(3.1416/16.);
  
  // allocate memory to store GRIDX*GRIDY grid of ATP concentration values
  grid = (double*)malloc(sizeof(double)*GRIDX*GRIDY);
  newgrid = (double*)malloc(sizeof(double)*GRIDX*GRIDY);

  // open file for output. we'll store statistics of each situation by timestep; we'll also (later) store snapshots of individual cases
  sprintf(masterfstr, "stats-%i-%i-%i-%i-%i.csv", expt, cellsize, DEPTH, NMITO, SUBDIV);
  fp = fopen(masterfstr, "w");
  fprintf(fp, "expt,mparam,kappa,delta,t,total.ATP,ex.molar,min.ATP,max.ATP,change.ATP,mean.ATP,var.ATP,vol.dm3,consumption,terminated\n");
  fclose(fp);

  // total cell volume in dm-3: (in metres) * (dm3 in 1 m3)
  vol = (cellsize*1e-6 * cellsize*1e-6 * DEPTH*1e-6) * (10*10*10);
	
  // kappa is the rate constant of ATP consumption
  for(kappa = 0.01; kappa <= 10; kappa *= 2)
    {
      // delta is the production rate per mitochondrion
      for(delta = 0.01; delta <= 10; delta *= 2)
	{
	  // output tracking info and initialise random seed for reproducibility
	  printf("Running expt %i, kappa %.2e, delta %.2e\n", expt, kappa, delta);
	  srand48(1);
	      
	  // initialises positions of mitochondria (x, y are coordinates)
	  for(m = 0; m < NMITO; m++)
	    {
	      // if expt == 1, cluster near centre; otherwise spead evenly
	      if(expt == 0 || expt == 1) { x[m] = RND*GRIDX; y[m] = RND*GRIDY; }
	      if(expt == 2 || expt == 3) { x[m] = (RND-0.5)*CENTRE+GRIDX/2; y[m] = (RND-0.5)*CENTRE+GRIDY/2; }
	      if(expt == 4 || expt == 5) { x[m] = RND*GRIDX; y[m] = RND*GRIDY; }
	      if(expt == 6 || expt == 7) { x[m] = RND*GRIDX; y[m] = RND*GRIDY; }
	      if(expt == 8 || expt == 9) { x[m] = RND*GRIDX; y[m] = RND*GRIDY; }
	    }

	  // initialise cell with zero ATP
	  for(i = 0; i < GRIDX; i++)
	    {
	      for(j = 0; j < GRIDY; j++)
		grid[i*GRIDY+j] = 0;
	    }

	  // initialise change trackers
	  lastt = -1; changeATP = 1; totalATP = lastATP = 0;

	  // initialise output files
	  if(expt >= 4)
	    {
	      sprintf(mitosnapfstr, "out-mitos-%i-%i-%i-%i-%.2f-%.2f-%i.txt", expt, cellsize, DEPTH, NMITO, kappa, delta, SUBDIV);
	      fp = fopen(mitosnapfstr, "w");
	      fclose(fp);
	    }
	  sprintf(snapfstr, "out-%i-%i-%i-%i-%.2f-%.2f-%i.txt", expt, cellsize, DEPTH, NMITO, kappa, delta, SUBDIV);
	  fp = fopen(snapfstr, "w");
	  fclose(fp);
			  
	  // loop until a time threshold or until equilibration criterion is met
	  for(t = 0; t < 1000.1 && changeATP > dt*1e-4*totalATP; t+=dt)
	    {
	      totalloss = totalgain = 0;
	      // actual PDE "solver": loop through each element of our discretised domain
	      for(i = 0; i < GRIDX; i++)
		{
		  for(j = 0; j < GRIDY; j++)
		    {
		      // find neighbours of this cell element -- using reflecting boundary conditions
		      l = (i == 0 ? grid[i*GRIDY+j] : grid[(i-1)*GRIDY+j]);
		      r = (i == GRIDX-1 ? grid[i*GRIDY+j] : grid[(i+1)*GRIDY+j]);
		      u = (j == 0 ? grid[i*GRIDY+j] : grid[i*GRIDY+j-1]);
		      d = (j == GRIDY-1 ? grid[i*GRIDY+j] : grid[i*GRIDY+j+1]);
		      h = grid[i*GRIDY+j];
			  
		      // finite difference solver with some parameters to capture diffusion
		      newgrid[i*GRIDY+j] = grid[i*GRIDY+j] + (D/(gridx*gridx))*(dt*(l + r - 2*h) + dt*(u + d - 2*h));
			  
		      // loss of ATP. in proportion to current concentration; sites of loss depend on model
		      if(expt % 2 == 0)
			{
			  // loss throughout cell
			  loss = dt* newgrid[i*GRIDY+j]*kappa;
			}
		      if(expt % 2 == 1)
			{
			  // loss only in central circular region
			  if((i-CENTRE)*(i-CENTRE) + (j-CENTRE)*(j-CENTRE) < (CENTRE*CENTRE/4))
			    loss = dt* newgrid[i*GRIDY+j]*kappa*CENTRESCALE;
			  else
			    loss = 0;
			}
		      // impose nonnegativity (crudely, and shouldn't happen)
		      if(newgrid[i*GRIDY+j]-loss < 0) newgrid[i*GRIDY+j] = 0;
		      else { newgrid[i*GRIDY+j] -= loss; totalloss += loss; }
		    }
		}

	      // move our mitochondria, if required
	      if(expt >= 4) {
		for(m = 0; m < NMITO; m++)
		  {
		    i = myround(x[m]); j = myround(y[m]);
		    // find neighbours of this mito -- using periodic boundary conditions
		    l = (i == 0 ? grid[0*GRIDY+j] : grid[(i-1)*GRIDY+j]);
		    r = (i == GRIDX-1 ? grid[i*GRIDY+j] : grid[(i+1)*GRIDY+j]);
		    u = (j == 0 ? grid[i*GRIDY+j] : grid[i*GRIDY+j-1]);
		    d = (j == GRIDY-1 ? grid[i*GRIDY+j] : grid[i*GRIDY+j+1]);
		    h = grid[i*GRIDY+j];

		    if(expt == 4 || expt == 5) {
		      // random diffusion with width mparam * [ATP] um / second
		      dx = gsl_ran_gaussian(mparam*h/(gridx*gridx));
		      dy = gsl_ran_gaussian(mparam*h/(gridx*gridx));
		    }
		    if(expt == 8 || expt == 9) {
		      // random diffusion with width mparam 
		      dx = gsl_ran_gaussian(mparam);
		      dy = gsl_ran_gaussian(mparam);
		    }
		    if(expt == 6 || expt == 7) {
		      gradx = (r-l)/2.;
		      grady = (u-d)/2.;
		      norm = sqrt(gradx*gradx + grady*grady);
		      if(norm == 0)
			{
			  // some small perturbation if concentration is completely flat (shouldn't happen)
			  dx = gsl_ran_gaussian(1./gridx);
			  dy = gsl_ran_gaussian(1./gridx);
			}
		      else {
			// step size U(0,1) * mparam * [x component of gradient] um / second
			double rnd = RND;
			dx = -(mparam/gridx)*rnd*gradx/norm;
			dy = (mparam/gridx)*rnd*grady/norm;
			//	if(t > 400)
			//  printf("%f %i: %f , l %f r %f d %f u %f, dx %f dy %f\n", t, m, h, l, r, d, u, dx, dy);
		      }
		    }
		    // we want to limit to 1 um/s maximum, but do this through empirical choice of the mparam values, not explicitly here
		    x[m] += dx*dt;
		    y[m] += dy*dt;
		    if(x[m] > 2*GRIDX || y[m] > 2*GRIDY || x[m] < -GRIDX || y[m] < -GRIDY) {
		      printf("Way out of bounds!\n");
		      return 0;
		    }
		    // reflecting boundaries
		    if(x[m] > GRIDX-1) x[m] = GRIDX-1;
		    if(x[m] < 0) x[m] = 0;
		    if(y[m] > GRIDY-1) y[m] = GRIDY-1;
		    if(y[m] < 0) y[m] = 0;
		  }
	      }
	      // the above took care of loss and diffusion, now take care of production terms
	      // loop through the NMITO mitos and add point ATP mass at each position
	      for(m = 0; m < NMITO; m++)
		{
		  newgrid[myround(x[m])*GRIDY+myround(y[m])] += delta*1e7*dt;
		  // for reference, consumption is 10^9 per cell per sec and we have (default) 10^2 mitos
		  // we're scaling through orders of magnitude for delta though so will see if we're in a reasonable region by tracking the cellular total 
		  totalgain += delta*1e7*dt;
		}

	      // finally, update new state of cell from buffer
	      totalATP = 0; totalATP2 = 0; minATP = 1e100; maxATP = 0; 
	      for(i = 0; i < GRIDX; i++)
		{
		  for(j = 0; j < GRIDY; j++)
		    {
		      grid[i*GRIDY+j] = newgrid[i*GRIDY+j];
		      if(grid[i*GRIDY+j] < minATP) minATP = grid[i*GRIDY+j];
		      if(grid[i*GRIDY+j] > maxATP) maxATP = grid[i*GRIDY+j];
		      // contribute to cellular ATP statistics (total/mean and variance)
		      totalATP += grid[i*GRIDY+j];
		      totalATP2 += grid[i*GRIDY+j]*grid[i*GRIDY+j];
		    }
		}
	      // record change in ATP for equilibration tracker
	      changeATP = totalATP - lastATP;
	      lastATP = totalATP;

	      // statistics of ATP: mean and variance
	      meanATP = totalATP/(GRIDX*GRIDY); meanATP2 = totalATP2/(GRIDX*GRIDY);
	      varATP = meanATP2 - meanATP*meanATP; 

	      // if we're in a new timestep
	      if((int)t != lastt)
		{
		  // output information, optionally
		  // printf("Total loss in last second %.2e\n", totalloss/dt);
		  //printf("vol is %.2e dm3\n", vol);
		  //printf("Total ATP = %.2e molecules (%.2e mol) (conc %.2eM)\n", totalATP, totalATP/6e23, totalATP/6e23 / vol);

		  // output stats snapshot
		  // ATP concentrations here are output as molecule number per simulation cell
		  // needs accounting for if we compare different simulation cell sizes or depths
		  fp = fopen(masterfstr, "a");
		  fprintf(fp, "%i,%e,%.2e,%.2e,%i,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,0\n", expt, mparam, kappa, delta, (int)t, totalATP, totalATP/6e23/vol, minATP, maxATP, changeATP, meanATP, varATP, vol, totalloss/dt);
		  fclose(fp);

		  if(expt >= 4)
		    {
		      fp = fopen(mitosnapfstr, "a");
		      for(m = 0; m < NMITO; m++)
			fprintf(fp, "%i %i %f %f\n", (int) t, m, x[m]*gridx, y[m]*gridx);
		      fclose(fp);
		    }
		  
		  // take full snapshots of early behaviour and subsequent changes
		  if((int)t < 5 || (int)t % 100 == 0)
		    {
		      fp = fopen(snapfstr, "a");
		      for(i = 0; i < GRIDX; i++)
			{
			  for(j = 0; j < GRIDY; j++)
			    {
			      fprintf(fp, "%i %.4f %.4f %.5f\n", (int) t, (double)i*gridx, (double)j*gridx, grid[i*GRIDY+j]);
			      total += grid[i*GRIDY+j];
			    }
			  fprintf(fp, "\n");
			}
		      fclose(fp);
		    }
		  lastt = ((int)t);
		}

	    }

	  // this run terminated, either through equilibration or time runout
	  printf("stopped at %e\n", t);
	  fp = fopen(masterfstr, "a");
	  fprintf(fp, "%i,%e,%.2e,%.2e,%i,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,1\n", expt, mparam, kappa, delta, (int)t, totalATP, totalATP/6e23/vol, minATP, maxATP, changeATP, meanATP, varATP, vol, totalloss/dt);
	  fclose(fp);

	  // output equilibrated state
	  fp = fopen(snapfstr, "a");
	  for(i = 0; i < GRIDX; i++)
	    {
	      for(j = 0; j < GRIDY; j++)
		{
		  // this outputs simply ATP molecules number per simulation cell
		  fprintf(fp, "%i %.4f %.4f %.5f\n", (int) t, (double)i*gridx, (double)j*gridx, grid[i*GRIDY+j]);
		  total += grid[i*GRIDY+j];
		}
	      fprintf(fp, "\n");
	    }
	  fclose(fp);

	  // final output details of mitos
	  sprintf(fstr, "mitos-%i-%i-%i-%i-%i.txt", expt, cellsize, DEPTH, NMITO, SUBDIV);
	  fp = fopen(fstr, "w");
	  for(m = 0; m < NMITO; m++)
	    fprintf(fp, "%.3f %.3f\n", x[m]*gridx, y[m]*gridx);
	  fclose(fp);
	}
    }

  return 0;
}

// 3D
// volume = 100*1e-6 * 100*1e-6 * d*1e-6 = d * 1e-14 m^3
//        = d * 1e-11 dm-3
// generally about 2e11 molecules = 3e-13 mol
// -> 3e-13/(d * 1e-11) mol dm-3
// say d = 10um, then 3e-3 mol dm-3 = 3mM
