// assume cell production is 10^9 ATP/s (https://www.molbiolcell.org/doi/10.1091/mbc.e14-09-1347)
// assume all OXPHOS: 100 mitos -> production (balanced) = 10^7 ATP/s/mito
// diffusion rate somewhere around (5?)*5*1.5e2 length^2 / s
// solution concentration ~ 10^5 ATP/length^2
// typical cellular ATP conc = 1-10mM (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8615055/)
// 1mM ATP = 10^-3 mol dm-3 = 10^-3 mol (10^15 um^3)-1 = 1e-3 * 6e23 / 1e15 = 6e5 molecules / um^3
// suggest picturing a 10um-thick cell reflected by the 2D profile we simulate

// Intracellular diffusion gradients of O2 and ATP https://journals.physiology.org/doi/epdf/10.1152/ajpcell.1986.250.5.C663

// lots of diffusion rates https://book.bionumbers.org/what-are-the-time-scales-for-diffusion-in-cells/
// generally cell is ~3 times slower than water
// fish muscle: ATP diffusion 2.5e-6 cm^2 / s = 2.5e2 um^2 / s (https://pubmed.ncbi.nlm.nih.gov/7547189/)
// water: 7e-6 cm^2 / s = 7e2 um^2 / s https://www.sciencedirect.com/science/article/abs/pii/0003986164902656
// so something like 2.5e2 um^2/s seems fine

#include <stdio.h>
#include <stdlib.h>

#define RND drand48()

int main(void)
{
  double *grid, *newgrid;
  int x[100], y[100];
  FILE *fp;
  char fstr[100];
  double loss, totalloss, totalgain;
  double h, l, r, u, d;
  int i, j;
  double t;
  double dt = 0.001;
  int lastt;
  double total;
  int expt;
  double totalATP, changeATP, minATP, maxATP, lastATP;
  double kappa, delta;
  double vol;
  
  // allocate memory to store 100x100 grid of ATP concentration values
  grid = (double*)malloc(sizeof(double)*100*100);
  newgrid = (double*)malloc(sizeof(double)*100*100);

  fp = fopen("stats.csv", "w");
  fprintf(fp, "expt,kappa,delta,t,total.ATP,ex.molar,min.ATP,max.ATP,change.ATP,consumption,terminated\n");
  fclose(fp);

  // total cell volume in dm-3: (in metres) * (dm3 in 1 m3)
  vol = (100*1e-6 * 100*1e-6 * 10*1e-6) * (10*10*10);
	
  // two different experiments: random mito spread or clustering near centre
  for(expt = 0; expt <= 2; expt++)
    {
      for(kappa = 0.01; kappa <= 10; kappa *= 2)
	{
	  for(delta = 0.01; delta <= 10; delta *= 2)
	    {
	      printf("Running expt %i, kappa %.2e, delta %.2e\n", expt, kappa, delta);
	      srand48(1);
	      // initialises positions of mitochondria (x, y are coordinates)
	      for(i = 0; i < 100; i++)
		{
		  if(expt != 1) { x[i] = RND*100; y[i] = RND*100; }
		  else { x[i] = RND*20+40; y[i] = RND*20+40; }
		  for(j = 0; j < 100; j++)
		    grid[i*100+j] = 0;
		}

	      // discrete time PDE solver (very crude)
	      lastt = -1; changeATP = 1; totalATP = lastATP = 0;
	      for(t = 0; t < 1000.1 && changeATP > dt*1e-4*totalATP; t+=dt)
		{
		  totalloss = totalgain = 0;
		  // actual PDE "solver": loop through each element of our discretised domain
		  for(i = 0; i < 100; i++)
		    {
		      for(j = 0; j < 100; j++)
			{
			  // find neighbours of this cell element -- using periodic boundary conditions
			  l = (i == 0 ? grid[1*100+j] : grid[(i-1)*100+j]);
			  r = (i == 99 ? grid[98*100+j] : grid[(i+1)*100+j]);
			  u = (j == 0 ? grid[i*100+1] : grid[i*100+j-1]);
			  d = (j == 99 ? grid[i*100+98] : grid[i*100+j+1]);
			  h = grid[i*100+j];
			  // finite difference solver with some parameters to capture diffusion
			  newgrid[i*100+j] = grid[i*100+j] + 2.5e2*(dt*(l + r - 2*h) + dt*(u + d - 2*h));
			  // loss of ATP. in proportion to current concentration; sites of loss depend on model
			  if(expt == 0 || expt == 1)
  			    loss = dt* newgrid[i*100+j]*kappa;
			  else
			    loss = ((i-50)*(i-50) + (j-50)*(j-50) < 25*25 ? dt* newgrid[i*100+j]*kappa : 0);
			  // impose nonnegativity (crudely and shouldn't happen
			  if(newgrid[i*100+j]-loss < 0) newgrid[i*100+j] = 0;
			  else { newgrid[i*100+j] -= loss; totalloss += loss; }
			}
		    }
		  // the above took care of loss and diffusion, now take care of production terms
		  // loop through the 100 mitos and add point ATP mass at each position
		  for(i = 0; i < 100; i++)
		    {
		      newgrid[x[i]*100+y[i]] += delta*1e7*dt; // consumption is 10^9 per cell per sec
		      totalgain += delta*1e7*dt;
		    }
		  // finally, update new state of cell from buffer
		  totalATP = 0; minATP = 1e100; maxATP = 0; 
		  for(i = 0; i < 100; i++)
		    {
		      for(j = 0; j < 100; j++)
			{
			  grid[i*100+j] = newgrid[i*100+j];
			  if(grid[i*100+j] < minATP) minATP = grid[i*100+j];
			  if(grid[i*100+j] > maxATP) maxATP = grid[i*100+j];
			  totalATP += grid[i*100+j];
			}
		    }
		  changeATP = totalATP - lastATP;
		  lastATP = totalATP;
	  
		  if((int)t != lastt)
		    {
		      // printf("Total loss in last second %.2e\n", totalloss/dt);
		      //printf("vol is %.2e dm3\n", vol);
		      //printf("Total ATP = %.2e molecules (%.2e mol) (conc %.2eM)\n", totalATP, totalATP/6e23, totalATP/6e23 / vol);

		      fp = fopen("stats.csv", "a");
		      fprintf(fp, "%i,%.2e,%.2e,%i,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,0\n", expt, kappa, delta, (int)t, totalATP, totalATP/6e23/vol, minATP, maxATP, changeATP, totalloss/dt);
		      fclose(fp);

		      if((int)t < 5 || (int)t % 100 == 0)
			{
			  sprintf(fstr, "out-%i-%.2f-%.2f-%i.txt", expt, kappa, delta, (int)t);
			  fp = fopen(fstr, "w");
			  for(i = 0; i < 100; i++)
			    {
			      for(j = 0; j < 100; j++)
				{
				  fprintf(fp, "%i %i %.5f\n", i, j, grid[i*100+j]);
				  total += grid[i*100+j];
				}
			      fprintf(fp, "\n");
			    }
			  fclose(fp);
			}
		      lastt = ((int)t);
		    }

		}
	      printf("stopped at %e\n", t);
	      fp = fopen("stats.csv", "a");
	      fprintf(fp, "%i,%.2e,%.2e,%i,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,1\n", expt, kappa, delta, (int)t, totalATP, totalATP/6e23/vol, minATP, maxATP, changeATP, totalloss/dt);
	      fclose(fp);

	      // final output details of mitos
	      sprintf(fstr, "mitos-%i.txt", expt);
	      fp = fopen(fstr, "w");
	      for(i = 0; i < 100; i++)
		fprintf(fp, "%i %i\n", x[i], y[i]);
	      fclose(fp);
	    }
	}
      //      printf("%.5f %.5f\n", totalloss/dt, totalgain/dt);
    }
  return 0;
}

// 3D
// volume = 100*1e-6 * 100*1e-6 * d*1e-6 = d * 1e-14 m^3
//        = d * 1e-11 dm-3
// generally about 2e11 molecules = 3e-13 mol
// -> 3e-13/(d * 1e-11) mol dm-3
// say d = 10um, then 3e-3 mol dm-3 = 3mM
