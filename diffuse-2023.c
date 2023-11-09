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
  double dt = 0.00001;
  int lastt;
  double total;
  int expt;

  // allocate memory to store 100x100 grid of ATP concentration values
  grid = (double*)malloc(sizeof(double)*100*100);
  newgrid = (double*)malloc(sizeof(double)*100*100);

  // two different experiments: random mito spread or clustering near centre
  for(expt = 0; expt <= 1; expt++)
    {
      // initialises positions of mitochondria (x, y are coordinates)
      for(i = 0; i < 100; i++)
	{
	  if(expt == 0) { x[i] = RND*100; y[i] = RND*100; }
	  else { x[i] = RND*20+40; y[i] = RND*20+40; }
	  for(j = 0; j < 100; j++)
	    grid[i*100+j] = 0;
	}

      // discrete time PDE solver (very crude)
      lastt = -1;
      for(t = 0; t < 2.1; t+=dt)
	{
	  // output snapshot of simulation periodically
	  if((int)t != lastt)
	    {
	      total = 0;
	      sprintf(fstr, "out-%i-%i.txt", expt, (int)t);
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
	      printf("%i %.5f\n", (int)t, total);
	      lastt = (int)t;
	    }
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
		  newgrid[i*100+j] = grid[i*100+j] + 5*5*1.5e2*(dt*(l + r - 2*h) + dt*(u + d - 2*h));
		  // uniform loss of ATP everywhere in cell
		  loss = dt*1e9/(100.*100.);
		  // impose nonnegativity (crudely)
		  if(newgrid[i*100+j]-loss < 0) newgrid[i*100+j] = 0;
		  else { newgrid[i*100+j] -= loss; totalloss += loss; }
		}
	    }
	  // the above took care of loss and diffusion, now take care of production terms
	  // loop through the 100 mitos and add point ATP mass at each position
	  for(i = 0; i < 100; i++)
	    {
	      newgrid[x[i]*100+y[i]] += 1e7*dt; // consumption is 10^9 per cell per sec
	      totalgain += 1e7*dt;
	    }
	  // finally, update new state of cell from buffer
	  for(i = 0; i < 100; i++)
	    {
	      for(j = 0; j < 100; j++)
		grid[i*100+j] = newgrid[i*100+j];
	    }
	}
      // final output details of mitos
      sprintf(fstr, "mitos-%i.txt", expt);
      fp = fopen(fstr, "w");
      for(i = 0; i < 100; i++)
	fprintf(fp, "%i %i\n", x[i], y[i]);
      fclose(fp);

      printf("%.5f %.5f\n", totalloss/dt, totalgain/dt);
    }
  return 0;
}
