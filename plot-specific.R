library(ggplot2)
library(ggpubr)
library(viridis)

gsize = 20
depth = 10

fname1.set = c("out-7-20-0.01-1.28-1.txt", "out-7-20-0.01-1.28-2.txt")
fname2.set = c("mitos-7-20-1.txt", "mitos-7-20-2.txt")
subdiv = c(1, 2)

g.plot = list()
for(i in 1:length(fname1.set)) {
  fname1 = fname1.set[i]
  fname2 = fname2.set[i]
  
  # ATP is output as molecule number in each simulation cell
  # so to scale up to dm^-3 we need (um3 volume associated with a simulation cell) * (um3 in 1 dm3)
  scale.sim.cell.um3 = 1/subdiv[i]
  scale.atp = (depth*scale.sim.cell.um3**2) * (1e-6 / 1e-1)**3
  scale.mol = 1./(6e23 * scale.atp) 

  atp.df = read.table(fname1, sep=" ", header=FALSE)
  colnames(atp.df) = c("t", "x", "y", "ATP")
  mt.df = read.csv(fname2, sep=" ", header=FALSE)
  colnames(mt.df) = c("x", "y")
  
  t.end = max(atp.df$t)
  this.df = atp.df[atp.df$t == t.end,]
  this.df$col = this.df$ATP*scale.mol*1e3
  g.plot[[i]] = ggplot() +
    geom_tile(data = this.df, aes(x=x, y=y, fill=col)) +
    scale_fill_viridis_c(option = "magma", limits = range(this.df$col) + c(-0.01, 0.01)) +
    labs(fill="[ATP]/mM") +     
    geom_point(data=mt.df, aes(x=x, y=y), color="white")
}

ggarrange(plotlist = g.plot)
