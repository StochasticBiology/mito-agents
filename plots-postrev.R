library(ggplot2)
library(ggpubr)
library(viridis)
library(dplyr)

expts.set = list( c(50, 10, 100, 2) )

for(this.expts in expts.set) {
  gsize = this.expts[1]
  depth = this.expts[2]
  nmito = this.expts[3]
  subdiv = this.expts[4] # number of simulation cells in 1um3
  
  flabels = c(0:9, -1, -3, -5, -7, -9)
  
  df = data.frame()
  fnames = paste0("stats-", flabels, "-", gsize, "-", depth, "-", nmito, "-", subdiv, ".csv")
  for(fname in fnames) {
    tmp.df = read.csv(fname)
    df = rbind(df, tmp.df)
  }
  df$expt[df$expt < 0] = 10 + (-df$expt[df$expt < 0]-1)/2
  
  print(paste0("----- Thinking about ", paste0(this.expts, collapse=" ")))
  
  expt.labels = c("Uniform mitos, Uniform κ", "Uniform mitos, Local κ", "Local mitos, Uniform κ", "Local mitos, Local κ",
                  "Diffusive* mitos, Uniform κ", "Diffusive* mitos, Local κ", "Directed mitos, Uniform κ", "Directed mitos, Local κ",
                  "Diffusive mitos, Uniform κ", "Diffusive mitos, Local κ",
                  "Uniform mitos, Fibrous κ", "Local mitos, Fibrous κ", "Diffusive* mitos, Fibrous κ", "Directed mitos, Fibrous κ", "Diffusive mitos, Fibrous κ")
  df$expt.label = expt.labels[df$expt+1]
  
  ###### RQ1. When is a reasonable fold range in ATP supported at equilibrium?
  # we are interested in whether bio-reasonable parameters can support high fold-changes in ATP conc
  
  # compute fold change across cell, and report for those instances that equilibrated
  # total.ATP is the total number of ATP molecules in the whole simulated cell
  # vol.dm3 gives the cellular volume in dm3
  # hence conc.ATP here will give the cellular ATP concentration in mol dm-3
  df$fold.range = df$max.ATP/df$min.ATP
  df$conc.ATP = df$total.ATP/df$vol.dm3 / 6e23
  df$equilibrated = ifelse(df$terminated == 1 & df$t < 1000, 1, 0)
  df$CV = sqrt(df$var.ATP)/df$mean.ATP
  
  
  ggplot(df[df$terminated == 1,], aes(x=kappa, y=delta, fill=log10(conc.ATP))) + 
    scale_x_continuous(transform="log2") + 
    scale_y_continuous(transform="log2") + 
    geom_tile() + facet_wrap(~expt)
  
  
  
  # did anything fail to equilibrate?
  length(which(df$terminated == 1 & df$equilibrated != 1))
  
  # timescales of equilibration
  hist(log10(df$t[df$equilibrated == 1]))
  max(df$fold.range[df$terminated == 1 & df$expt == 0])
  max(df$fold.range[df$terminated == 1 & df$expt == 1])
  
  # bio-reasonable values for consumption are around 10^9 ATP/cell/s; total ATP around 6e10 ATP/cell
  df.legit = df[df$terminated==1 & 
                  df$consumption > 1e8 & df$consumption < 1e10 & 
                  df$conc.ATP > 5e-4,]
  df.legit
  # so low(er) ATP concentration is usually necessary for high fold range
  # this typically gives higher consumption values, but not outrageously so
  
  df[df$terminated==1 & 
       df$conc.ATP > 5e-4,]
  
  order.df = df.legit[order(-df.legit$fold.range), ]
  order.df = order.df[order.df$equilibrated == 1,]
  head(order.df[order.df$expt==0,])
  head(order.df[order.df$expt==1,])
  head(order.df[order.df$expt==2,])
  head(order.df[order.df$expt==3,])
  
  good.2 = df[df$expt==2 & df$equilibrated==1,]
  good.2 = good.2[order(-good.2$fold.range),]
  head(good.2[good.2$consumption > 5e8 & good.2$consumption < 5e9 & good.2$total.ATP < 5e9,] )
  good.3 = df[df$expt==3 & df$equilibrated==1,]
  good.3 = good.3[order(-good.3$fold.range),]
  head(good.3[good.3$consumption > 5e8 & good.3$consumption < 5e9 & good.3$total.ATP < 5e9,])
  # model 2 and 3 instances with higher fold ranges generally have lower ATP concentrations
  
  df.legit[df.legit$expt==2,]
  # these plots aren't very helpful
  if(FALSE) {
    # offset for plotting
    dy = 0.25
    # what combinations support fold range at reasonable ATP levels?
    p.1 = ggplot(df[df$terminated==1,], aes(x=total.ATP, y=fold.range)) + 
      geom_point(size=1) + 
      geom_text(aes(x=total.ATP,y=fold.range,label=kappa), nudge_y=dy, size=2) + 
      geom_text(aes(x=total.ATP,y=fold.range,label=delta), nudge_y=-dy, size=2) +
      scale_x_continuous(transform = "log10") + facet_wrap(~expt, scales = "free")
    
    p.2 = ggplot(df[df$terminated==1,], aes(x=consumption, y=fold.range)) + 
      geom_point() + 
      geom_text(aes(x=consumption,y=fold.range,label=kappa), nudge_y=dy) + 
      geom_text(aes(x=consumption,y=fold.range,label=delta), nudge_y=-dy) +
      scale_x_continuous(transform = "log10") + facet_wrap(~expt, scales = "free")
    
    # higher kappa gives higher fold range
    # but need to balance with delta to keep consumption and concentration reasonable
    # ATP 1e11 with a 10um-thick cell gives ~1mM
    # ATP 1e10 with a 1um-thick cell gives ~1mM
    ggarrange(p.1, p.2, nrow=2)
  }
  
  p.3 = ggplot() +
    geom_rect(data = data.frame(xmin=5e8, xmax=5e9, ymin=5e-4, ymax=1e-2), 
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill = "lightblue") +
    geom_point(data = df[df$terminated==1 & df$expt<4,], 
               aes(x=consumption, y=conc.ATP, size=log10(fold.range), color=log10(fold.range))) +
    scale_x_continuous(transform = "log10") + 
    scale_y_continuous(transform = "log10") + facet_wrap(~expt.label, scales = "free")
  
  p.3 + scale_color_viridis()
  
  p.3.zoom = ggplot() +
    geom_point(data = df[df$terminated==1 & 
                           df$conc.ATP > 5e-4 & df$conc.ATP < 1e-2 &
                           df$consumption > 1e8,], 
               aes(x=consumption, y=conc.ATP, size=fold.range, color=fold.range)) +
    scale_x_continuous(transform = "log10") + 
    scale_y_continuous(transform = "log10") + facet_wrap(~expt.label, scales = "free")
  
  p.3.zoom + scale_color_viridis()
  
   plot.order = c("Uniform mitos, Uniform κ",  "Local mitos, Uniform κ", "Diffusive mitos, Uniform κ", "Diffusive* mitos, Uniform κ", "Directed mitos, Uniform κ", 
                 "Uniform mitos, Local κ", "Local mitos, Local κ", "Diffusive mitos, Local κ", "Diffusive* mitos, Local κ", "Directed mitos, Local κ",
                 "Uniform mitos, Fibrous κ", "Local mitos, Fibrous κ", "Diffusive mitos, Fibrous κ", "Diffusive* mitos, Fibrous κ", "Directed mitos, Fibrous κ")
  
   df$expt.label.red = gsub(" mitos", "", df$expt.label)
   df$expt.label.red = gsub("Local", "Central", df$expt.label.red)
   plot.order.red = gsub(" mitos", "", plot.order)
   plot.order.red = gsub("Local", "Central", plot.order.red)
  p.2.zoom.a = ggplot() +
    geom_point(data = df[df$terminated==1 & 
                           df$conc.ATP > 5e-4 & df$conc.ATP < 1e-2 &
                           df$consumption > 1e8,], 
               aes(x=conc.ATP*1e3, y=consumption, size=fold.range, fill=fold.range), pch=21) +
    scale_x_continuous(transform = "log10") + 
    scale_y_continuous(transform = "log10") + 
    scale_fill_viridis(option="inferno") +
    labs(y="ATP consumption / "~cell^{-1}~s^{-1}, x="[ATP] / mM", fill="Fold range", size="Fold range") +
    facet_wrap(~factor(expt.label.red, levels=plot.order.red), nrow=3, ncol=5) +
    theme_light()
  
  p.3.zoom.a = ggplot() +
    geom_point(data = df[df$terminated==1 & 
                           df$conc.ATP > 5e-4 & df$conc.ATP < 1e-2 &
                           df$consumption > 1e8,], 
               aes(x=conc.ATP*1e3, y=consumption, size=CV, fill=CV), pch=21) +
    scale_x_continuous(transform = "log10") + 
    scale_y_continuous(transform = "log10") + 
    scale_fill_viridis(option="inferno") +
    labs(y="ATP consumption / "~cell^{-1}~s^{-1}, x="[ATP] / mM") +
    facet_wrap(~factor(expt.label.red, levels=plot.order.red), nrow=3, ncol=5) +
    theme_light()
  
  # so model 0 shows only limited maximal values; model 1 more; model 2 very high
  
  ggplot(df.legit, aes(x=log10(CV), fill=factor(expt.label))) + geom_histogram(position="dodge")
  ggplot(df.legit, aes(x=log10(CV), fill=factor(expt.label))) + geom_density(alpha=0.3) + facet_wrap(~expt.label)
  
  p.3.cv = ggplot() +
    geom_rect(data = data.frame(xmin=1e8, xmax=1e10, ymin=5e-4, ymax=1e-2), 
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill = "lightblue") +
    geom_point(data = df[df$terminated==1 & df$expt<4,], 
               aes(x=consumption, y=conc.ATP, size=CV, color=CV)) +
    scale_x_continuous(transform = "log10") + 
    scale_y_continuous(transform = "log10") + facet_wrap(~expt.label, scales = "free")
  
  tmp.df = df[df$terminated==1 & df$expt<4,]
  p.3.cv + scale_color_gradientn(
    colors = c("black", "blue", "white", "red"),
    values = scales::rescale(c(min(tmp.df$CV), 0.5, 0.5, max(tmp.df$CV)))
  ) 
  
  p.3.a = ggplot() +
    geom_rect(data = data.frame(xmin=1e8, xmax=1e10, ymin=5e9, ymax=5e11), 
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill = "lightblue") +
    geom_point(data = df[df$terminated==1 & df$expt>=4,], 
               aes(x=consumption, y=total.ATP, size=log10(fold.range), color=log10(fold.range))) +
    scale_x_continuous(transform = "log10") + 
    scale_y_continuous(transform = "log10") + facet_wrap(~expt.label, scales = "free")
  
  p.3.a
  
  # check for convergence
  ggplot(df, aes(x=t, y=total.ATP, color=paste(kappa,delta))) + geom_point() + scale_x_continuous(transform = "log10")
  ggplot(df, aes(x=t, y=consumption, color=paste(kappa,delta))) + geom_point() + scale_x_continuous(transform = "log10")
  ggplot(df, aes(x=t, y=fold.range, color=paste(kappa,delta))) + geom_point() + 
    scale_x_continuous(transform = "log10") + scale_y_continuous(transform = "log10")
  
  #### snapshots of particular instances
  
  # which ones are most comparable and reasonable
  df[df$terminated==1 & df$conc.ATP > 2e-3 & df$conc.ATP < 4e-3,]
  
  nexpt = length(unique(df$expt.label))
  if(gsize == 50) {
    kappas = rep(0.01, nexpt)
    deltas = rep(0.64, nexpt)
    kappas = rep(0.16, nexpt)
    deltas = rep(5.12, nexpt)
    #was deltas = rep(c(0.32, 5.12), 4)
  }
  if(gsize == 20) {
    kappas = rep(0.01, nexpt)
    deltas = rep(0.08, nexpt)
    kappas = rep(0.64, nexpt)
    deltas = rep(2.56, nexpt)
    # was  deltas = rep(c(0.08, 0.32), 4)
  }
  
  expt = 5
  fname2 = paste0("out-mitos-", expt-1, "-", gsize, "-", depth, "-", nmito, "-", kappas[expt], "-", deltas[expt], "-", subdiv, ".txt", collapse="")
  mt.df = read.csv(fname2, sep=" ", header=FALSE)
  colnames(mt.df) = c("t", "mito", "x", "y")
  mt.df[mt.df$mito==0,]
  dyn.list = list()
  
  # ATP is output as molecule number in each simulation cell
  # so to scale up to dm^-3 we need (um3 volume associated with a simulation cell) * (um3 in 1 dm3)
  scale.sim.cell.um3 = 1/subdiv
  scale.atp = (depth*scale.sim.cell.um3**2) * (1e-6 / 1e-1)**3
  scale.mol = 1./(6e23 * scale.atp) 
  
  eq.df = df[df$terminated==1 & df$conc.ATP > 5e-4,]
  mito.speed.dist = list()
  for(expt in 1:nexpt) {
    tmp.sub = df[df$expt == expt-1 & df$kappa == kappas[expt] & df$delta == deltas[expt],]
    final.sub = tmp.sub[tmp.sub$t == max(tmp.sub$t),]
    print(paste0("Expt ", expt, " [", flabels[expt], "] consumption ", final.sub$consumption[1]))
    
    dyn.list[[expt]] = list()
    fname1 = paste0("out-", flabels[expt], "-", gsize, "-", depth, "-", nmito, "-", kappas[expt], "-", deltas[expt], "-", subdiv, ".txt", collapse="")
    
    atp.df = read.table(fname1, sep=" ", header=FALSE)
    colnames(atp.df) = c("t", "x", "y", "ATP")

    if(!(expt %in% c(1:4, 11, 12))) {
      fname2 = paste0("out-mitos-", flabels[expt], "-", gsize, "-", depth, "-", nmito, "-", kappas[expt], "-", deltas[expt], "-", subdiv, ".txt", collapse="")
      mt.df = read.csv(fname2, sep=" ", header=FALSE)
      colnames(mt.df) = c("t", "mito", "x", "y")
    } else {
      fname2 = paste0("mitos-", flabels[expt], "-", gsize, "-", depth, "-", nmito, "-", subdiv, ".txt", collapse="")
      mt.df = read.csv(fname2, sep=" ", header=FALSE)
      colnames(mt.df) = c("x", "y")
    }
    t.set = c(max(atp.df$t))
    for(i in 1:length(t.set)) {
      this.df = atp.df[atp.df$t == t.set[i],]
      this.df$col = this.df$ATP*scale.mol*1e3
      dyn.list[[expt]][[i]] = ggplot() +
        geom_tile(data = this.df, aes(x=x, y=y, fill=col)) +
        scale_fill_viridis_c(option = "magma", limits = range(this.df$col) + c(-0.01, 0.01)) +
        labs(fill="[ATP]/mM") #+
      #ggtitle(paste0("    ", t.set[i], collapse=""))
      if(!(expt %in% c(1:4, 11, 12))) {
        if(expt == 5 | expt == 6 || expt == 13) { 
          plot.period = 10 
        } else {
          plot.period = Inf
        }
        dyn.list[[expt]][[i]] = dyn.list[[expt]][[i]] +     
          geom_path(size=0.2,alpha=0.5,data=mt.df[mt.df$t <= t.set[i] & mt.df$t > t.set[i]-plot.period,], 
                    aes(x=x, y=y, group=factor(mito)), color="white") + 
          geom_point(data=mt.df[mt.df$t == t.set[i],], aes(x=x, y=y), 
                     shape = 21, fill = NA, color = "white", size = 0.8, stroke = 0.5)
        
        df_with_speed <- mt.df %>%
          arrange(mito, t) %>%  # Ensure data is sorted by mito and time
          group_by(mito) %>%    # Group by individual mitochondrion
          mutate(
            dx = x - lag(x),    # Change in x
            dy = y - lag(y),    # Change in y
            dt = t - lag(t),    # Time difference
            speed = sqrt(dx^2 + dy^2) / dt  # Calculate speed
          ) %>%
          ungroup()
        mito.speed.dist[[expt]] = df_with_speed$speed
      } else {
        dyn.list[[expt]][[i]] = dyn.list[[expt]][[i]] +     
          geom_point(data=mt.df, aes(x=x, y=y), 
                     shape = 21, fill = NA, color = "white", size = 0.8, stroke = 0.5)
      }
    }
  }
  
  g.speeds = ggarrange(ggplot(data.frame(x=mito.speed.dist[[5]]), aes(x=x)) + geom_histogram(),
                       ggplot(data.frame(x=mito.speed.dist[[6]]), aes(x=x)) + geom_histogram(),
                       ggplot(data.frame(x=mito.speed.dist[[9]]), aes(x=x)) + geom_histogram(),
                       ggplot(data.frame(x=mito.speed.dist[[10]]), aes(x=x)) + geom_histogram(),
                       ggplot(data.frame(x=mito.speed.dist[[7]]), aes(x=x)) + geom_histogram(),
                       ggplot(data.frame(x=mito.speed.dist[[8]]), aes(x=x)) + geom_histogram(),
                       labels=c("A", "B", "C", "D", "E", "F"))
  
  g.static = ggarrange(
    ggarrange(plotlist = dyn.list[[1]], nrow=1),
    ggarrange(plotlist = dyn.list[[2]], nrow=1),
    ggarrange(plotlist = dyn.list[[11]], nrow=1),
    ggarrange(plotlist = dyn.list[[3]], nrow=1),
    ggarrange(plotlist = dyn.list[[4]], nrow=1),
    ggarrange(plotlist = dyn.list[[12]], nrow=1),
    nrow=2, ncol=3
  )
  
  g.static.labs = ggarrange(
   ggarrange(plotlist = dyn.list[[1]], nrow=1),
    ggarrange(plotlist = dyn.list[[2]], nrow=1),
    ggarrange(plotlist = dyn.list[[11]], nrow=1),
    ggarrange(plotlist = dyn.list[[3]], nrow=1),
    ggarrange(plotlist = dyn.list[[4]], nrow=1),
    ggarrange(plotlist = dyn.list[[12]], nrow=1),
    nrow=2, ncol=3,
   labels=c("Ai", "ii", "iii", "Bi", "ii", "iii"),
   label.x = -0.02
  )
  
  # expt.labels
  
  g.dynamic = ggarrange(
    ggarrange(plotlist = dyn.list[[5]], nrow=1),
    ggarrange(plotlist = dyn.list[[9]], nrow=1),
    ggarrange(plotlist = dyn.list[[7]], nrow=1),
    
    ggarrange(plotlist = dyn.list[[6]], nrow=1),
    ggarrange(plotlist = dyn.list[[10]], nrow=1),
    ggarrange(plotlist = dyn.list[[8]], nrow=1),
    
    ggarrange(plotlist = dyn.list[[13]], nrow=1),
    ggarrange(plotlist = dyn.list[[15]], nrow=1),
    ggarrange(plotlist = dyn.list[[14]], nrow=1),
    nrow=3, ncol=3
  )
  
  g.dynamic.labs = ggarrange(
    ggarrange(plotlist = dyn.list[[9]], nrow=1),
    ggarrange(plotlist = dyn.list[[10]], nrow=1),
    ggarrange(plotlist = dyn.list[[15]], nrow=1),
    
    ggarrange(plotlist = dyn.list[[5]], nrow=1),
    ggarrange(plotlist = dyn.list[[6]], nrow=1),
    ggarrange(plotlist = dyn.list[[13]], nrow=1),
    
    ggarrange(plotlist = dyn.list[[7]], nrow=1),
    ggarrange(plotlist = dyn.list[[8]], nrow=1),
    ggarrange(plotlist = dyn.list[[14]], nrow=1),
    nrow=3, ncol=3,
    labels=c("Ai", "ii", "iii", "Bi", "ii", "iii", "Ci", "ii", "iii"),
    label.x = -0.02
  )
  
  #dyn.list[[1]][[1]] + scale_fill_viridis_c(limits = range(data$value) + c(-0.05, 0.05))
  ######
  sf = 2
  
  fname = paste0("speeds-", gsize, "-", depth, "-", nmito, "-", subdiv, ".png", collapse="")
  png(fname, width=600*sf, height=400*sf, res=72*sf)
  print(g.speeds)
  dev.off()
  
  fname = paste0("cv-zoom-", gsize, "-", depth, "-", nmito, "-", subdiv, ".png", collapse="")
  png(fname, width=600*sf, height=400*sf, res=72*sf)
  print(p.3.zoom.a)
  dev.off()
  
  fname = paste0("fr-zoom-", gsize, "-", depth, "-", nmito, "-", subdiv, ".png", collapse="")
  png(fname, width=600*sf, height=400*sf, res=72*sf)
  print(p.2.zoom.a + scale_color_viridis())
  dev.off()
  
  fname = paste0("snaps-static-", gsize, "-", depth, "-", nmito, "-", subdiv, ".png", collapse="")
  png(fname, width=900*sf, height=400*sf, res=72*sf)
  print(g.static.labs)
  dev.off()
  
  fname = paste0("snaps-dynamic-", gsize, "-", depth, "-", nmito, "-", subdiv, ".png", collapse="")
  png(fname, width=800*sf, height=600*sf, res=72*sf)
  print(g.dynamic.labs)
  dev.off()
  
  fname = paste0("all-", gsize, "-", depth, "-", nmito, "-", subdiv, ".png", collapse="")
  png(fname, width=900*sf, height=600*sf, res=72*sf)
  print(ggarrange(
    p.3.zoom.a + scale_color_viridis(),
    ggarrange(g.static, g.dynamic, nrow=1, labels=c("B", "C")),
    nrow=2, 
    labels = c("A", "")
  ))
  dev.off()
}
