library(ggplot2)
library(ggpubr)

df.0 = read.csv("stats-0.csv")
df.1 = read.csv("stats-1.csv")
df.2 = read.csv("stats-2.csv")
df.3 = read.csv("stats-3.csv")
df = rbind(df.0, df.1, df.2, df.3)

###### RQ1. When is a reasonable fold range in ATP supported at equilibrium?
# we are interested in whether bio-reasonable parameters can support high fold-changes in ATP conc

# compute fold change across cell, and report for those instances that equilibrated
df$fold.range = df$max.ATP/df$min.ATP
df$equilibrated = ifelse(df$terminated == 1 & df$t < 1000, 1, 0)

# did anything fail to equilibrate?
length(which(df$terminated == 1 & df$equilibrated != 1))

# bio-reasonable values for consumption are around 10^9 ATP/cell/s; total ATP around 6e10 ATP/cell
df.legit = df[df$terminated==1 & 
     df$consumption > 5e8 & df$consumption < 5e9 & 
     df$total.ATP > 5e10 & df$fold.range > 1.,]
df.legit
# so low(er) ATP concentration is usually necessary for high fold range
# this typically gives higher consumption values, but not outrageously so

df[df$terminated==1 & 
     df$total.ATP > 1e11 & df$fold.range > 1.,]

order.df = df.legit[order(-df.legit$fold.range), ]
order.df = order.df[order.df$equilibrated == 1,]
head(order.df[order.df$expt==0,])
head(order.df[order.df$expt==1,])
head(order.df[order.df$expt==2,])

good.2 = df[df$expt==2 & df$equilibrated==1,]
good.2[good.2$consumption > 5e8 & good.2$consumption < 5e9 & good.2$total.ATP < 5e9,] 
# model 2 instances with higher fold ranges generally have lower ATP concentrations

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
  geom_rect(data = data.frame(xmin=5e8, xmax=5e9, ymin=5e9, ymax=5e11), 
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill = "lightblue") +
  geom_point(data = df[df$terminated==1,], 
             aes(x=consumption, y=total.ATP, size=log10(fold.range), color=log10(fold.range), label=kappa)) +
  scale_x_continuous(transform = "log10") + 
  scale_y_continuous(transform = "log10") + facet_wrap(~expt, scales = "free")

p.3
# so model 0 shows only limited maximal values; model 1 more; model 2 very high

# check for convergence
ggplot(df, aes(x=t, y=total.ATP, color=paste(kappa,delta))) + geom_point() + scale_x_continuous(transform = "log10")
ggplot(df, aes(x=t, y=consumption, color=paste(kappa,delta))) + geom_point() + scale_x_continuous(transform = "log10")
ggplot(df, aes(x=t, y=fold.range, color=paste(kappa,delta))) + geom_point() + 
  scale_x_continuous(transform = "log10") + scale_y_continuous(transform = "log10")

#### snapshots of particular instances

p.list = list()

# when do the individual experiments terminate?
terms = df$t[df$kappa==0.01 & df$delta==0.01 & df$terminated==1]

for(expt in 1:4) {
  if(expt == 1) {
    snaps = c("out-0-0.01-0.01-1.txt", "out-0-0.01-0.01-2.txt", "out-0-0.01-0.01-4.txt", paste0("out-0-0.01-0.01-", terms[1], ".txt", collapse=""))
    snaps.m = rep("mitos-0.txt", 4)
  } else if(expt == 2) {
    snaps = c("out-1-0.01-0.01-1.txt", "out-1-0.01-0.01-2.txt", "out-1-0.01-0.01-4.txt", paste0("out-1-0.01-0.01-", terms[2], ".txt", collapse=""))
    snaps.m = rep("mitos-1.txt", 4)
  } else if(expt == 3) {
    snaps = c("out-2-0.01-0.01-1.txt", "out-2-0.01-0.01-2.txt", "out-2-0.01-0.01-4.txt", paste0("out-2-0.01-0.01-", terms[3], ".txt", collapse=""))
    #snaps = c("out-2-0.08-5.12-1.txt", "out-2-0.08-5.12-4.txt", "out-2-2.56-0.01-100.txt", "out-2-2.56-0.01-400.txt")
    snaps.m = rep("mitos-2.txt", 4)
  } else if(expt == 3) {
    snaps = c("out-3-0.01-0.01-1.txt", "out-3-0.01-0.01-2.txt", "out-3-0.01-0.01-4.txt", paste0("out-3-0.01-0.01-", terms[4], ".txt", collapse=""))
    #snaps = c("out-2-0.08-5.12-1.txt", "out-2-0.08-5.12-4.txt", "out-2-2.56-0.01-100.txt", "out-2-2.56-0.01-400.txt")
    snaps.m = rep("mitos-3.txt", 4)
  }
  
  p.list[[expt]] = list()
  for(i in 1:length(snaps)) {
    t.snap = read.table(snaps[i], sep=" ", header=FALSE)
    t.snap.m = read.table(snaps.m[i], sep=" ", header=FALSE)
    p.list[[expt]][[i]] = ggplot() +
      geom_tile(data=t.snap, aes(x=V1,y=V2,fill=V3)) +
      scale_fill_continuous(limits=c(0,NA)) +
      geom_point(data=t.snap.m, aes(x=V1,y=V2), color="white", size=0.5) +
      labs(x="", y="", fill="[ATP]")
  }
}

ggarrange(ggarrange(plotlist=p.list[[1]], nrow=1),
          ggarrange(plotlist=p.list[[2]], nrow=1),
          ggarrange(plotlist=p.list[[3]], nrow=1), 
          ggarrange(plotlist=p.list[[4]], nrow=1), nrow=4)

