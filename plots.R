library(ggplot2)
library(ggpubr)

df.0 = read.csv("stats-0.csv")
df.1 = read.csv("stats-1.csv")
df.2 = read.csv("stats-2.csv")
df.3 = read.csv("stats-3.csv")
df.4 = read.csv("stats-4.csv")
df.5 = read.csv("stats-5.csv")
df.6 = read.csv("stats-6.csv")
df.7 = read.csv("stats-7.csv")
df = rbind(df.0, df.1, df.2, df.3,
           df.4, df.5, df.6, df.7)

expt.labels = c("Uni M, Uni C", "Uni M, Clu C", "Clu M, Uni C", "Clu M, Clu C",
                "Mo1 M, Uni C", "Mo1 M, Clu C", "Mo2 M, Uni C", "Mo2 M, Clu C")
df$expt.label = expt.labels[df$expt+1]

###### RQ1. When is a reasonable fold range in ATP supported at equilibrium?
# we are interested in whether bio-reasonable parameters can support high fold-changes in ATP conc

# compute fold change across cell, and report for those instances that equilibrated
df$fold.range = df$max.ATP/df$min.ATP
df$equilibrated = ifelse(df$terminated == 1 & df$t < 1000, 1, 0)

# did anything fail to equilibrate?
length(which(df$terminated == 1 & df$equilibrated != 1))

# timescales of equilibration
hist(log10(df$t[df$equilibrated == 1]))
max(df$fold.range[df$terminated == 1 & df$expt == 0])
max(df$fold.range[df$terminated == 1 & df$expt == 1])

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
  geom_rect(data = data.frame(xmin=5e8, xmax=5e9, ymin=5e9, ymax=5e11), 
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill = "lightblue") +
  geom_point(data = df[df$terminated==1 & df$expt<4,], 
             aes(x=consumption, y=total.ATP, size=log10(fold.range), color=log10(fold.range))) +
  scale_x_continuous(transform = "log10") + 
  scale_y_continuous(transform = "log10") + facet_wrap(~expt.label, scales = "free")

p.3
# so model 0 shows only limited maximal values; model 1 more; model 2 very high

p.3.a = ggplot() +
  geom_rect(data = data.frame(xmin=5e8, xmax=5e9, ymin=5e9, ymax=5e11), 
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

p.list = list()

# when do the individual experiments terminate?
set.snap = df[df$kappa==0.04 & df$delta==2.56 & df$terminated==1,]
terms = set.snap$t

for(expt in 1:4) {
  if(expt == 1) {
    snaps = c("out-0-0.04-2.56-1.txt", "out-0-0.04-2.56-2.txt", "out-0-0.04-2.56-4.txt", paste0("out-0-0.04-2.56-", terms[1], ".txt", collapse=""))
    snaps.m = rep("mitos-0.txt", 4)
  } else if(expt == 2) {
    snaps = c("out-1-0.04-2.56-1.txt", "out-1-0.04-2.56-2.txt", "out-1-0.04-2.56-4.txt", paste0("out-1-0.04-2.56-", terms[2], ".txt", collapse=""))
    snaps.m = rep("mitos-1.txt", 4)
  } else if(expt == 3) {
    snaps = c("out-2-0.04-2.56-1.txt", "out-2-0.04-2.56-2.txt", "out-2-0.04-2.56-4.txt", paste0("out-2-0.04-2.56-", terms[3], ".txt", collapse=""))
    snaps.m = rep("mitos-2.txt", 4)
  } else if(expt == 4) {
    snaps = c("out-3-0.04-2.56-1.txt", "out-3-0.04-2.56-2.txt", "out-3-0.04-2.56-4.txt", paste0("out-3-0.04-2.56-", terms[4], ".txt", collapse=""))
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
          ggarrange(plotlist=p.list[[4]], nrow=1),
          labels = paste(expt.labels, round(set.snap$fold.range, digits=2)), 
          nrow=4)


#

mt.df = read.csv("out-mitos-6-0.01-0.01.txt", sep=" ", header=FALSE)
colnames(mt.df) = c("t", "mito", "x", "y")
ggplot(mt.df[mt.df$t < 20,], aes(x=x, y=y, color=factor(mito))) + geom_path(size=0.5) + theme(legend.position="none")

mt.df = read.csv("out-mitos-7-0.08-0.01.txt", sep=" ", header=FALSE)
colnames(mt.df) = c("t", "mito", "x", "y")
ggplot() +
  geom_path(data=mt.df[mt.df$t < 100,], aes(x=x, y=y, color=factor(mito))) + 
  geom_point(data=mt.df[mt.df$t == max(mt.df$t),], aes(x=x, y=y, color=factor(mito))) + 
  theme(legend.position="none")
mt.df[mt.df$mito==1,]

atp.df = read.table("out-4-0.01-1.28.txt", sep=" ", header=FALSE)
colnames(atp.df) = c("t", "x", "y", "ATP")
mt.df = read.csv("out-mitos-4-0.01-1.28.txt", sep=" ", header=FALSE)
colnames(mt.df) = c("t", "mito", "x", "y")
ggplot() +
  geom_tile(data = atp.df[atp.df$t == max(atp.df$t),], aes(x=x, y=y, fill=ATP)) +
  geom_point(size=1,alpha=0.1,data=mt.df[mt.df$t < 1000,], aes(x=x, y=y, color=factor(mito))) + 
  geom_point(data=mt.df[mt.df$t == max(mt.df$t),], aes(x=x, y=y, color=factor(mito))) + 
  theme(legend.position="none")

dyn.list = list()
for(expt in 4:7) {
fname1 = paste0("out-", expt, "-0.04-2.56.txt", collapse="")
fname2 = paste0("out-mitos-", expt, "-0.04-2.56.txt", collapse="")
atp.df = read.table(fname1, sep=" ", header=FALSE)
colnames(atp.df) = c("t", "x", "y", "ATP")
mt.df = read.csv(fname2, sep=" ", header=FALSE)
colnames(mt.df) = c("t", "mito", "x", "y")
dyn.list[[expt-3]] = ggplot() +
  geom_tile(data = atp.df[atp.df$t == max(atp.df$t),], aes(x=x, y=y, fill=ATP)) +
  geom_path(size=1,alpha=0.1,data=mt.df[mt.df$t < 1000,], aes(x=x, y=y, group=factor(mito)), color="white") + 
  geom_point(data=mt.df[mt.df$t == max(mt.df$t),], aes(x=x, y=y), color="white") + 
  theme(legend.position="none")

}
ggarrange(plotlist = dyn.list, labels=c("A", "B", "C", "D"))

if(FALSE) {
  ggplot() +
  geom_tile(data = atp.df[atp.df$t == max(atp.df$t),], aes(x=x, y=y, fill=ATP)) +
  geom_path(size=1,alpha=0.1,data=mt.df[mt.df$t < 1000,], aes(x=x, y=y, color=factor(mito))) + 
  geom_point(data=mt.df[mt.df$t == max(mt.df$t),], aes(x=x, y=y, color=factor(mito))) + 
  theme(legend.position="none")
}

ggplot() +
  geom_tile(data = atp.df[atp.df$t == max(atp.df$t),], aes(x=x, y=y, fill=ATP)) +
  geom_path(size=1,alpha=0.1,data=mt.df[mt.df$t < 1000,], aes(x=x, y=y, group=factor(mito)), color="white") + 
  geom_point(data=mt.df[mt.df$t == max(mt.df$t),], aes(x=x, y=y), color="white") + 
  theme(legend.position="none")

ggplot() + geom_path(size=1,alpha=1,data=mt.df[mt.df$t < 0,], aes(x=x, y=y, color=factor(mito))) 

p.list[[expt]][[i]] = ggplot() +
  geom_tile(data=t.snap, aes(x=V1,y=V2,fill=V3)) +
  scale_fill_continuous(limits=c(0,NA)) +
  geom_point(data=t.snap.m, aes(x=V1,y=V2), color="white", size=0.5) +
  labs(x="", y="", fill="[ATP]")
