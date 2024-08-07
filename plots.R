library(ggplot2)
library(ggpubr)

df = read.csv("stats.csv")

### when is a reasonable fold range supported at equilibrium?

df$fold.range = df$max.ATP/df$min.ATP
df[df$terminated==1,]

df[df$terminated==1 & df$expt == 2,]
df[df$terminated==1 & 
     df$consumption > 5e8 & df$consumption < 5e9 & 
     df$total.ATP > 5e10 & df$fold.range > 1.,]
# so low(er) ATP concentration is usually necessary for high fold range
# this typically gives higher consumption values, but not outrageously so

df[df$terminated==1 & 
     df$total.ATP > 1e11 & df$fold.range > 1.,]

dy = 0.25
# what combinations support fold range at reasonable ATP levels?
p.1 = ggplot(df[df$terminated==1,], aes(x=total.ATP, y=fold.range)) + 
  geom_point() + 
  geom_text(aes(x=total.ATP,y=fold.range,label=kappa), nudge_y=dy) + 
  geom_text(aes(x=total.ATP,y=fold.range,label=delta), nudge_y=-dy) +
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

p.3 = ggplot(df[df$terminated==1,], aes(x=consumption, y=total.ATP, size=log10(fold.range), label=kappa)) + 
  geom_point() + 
  scale_x_continuous(transform = "log10") + 
  scale_y_continuous(transform = "log10") + facet_wrap(~expt, scales = "free")

p.3

ggplot(df, aes(x=t, y=total.ATP, color=paste(kappa,delta))) + geom_point() + scale_x_continuous(transform = "log10")
ggplot(df, aes(x=t, y=consumption, color=paste(kappa,delta))) + geom_point() + scale_x_continuous(transform = "log10")
ggplot(df, aes(x=t, y=fold.range, color=paste(kappa,delta))) + geom_point() + 
  scale_x_continuous(transform = "log10") + scale_y_continuous(transform = "log10")

#### instances

snaps = c("out-0-0.01-0.01-1.txt", "out-0-0.01-0.01-4.txt", "out-0-0.01-0.01-100.txt", "out-0-0.01-0.01-400.txt")
snaps.m = rep("mitos-0.txt", 4)

snaps = c("out-1-0.01-0.01-1.txt", "out-1-0.01-0.01-4.txt", "out-1-0.01-0.01-100.txt", "out-1-0.01-0.01-400.txt")
snaps.m = rep("mitos-1.txt", 4)

snaps = c("out-2-0.01-0.01-1.txt", "out-2-0.01-0.01-4.txt", "out-2-0.01-0.01-100.txt", "out-2-0.01-0.01-400.txt")
snaps = c("out-2-0.08-5.12-1.txt", "out-2-0.08-5.12-4.txt", "out-2-2.56-0.01-100.txt", "out-2-2.56-0.01-400.txt")
snaps.m = rep("mitos-2.txt", 4)

p.list = list()
for(i in 1:length(snaps)) {
  t.snap = read.table(snaps[i], sep=" ", header=FALSE)
  t.snap.m = read.table(snaps.m[i], sep=" ", header=FALSE)
  p.list[[i]] = ggplot() +
    geom_tile(data=t.snap, aes(x=V1,y=V2,fill=V3)) +
    scale_fill_continuous(limits=c(0,NA)) +
    geom_point(data=t.snap.m, aes(x=V1,y=V2), color="white", size=0.5)
}
ggarrange(plotlist=p.list)
