set multiplot
set size 1,0.5
set view map
unset key
set origin 0,0.5
splot "out-0-1.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0,0
splot "out-1-2.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
