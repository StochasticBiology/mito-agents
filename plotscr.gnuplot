set multiplot
set size 0.33,0.5
set view map
unset key
set cbrange [0:1.3e7]

set origin 0,0.5
splot "out-0-5.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.33,0.5
splot "out-0-200.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.66,0.5
splot "out-0-1000.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"

set origin 0,0.
splot "out-1-5.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.33,0.
splot "out-1-200.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.66,0.
splot "out-1-1000.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"


