set multiplot
set size 0.33,0.33
set view map
unset key
set cbrange [0:1.3e5]

set origin 0,0.66
splot "out-0-4.00e-02-4.00e-01-2.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.33,0.66
splot "out-0-4.00e-02-4.00e-01-4.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.66,0.66
splot "out-0-4.00e-02-4.00e-01-100.txt" u 1:2:3 w pm3d, "mitos-0.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"

set origin 0,0.33
splot "out-1-4.00e-02-4.00e-01-2.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.33,0.33
splot "out-1-4.00e-02-4.00e-01-4.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.66,0.33
splot "out-1-4.00e-02-4.00e-01-100.txt" u 1:2:3 w pm3d, "mitos-1.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"

set origin 0,0.
splot "out-2-4.00e-02-4.00e-01-2.txt" u 1:2:3 w pm3d, "mitos-2.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.33,0.
splot "out-2-4.00e-02-4.00e-01-4.txt" u 1:2:3 w pm3d, "mitos-2.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"
set origin 0.66,0.
splot "out-2-4.00e-02-4.00e-01-100.txt" u 1:2:3 w pm3d, "mitos-2.txt" u 1:2:(103) w p pt 6 lc rgbcolor "#FFFFFF"


