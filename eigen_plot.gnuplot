set terminal x11 size 1600,800

# First Figure
set multiplot layout 4,1
set title "Eigen Plot 1"


datafile = "output/eigen_data.txt"
set grid
# Specify the range of columns for the first figure
plot datafile using 2:3 with linespoints pointtype 1 pointsize 1 title "Trajectory"
    # datafile using 8:9 with linespoints pointtype 1 pointsize 1 title "Ref Trajectory"


# Position the legend at the top-center for the first figure
set key top left




# Second Figure

set xlabel "X-Axis Label"
set ylabel "Y-Axis Label"

# Specify the range of columns for the second figure
plot datafile using 1:7 with linespoints pointtype 1 pointsize 1 title "DELTA"

# Position the legend at the top-center for the second figure
set key top left

# Third Figure

set xlabel "X-Axis Label"
set ylabel "Y-Axis Label"

# Specify the range of columns for the second figure
plot datafile using 1:5 with linespoints pointtype 1 pointsize 1 title "V1"
#     datafile using 1:10 with linespoints pointtype 1 pointsize 1 title "Ref v"


# Position the legend at the top-center for the second figure
set key top left



# Fourth Figure

set xlabel "X-Axis Label"
set ylabel "Y-Axis Label"

# Specify the range of columns for the second figure



# Position the legend at the top-center for the second figure
set key top left
unset multiplot
pause -1
