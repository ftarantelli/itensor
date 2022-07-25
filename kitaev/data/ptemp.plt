reset

set terminal postscript eps enhanced color dashed font 'Arial, 23'
set output 'graph.eps'

#-1-
# 0.5**6*1/2/(3/4)**(1/8)/6

#set terminal pngcairo font "Arial, 12"
#set terminal png font "Arial, 12"
#set output 'graph.png'

#set terminal 'jpeg' font arial 12 
#set output 'graph.jpeg'

#set for [i=1:9] linetype i dashtype i    {/Arial=20 u} = 0.01   {/Arial=20 uL^2} = 6000
#set style data linespoints

set multiplot layout 2,1 #rowsfirst title "{/:Bold=15 title}"\
              #margins screen MP_LEFT=0., MP_RIGHT=0., MP_BOTTOM=0., MP_TOP=0. spacing screen MP_xGAP=1., MP_yGAP=1.       # set xlabel after a plot to have not repetition

#set multiplot


array length[5] = [210, 400, 600, 800, 1000]
array dfile1[5]
array dfile2[5]

extr1 = -3
extr2 = 3
var = 2
#verse = 1
chr = 'w'
aytics = 0.2

if(ARG3) {extr1 = real(ARG3);}
if(ARG4) {extr2 = real(ARG4);}
if(ARG5) {var = int(ARG5);}
#if(ARG6) {verse = int(ARG6);}
if(ARG6) {if(int(ARG6) == 1) {chr = 'k';}}

do for [N=1:5] {
dfile1[N] = sprintf("wireABCu%s%s%sl%d.dat1", ARG1, chr, ARG2, length[N])
dfile2[N] = sprintf("wireABCu%s%s%sl%d.dat2", ARG1, chr, ARG2, length[N])
#wireABCu0w1l400.dat1
}

set xrange [extr1:extr2]

if(ARG7) {
	if(int(ARG7) == 1) {set key top right;};
	if(int(ARG7) == 2) {set key at graph 0.4, graph 0.97;};
	if(int(ARG7) == 3) {set key bottom left;};
	if(int(ARG7) == 4) {set key bottom right;};
}
set key font 'Arial, 30'

#set key opaque
#set key outside

#set lmargin at screen 0.15#{{at screen} <margin>} #set rmargin #set tmargin #set bmargin #set margins
set lmargin at screen 0.14
set bmargin at screen 0.73
set rmargin at screen 0.965
#set tmargin at screen 1

#set xtics axis                    # start, incr, end
set mxtics 4                                                # number of mtics in 1 interval
set xtics scale 2.3, 1.                                    # dimension major, minor tics   1=default
#set xtics border offset 0,0.5 -5,1,5
#set xtics offset 0,0.5

set mytics 4                                               
set ytics scale 2.3, 1.

#set  border lt 1 ls 2 lc 1 lw 2 dt 2

#set xtics 1
#set ytics 0.5

#set ytics add ("[+0.135]" 0.135) 
set xtics add ("0" 0) 
set ytics add ("0" 0)

set ytics font "Arial,20"
set xtics  font "Arial,0.01"


#set format x "%2.0t{/Symbol \264}10^{%L}"
#set format y "%2.0t{/Symbol \264}10^{%L}"
#set format x "10^{%L}"

#set format y "%.1f"   # number format in axis
#set format x "%1.0t {/Bold=18 ~ {.4\.}}10^{%L}"
#set format y "%1.0t {/Bold=18 ~ {.4\.}}10^{%L}"

set ylabel font 'Arial, 37'
set xlabel font 'Arial, 37'
#set format y "%.1f"   # number format in axis

#set grid 

#set title "Mean value of particles number N with {/Symbol=25 k}=1, {/Symbol=25 g}=1 and {/Symbol=25 d}=0"
#set title "Quench : {/Symbol=25 k} _i = 0 {/Symbol -->}  {/Symbol=25 k} = 3,   {/Symbol=25 g} = 100"
#set title " {/Arial=20 u} = 0.01 "
#set title "{/Symbol=25 k} = 0,  {/Symbol=25 g} = 20"#           {/Symbol=25 m}=1.8 and 
#set title "{/Symbol=25 d} = 1 and {/Symbol=25 m} = 2"
#set xrange [-3:3]
#set yrange [-1.8:0.8] 
#set xlabel "{/Symbol=45 q}" offset 0.0,0.5
#set xlabel "{/Symbol Q}" offset 0, 0.4
#set xlabel "{/Symbol=25 k}"
#set xlabel 'w' offset 0,0.8
if(var == 2) {set ylabel "L^{y_l} M^{/Roman (a)}" offset 1.4, 0;}	# L^{/Arial=20 y_{/Symbol f}} L^{y_{/Symbol \152}}
if(var == 3) {set ylabel "L C@_{L/3}^{/Roman (a)}" offset 1.4, 0;}
if(var == 4) {set ylabel "A^{/Roman (a)}" offset 1.4, 0;}
if(var == 5) {set ylabel "L E@^{/Roman (a)}_s" offset 1.4, 0;}

#set ylabel "{/Times=50 J}" offset 0.8,0
#set ylabel "L [ n_{L/2}({/Symbol=22 q}) - n_{L/2}(0) ]"
#set ylabel "L < C^+_{L/3} C_{2L/3} + h.c. >"
#set ylabel "~G{.6\\~}_c" offset 0.3, 0
#set ylabel "L G_c" offset 0.,0
#set ylabel "N / N_o" offset 0.8,0

#set logscale x 10
#set logscale y 10

#set arrow from -2, graph 0 to -2, graph 1 nohead lt 6 dt '_ '
#set arrow from graph 0, first 0.25 to graph 1, first 0.25 nohead lt 4 lw 2 dt '_._. ' 
#set arrow from graph 0, first -0.05  to graph 1, first -0.05 nohead lt 7 lw 1 dt '_._.'
#set arrow from graph 0, first -0.5  to graph 1, first -0.5 nohead lt 7 lw 1 dt '_._.'
#set arrow from graph 0, first -1.5  to graph 1, first -1.5 nohead lt 7 lw 1 dt '_._.'

set label font "Arial,30"

#set label 2 at first 0.01, first 0.55 "y = 0.5" tc lt 7 # graph, screen or first
#set label 3 at graph 0.93, first -0.6 "- 0.5"
#set label 4 at graph 0.955, first -1.1 "- 1"
#set label 1 at graph 0.05, first 0.45 "D@^{} _{}" tc lt 4
#set label 1 at graph 0.02, graph 0.09 "~{/Symbol m}{.9/Bold=35\\_} = 0" tc lt 8
#set label 2 at graph 0.85, graph 0.05 "{/Arial=30 w = 1}" tc lt 8


array expn[5] = [ 1., 1., 1., 1., 1. ]

#stats dfile2[1] using ($1):($2*6.**expn[var-1])  #nooutput #name 'ycol'
xttics = (extr2-extr1)/4.
if(ARG9){ xttics = real(ARG9);}
#yttics = (STATS_max_y - STATS_min_y)/4.
set xtics xttics


if(var == 2){
set ytics 100.
plot  dfile1[2] using ($1 ):($5*length[2]**expn[var-1]) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile1[5] using ($1 ):($5*length[5]**expn[var-1]) t sprintf("L = %d", length[5]) w l lw 5. lt 8 dt 1,\
}#      dfile2[6] using ($1 ):($2*16.**expn[var-1]) t 'L = 16' w l lw 5. lt 8 dt 1,\
#      dfile2[7] using ($1 ):($2*18.**expn[var-1]) t 'L = 18' w l lw 5. lt 6 dt 1,\
#}
if(var == 3){
set ytics 0.5
plot dfile1[2] using ($1 ):($2*length[2]**expn[var-1]) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile1[5] using ($1 ):($2*length[5]**expn[var-1]) t sprintf("L = %d", length[5]) w l lw 5. lt 8 dt 1,\
      }#      dfile2[7] using ($1 ):($3*18.**expn[var-1]) t 'L = 18' w l lw 5. lt 6 dt 1,\
#}
if(var == 4){
set ytics aytics
plot dfile1[2] using ($1 ):(exp(-$4)) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile1[5] using ($1 ):(exp(-$4)) t sprintf("L = %d", length[5]) w l lw 5. lt 8 dt 1,\
}

if(var == 5){
set ytics 400
plot dfile1[2] using ($1 ):($7*length[2]**expn[var-1]) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile1[5] using ($1 ):($7*length[5]**expn[var-1]) t sprintf("L = %d", length[5]) w l lw 3. lt 8 dt 1,\
      }
      

set bmargin at screen 0.129
set tmargin at screen 0.69

set xtics  font "Arial,20"
set xlabel "{/Symbol Q}" offset 0, 0.4

#set xrange [-extr2:-extr1]
#unset key


if(var == 2) {set ylabel "L^{y_l} M^{/Roman (b)}" offset 1.4, 0;}	# L^{/Arial=20 y_{/Symbol f}} L^{y_{/Symbol \152}}
if(var == 3) {set ylabel "L C@_{L/3}^{/Roman (b)}" offset 0.4, 0;}
if(var == 4) {set ylabel "A^{/Roman (b)}" offset 1.4, 0;}
if(var == 5) {set ylabel "L E@^{/Roman (b)}_s" offset 1.4, 0;}

set xtics  font "Arial,20"
set xlabel "{/Symbol Q}" offset 0, 0.4

#set key horizontal

if(ARG8) {
	if(int(ARG8) == 1) {
		set label 1 at graph 0.75,graph 0.9 "{/Arial=30 w_{/ZapfDingbats=15 \110} = 0.01}" tc lt 8;
		if(ARG6) {if(int(ARG6) == 1) {set label 1 at graph 0.75,graph 0.9 "{/Arial=30 {/Symbol Q}_{/ZapfDingbats=15 \110} = 10" tc lt 8;}};
		set label 2 at graph 0.75,graph 0.8 "{/Arial=30 {/Symbol U} = 0.001 }" tc lt 8;};

	if(int(ARG8) == 2) {
		set label 1 at graph 0.08, graph 0.9 "{/Arial=30 w_{/ZapfDingbats=15 \110} = 0.01}" tc lt 8;
		if(ARG6) {if(int(ARG6) == 1) {set label 1 at graph 0.08, graph 0.9 "{/Arial=30 {/Symbol Q}_{/ZapfDingbats=15 \110} = 10" tc lt 8;}};
		set label 2 at graph 0.08, graph 0.8 "{/Arial=30 {/Symbol U} = 0.001 }" tc lt 8;};

	if(int(ARG8) == 3) {
		set label 1 at graph 0.08, graph 0.15 "{/Arial=30 w_{/ZapfDingbats=15 \110} = 0.01}" tc lt 8;
		if(ARG6) {if(int(ARG6) == 1) {set label 1 at graph 0.08, graph 0.15 "{/Arial=30 {/Symbol Q}_{/ZapfDingbats=15 \110} = 10" tc lt 8;}};
		set label 2 at graph 0.08, graph 0.05 "{/Arial=30 {/Symbol U} = 0.001 }" tc lt 8;};

	if(int(ARG8) == 4) {
		set label 1 at graph 0.75,graph 0.15 "{/Arial=30 w_{/ZapfDingbats=15 \110} = 0.01}" tc lt 8;
		if(ARG6) {if(int(ARG6) == 1) {set label 1 at graph 0.75,graph 0.15 "{/Arial=30 {/Symbol Q}_{/ZapfDingbats=15 \110} = 10" tc lt 8;}};
		set label 2 at graph 0.75,graph 0.05 "{/Arial=30 {/Symbol U} = 0.001 }" tc lt 8;};


}

#set key outside



if(var == 2){
set ytics 100.
plot  dfile2[1] using ($1 ):($5*length[1]**expn[var-1]) t sprintf("L = %d  ", length[1]) w l lw 5. lt rgb 'dark-green' dt 4 ,\
      dfile2[2] using ($1 ):($5*length[2]**expn[var-1]) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile2[3] using ($1 ):($5*length[3]**expn[var-1]) t sprintf("L = %d  ", length[3]) w l lw 5. lt rgb 'red' dt '.',\
      dfile2[4] using ($1 ):($5*length[4]**expn[var-1]) t sprintf("L = %d  ", length[4]) w l lw 5. lt 1 dt 6,\
      dfile2[5] using ($1 ):($5*length[5]**expn[var-1]) t sprintf("L = %d", length[5]) w l lw 5. lt 8 dt 1,\
}#      dfile2[6] using ($1 ):($2*16.**expn[var-1]) t 'L = 16' w l lw 5. lt 8 dt 1,\
#      dfile2[7] using ($1 ):($2*18.**expn[var-1]) t 'L = 18' w l lw 5. lt 6 dt 1,\
#}
if(var == 3){
set ytics 1
set key bottom right
plot  dfile2[1] using ($1 ):($2*length[1]**expn[var-1]) t sprintf("L = %d  ", length[1]) w l lw 5. lt rgb 'dark-green' dt 4 ,\
      dfile2[2] using ($1 ):($2*length[2]**expn[var-1]) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile2[3] using ($1 ):($2*length[3]**expn[var-1]) t sprintf("L = %d  ", length[3]) w l lw 5. lt rgb 'red' dt '.',\
      dfile2[4] using ($1 ):($2*length[4]**expn[var-1]) t sprintf("L = %d  ", length[4]) w l lw 5. lt 1 dt 6,\
      dfile2[5] using ($1 ):($2*length[5]**expn[var-1]) t sprintf("L = %d", length[5]) w l lw 5. lt 8 dt 1,\
}#      dfile2[7] using ($1 ):($3*18.**expn[var-1]) t 'L = 18' w l lw 5. lt 6 dt 1,\
#}
if(var == 4){
set ytics aytics
plot  dfile2[1] using ($1 ):(exp(-$4)) t sprintf("L = %d  ", length[1]) w l lw 5. lt rgb 'dark-green' dt 4 ,\
      dfile2[2] using ($1 ):(exp(-$4)) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile2[3] using ($1 ):(exp(-$4)) t sprintf("L = %d  ", length[3]) w l lw 5. lt rgb 'red' dt '.',\
      dfile2[4] using ($1 ):(exp(-$4)) t sprintf("L = %d  ", length[4]) w l lw 5. lt 1 dt 6,\
      dfile2[5] using ($1 ):(exp(-$4)) t sprintf("L = %d", length[5]) w l lw 5. lt 8 dt 1,\
}#      dfile2[6] using ($1 ):(exp(-$4)) t 'L = 16' w l lw 5. lt 8 dt 1,\
#      dfile2[7] using ($1 ):($4*18.**expn[var-1]) t 'L = 18' w l lw 5. lt 6 dt 1,\
#} 

if(var == 5){
set key horizontal
set key top center
set ytics 400
plot  dfile2[1] using ($1 ):($7*length[1]**expn[var-1]) t sprintf("L = %d  ", length[1]) w l lw 5. lt rgb 'dark-green' dt 4 ,\
      dfile2[2] using ($1 ):($7*length[2]**expn[var-1]) t sprintf("L = %d  ", length[2]) w l lw 5. lt rgb 'blue' dt 5,\
      dfile2[3] using ($1 ):($7*length[3]**expn[var-1]) t sprintf("L = %d  ", length[3]) w l lw 5. lt rgb 'red' dt '.',\
      dfile2[4] using ($1 ):($7*length[4]**expn[var-1]) t sprintf("L = %d  ", length[4]) w l lw 5. lt 1 dt 6,\
      dfile2[5] using ($1 ):($7*length[5]**expn[var-1]) t sprintf("L = %d", length[5]) w l lw 3. lt 8 dt 1,\
}










unset multiplot
