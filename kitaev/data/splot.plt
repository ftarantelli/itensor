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

#set for [i=1:9] linetype i dashtype i    {/Arial=20 u} =  0.01   {/Arial=20 uL^2} =  6000
#set style data linespoints

#set multiplot layout 2,2 rowsfirst title "{/:Bold=15 title}"\
              #margins screen MP_LEFT=0., MP_RIGHT=0., MP_BOTTOM=0., MP_TOP=0. spacing screen MP_xGAP=1., MP_yGAP=1.       # set xlabel after a plot to have not repetition

set multiplot

set key spacing 1.2
set key top left          # top, bottom, right, left, center
#set key at first 7.7, first 0.3
#set key width 0.1 height 0.1
#set key horizontal 
#set key maxcols 4 #maxrows
#set key samplen 3        # length of line's type
set key font 'Arial, 30'

#set lmargin at screen 0.15#{{at screen} <margin>} #set rmargin #set tmargin #set bmargin #set margins
set lmargin at screen 0.14
set bmargin at screen 0.129
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

set ytics 400
#set ytics 1

#set ytics add ("[+0.135]" 0.135) 
set xtics add ("0" 0) 
set ytics add ("0" 0)

set ytics font "Arial,20"
set xtics  font "Arial,20"

#set format x "%2.0t{/Symbol \264}10^{%L}"
#set format y "%2.0t{/Symbol \264}10^{%L}"
#set format x "10^{%L}"

#set format y "%.1f"   # number format in axis
#set format x "%1.0t {/Bold=18 ~ {.4\.}}10^{%L}"
#set format y "%1.0t {/Bold=18 ~ {.4\.}}10^{%L}"

set ylabel font 'Arial, 35'
set xlabel font 'Arial, 35'
#set format y "%.1f"   # number format in axis

#set grid 

#set title "Mean value of particles number N with {/Symbol=25 k}=1, {/Symbol=25 g}=1 and {/Symbol=25 d}=0"
#set title "Quench : {/Symbol=25 k} _i = 0 {/Symbol -->}  {/Symbol=25 k} =  3,   {/Symbol=25 g} =  100"
#set title " {/Arial=20 u} =  0.01 "
#set title "{/Symbol=25 k} =  0,  {/Symbol=25 g} =  20"#           {/Symbol=25 m}=1.8 and 
#set title "{/Symbol=25 d} =  1 and {/Symbol=25 m} =  2"
set xrange [-140:140]
#set yrange [-150:] 
#set xlabel "{/Symbol=45 q}" offset 0.0,0.5
set xlabel "{/Symbol Q}" offset 0, 0.4
#set xlabel "{/Symbol=25 k}"
#set xlabel 'w' offset 0,0.8
set ylabel "L ({/Symbol=40 r - r}_c)" offset 1.5, 0		# L^{/Arial=20 y_{/Symbol f}}
#set ylabel "A" offset 2.,0
#set ylabel "L [ n_{L/2}({/Symbol=22 q}) - n_{L/2}(0) ]"
#set ylabel "L < C^+_{L/3} C_{2L/3} + h.c. >"
#set ylabel "~G{.6\\~}_c" offset 0.3, 0
#set ylabel "L G_c" offset 0.,0
#set ylabel "N / N_o" offset 0.8,0

#set logscale x 10
#set logscale y 10

#set arrow from -0.7532999997, graph 0 to -0.7532999997,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from -0.5062499994, graph 0 to -0.5062499994,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from -0.2565499992, graph 0 to -0.2565499992,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from -0.0069999989, graph 0 to -0.0069999989,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from 0.24285000084, graph 0 to 0.24285000084,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from 0.49255000065, graph 0 to 0.49255000065,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from 0.74040000047, graph 0 to 0.74040000047,graph 1 nohead lt 8 dt 2 lw 3
#set arrow from graph 0, first -1.5  to graph 1, first -1.5 nohead lt 7 lw 1 dt '_._.'

set label font "Arial,30"

#set label 1 at graph 0.82, graph 0.88 "{/Arial=30 L = 400}" tc lt 8

set label 3 at graph 0.75, graph 0.08 "{/Arial=30 L = 2000}" tc lt 8
set label 2 at graph 0.75, graph 0.18 "{/Arial=30 {/Symbol U} =  0.001}" tc lt 8

#plot  "wireABCu10k1l50.dat1" using ($1):($5*2000.0*50.**(0.)) t 'L = 50  ' w l lw 5. lt 8 dt 7,\

#plot   "wireABCu0k200l200.dat1" using ($1):($5*2000.0*200.**(0.)) t 'L = 200' w l lw 5. lt rgb 'dark-green' dt 1,\
#	"wireABCu0k200l200.dat1" using ($1):($5*2000.0*200.**(0.)) t '' w l lw 4. lt rgb 'dark-green' dt '.',\
#     "wireABCu0k200l400.dat1" using ($1):($5*2000.0*400.**(0.)) t 'L = 400' w l lw 4. lt rgb 'blue' dt 1,\
#	"wireABCu0k200l400.dat1" using ($1):($5*2000.0*400.**(0.)) t '' w l lw 4. lt rgb 'blue' dt '.',\
#      "wireABCu0k200l2000.dat1" using ($1):($5*2000.0*600.**(0.)) t 'L = 600' w l lw 5. lt 7 dt 1,\
#      "wireABCu0k200l2000.dat1" using ($1):($5*2000.0*600.**(0.)) t '' w l lw 4. lt 7 dt '.',\

#plot      "wireABCu0k8l2000.dat1" using ($1):($5*2000.0*800.**(0.)) t '{/Symbol Q}_{/ZapfDingbats=15 \110} = 8    ' w l lw 4. lt 1 dt 1,\
#      "wireABCu0k8l2000.dat1" using ($1):($5*2000.0*800.**(0.)) t '' w l lw 4. lt 1 dt '.',\

plot      "wireABCu0k20l2000.dat1" using ($1):($5*2000.0*1000.**(0.)) t '{/Symbol Q}_{/ZapfDingbats=15 \110} = 20  ' w l lw 5. lt rgb 'blue' dt '.',\
      "wireABCu0k20l2000.ddat2" using ($1):($5*2000.0*1000.**(0.)) t '' w l lw 4. lt 1 dt '.',\
      "wireABCu0k40l2000.dat1" using ($1):($5*2000.0*2000.**(0.)) t '{/Symbol Q}_{/ZapfDingbats=15 \110} = 40  ' w l lw 5. lt rgb 'cyan' dt 4,\
      "wireABCu0k40l2000.ddat2" using ($1):($5*2000.0*2000.**(0.)) t '' w l lw 4. lt rgb 'blue' dt '.',\
      "wireABCu0k80l2000.dat1" using ($1):($5*2000.0*4000.**(0.)) t '{/Symbol Q}_{/ZapfDingbats=15 \110} = 80  ' w l lw 5. lt 7 dt 5,\
      "wireABCu0k80l2000.ddat2" using ($1):($5*2000.0*4000.**(0.)) t '' w l lw 4. lt 7 dt '.',\
      "wireABCu0k100l2000.dat1" using ($1):($5*2000.0*600.**(0.)) t '{/Symbol Q}_{/ZapfDingbats=15 \110} = 100' w l lw 5. lt 1 dt 6,\
      "wireABCu0k100l2000.ddat2" using ($1):($5*2000.0*600.**(0.)) t '' w l lw 4. lt rgb 'dark-green' dt '.',\
      "wireABCu0k140l2000.dat1" using ($1):($5*2000.0*600.**(0.)) t '{/Symbol Q}_{/ZapfDingbats=15 \110} = 140' w l lw 2.5 lt 8 dt 1,\
      "wireABCu0k140l2000.ddat2" using ($1):($5*2000.0*600.**(0.)) t '' w l lw 4. lt 8 dt '.',\
#      "wireABCu0k200l20000.dat1" using ($1):(log(1-$5*2000.0)*10000.**(0.)) t 'L = 10000' w l lw 4. lt 8 dt 1,\
#      "wireABCu1k20l400.dat1" using ($1):($5*2000.0*400.**(0.)) t 'C@^o_{L/3}' w l lw 5. lt 1 dt 1,\
#      "wireABCu1k20l400.dat1" using ($1):($5*2000.0*400.**(0.)) t '' w l lw 4. lt 1 dt '.',\
#      "wireABCu0k200l20000.dat1" using ($1/500):(log(1-$5*2000.0)*500.**(0.)) t 'L = 12' w l lw 5. lt 7 dt 7,\
#      "data/kzUp1Sig-1L12p1A2.dat1" using ($1/12):(log(1-$5*2000.0)*4.**(0.)) t 'L = 12' w l lw 5. lt 8 dt 7,\
#      "data/tPyL8g1.dat1" using ($1):($5*2000.0) t 'L=4' pt 7 ps 0.8 lt 7,\
#      "data/CtL6g1.dat1" using ($1):($5*2000.0) t 'L=4' w l lw 5. lt rgb 'blue' dt 1 ,\
#do for [N=1:5] {
#plot func(N, x)
#pause -1
#}

#set size 0.525,0.445
#set origin 0.454,0.12
unset margin

set size 0.465,0.405
set origin 0.07,0.5

unset xlabel
unset ylabel
unset title

set ytics font "Arial,23"
set xtics  font "Arial,23"
#set xtics 0.02#, _, _  # offset _, _
unset key
#set xrange
unset label
unset xrange
unset yrange

set xlabel font "Arial,27"
#set ylabel font "Arial,27"
set label font "Arial,25"

#set format x "%.3f"

#set xrange [0.00:0.008]
set xlabel "1/L" offset 0, 0.8
#set ylabel "D({/Arial=9 L/2, t}) - D({/Arial=9 L/2, 0})" offset 3.5,0
#set xlabel "t L^{- 3}" offset 0,0.8
#set ylabel "L < C^+_{3L/8} C_{5L/8} + h.c. >" offset 2.5,0

set label 1 at screen 0.415, screen 0.825 "{/Symbol=28 q} =  0.5"

#f(x)=a*x + b 
#fit f(x) "dissipation-.dat1" u (1./$1):2 via a,b


#plot  f(x) notitle w l dashtype '_' lt 1, "dissipation-.dat1" using (1./$1):2 t '' pt 7 ps 2 # lt 8
#plot  "k-dissipationL30.dat1" using 1:2 notitle 

#set label 1 at screen 0.685, screen 0.58 "Asymptotic behavior" tc lt 4
#set arrow from graph 0, first -0.5 to graph 1, first -0.5 nohead lt 3 lw 2 dt ' _ _ '

#plot  "k-dissipationL25.dat1" using (1./$1):2 notitle w l lw 1.0 lt 8 dt '- - ',\
#      "k-dissipationL20.dat1" using (1./$1):2 notitle w l lw 1.0 lt 7 dt '- - ',\

#plot  "k-dissipationL20.dat1" using 1:2 t 'L = 16,  {/Arial=20 u} =  0.01' w l lw 1.0 lt 8 dt 3,\
#      "k-dissipationL25.dat1" using 1:2 t 'L = 32,  {/Arial=20 u} =  0.01' w l lw 1.0 lt 7 dt 3,\

unset multiplot
