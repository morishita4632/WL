# set terminal pdfcairo color enhanced size 4in, 4in
# set output "../fig/T_Js.pdf"
set xlabel "J_1"
set ylabel "J_2"
set size square
plot "../data/T_Js.dat" u 1:2 lc 1 ps 0.05 notitle,\
     "../data/T_Js.dat" u 1:3 lc 1 ps 0.05 notitle,\
     "../data/T_Js.dat" u 2:3 lc 1 ps 0.05 notitle,\
     "../data/T_Js.dat" u 2:1 lc 1 ps 0.05 notitle,\
     "../data/T_Js.dat" u 3:1 lc 1 ps 0.05 notitle,\
     "../data/T_Js.dat" u 3:2 lc 1 ps 0.05 notitle