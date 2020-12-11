set terminal pdfcairo color enhanced size 4in, 4in
set output "../fig/H_Js.pdf"
set xlabel "J_1"
set ylabel "J_2"
set zlabel "J_3"
set view equal xyz
set view 45,135,1,1
set ticslevel 0
set size square
splot "../data/H_Js.dat" u 1:2:3 lc 1 ps 0.1 notitle