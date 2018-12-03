# Output dans un fichier pdf
set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 10,6 # inches

set palette rgbformulae 7,5,15
set isosample 1000,1000
set view 20,45
set zrange[0:1]
set ztics 0,0.5,1

#set output "p5.pdf"
#splot [0:1] [0:1] 16*x*y*(x-1)*(y-1) with pm3d

#set output "p7.pdf"
#splot [0:1] [0:1] 4*(x-0.5)*(y-0.5)*(x-1)*y with pm3d

#set output "p4.pdf"
#splot [0:1] [0:1] -8*(x-0.5)*y*(x-1)*(y-1) with pm3d

#set output "p1.pdf"
#splot [0:1] [0:1] 4*(x-0.5)*(y-0.5)*(x-1)*(y-1) with pm3d

set output "base_de_Q1.pdf"
splot [0:2][0:2] x<1 && y<1 ? x*y :  x<1 && y>1 ? -x*(y-2) : x>1 && y < 1 ? -(x-2)*y :  x>1 && y>1 ? (x-2)*(y-2) : 0 with pm3d
