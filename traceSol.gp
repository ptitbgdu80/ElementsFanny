# Output dans un fichier pdf
set terminal pdfcairo enhanced color font "Arial, 12" linewidth 2 fontscale 1.0 \
    size 10,6 # inches

set palette rgbformulae 7,5,15

set view 60,45
#set hidden3d
#set zrange[0:1]
set ztics
unset key

set output "resultat_p.pdf"
plot "resultat_p.txt" u 1:2:3 with image

set palette rgbformulae 22,13,-31

set output "resultat_u.pdf"
plot[-0.01:1.1][-0.01:1.1] "resultat_u.txt" u 1:2:3:4:5 with vectors head size 0.1,5,30 filled lc palette
