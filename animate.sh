mkdir temp
convert $(ls -tr *.png) -delay 5 -morph 5 temp/%05d.jpg
avconv -r 25 -qscale 2 -i temp/%05d.jpg output.mp4
rm -r temp
