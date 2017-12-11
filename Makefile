#
# Makefile
# weihao, 2017-12-11 16:37
#

ana: analysis.cpp mat.cpp
	g++ -g analysis.cpp mat.cpp -o ana.out

rec: recognize.cpp mat.cpp
	g++ -g recognize.cpp mat.cpp -o rec.out

# vim:ft=make
#
