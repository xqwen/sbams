main: main.o controller.o parser.o MVLR.o
	g++  -O3 main.o parser.o controller.o MVLR.o -lm -L /usr/local/lib -lgsl -lgslcblas -o sbams_mvlr
main.o: main.cc
	g++ -c main.cc
parser.o: parser.cc parser.h
	g++ -c parser.cc
controller.o: controller.cc controller.h
	g++ -c controller.cc
MVLR.o: MVLR.h MVLR.cc
	g++ -c MVLR.cc
clean:
	rm *.o sbams_mvlr
