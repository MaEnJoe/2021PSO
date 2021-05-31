CXX = g++
CXXFLAG = -std=c++17 -fopenmp -O3 -lm
CLINFLAG = -fopenmp
OBJS = Particle.o \
	   main.o

PROGRAMNAME = a.out

all:${PROGRAMNAME}

run : ${PROGRAMNAME}
	./$<

${PROGRAMNAME}: ${OBJS}
	${CXX} ${CLINFLAG} ${OBJS} -o $@

${OBJS}: main.cpp Particle.cpp Particle.h utility.h
	${CXX} ${CXXFLAG} -c Particle.cpp main.cpp

clean:
	rm -f *.o *~ *.out PSO

fs : fs.cpp
	g++ fs.cpp -o fs
