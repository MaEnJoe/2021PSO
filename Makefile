CXX = g++
CXXFLAG = -std=c++11
OBJS = Particle.o \
	   main.o

PROGRAMNAME = a.out

all:${PROGRAMNAME}

${PROGRAMNAME}: ${OBJS}
	${CXX} ${OBJS} -o $@

${OBJS}: main.cpp Particle.cpp Particle.h utility.h
	${CXX} ${CXXFLAG} -c Particle.cpp main.cpp

clean:
	rm -f *.o *~ *.out PSO



