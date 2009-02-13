CPP = g++
INC = 
CPPFLAGS = ${INC} -Wall -Wshadow -W -O3 -g
LIBS = 

PROGS = sample_data

all:	${PROGS}
clean:
	-rm -f ${PROGS} *.o

sample_data: sample_data.cc	
	$(CPP) $(CPPFLAGS) -o $@ sample_data.cc ${LIBS}
#measure_crosstalk: measure_crosstalk.cc	
#	$(CPP) $(CPPFLAGS) -o $@ measure_crosstalk.cc ${LIBS}

sample_data.o : sample_data.cc
#measure_crosstalk.o : measure_crosstalk.cc
