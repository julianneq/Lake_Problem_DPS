# Makefile for traditional lake problem
CXX = mpiCC
CXXFLAGS = -c -Wall
INC = -I /opt/boost/1.54.0/include
UNAME_S = $(shell uname -s)

ifneq (, $(findstring SunOS, $(UNAME_S)))
    LIBS += -lnsl -lsocket -lresolv
endif

all: LakeITparallel

LakeITparallel: lakeIT.o moeaframework.o borgms.o mt19937ar.o
	$(CXX) lakeIT.o moeaframework.o borgms.o mt19937ar.o -o LakeITparallel

lakeIT.o: main-lake-IT.cpp ./../borg/moeaframework.h ./../boostutil.h ./../borg/borgms.h
	$(CXX) $(CXXFLAGS) $(INC) main-lake-IT.cpp -o lakeIT.o

moeaframework.o: ./../borg/moeaframework.c ./../borg/moeaframework.h
	$(CXX) $(CXXFLAGS) ./../borg/moeaframework.c

borgms.o: ./../borg/borgms.c ./../borg/borgms.h \
	./../borg/borg.h
	$(CXX) $(CXXFLAGS) ./../borg/borgms.c $(LIBS)

mt19937ar.o: ./../borg/mt19937ar.c ./../borg/mt19937ar.h
	$(CXX) $(CXXFLAGS) ./../borg/mt19937ar.c $(LIBS)

clean:
	rm -rf *.o
	rm LakeITparallel
