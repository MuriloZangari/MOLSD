# makefile moeadeda

GCC = g++

CC_FLAGS =  -c  -Wno-deprecated -fno-rtti $(DEBUGGING) -Wall $(OPTIMIZE) $(PROFILING) 
LN_FLAGS = $(DEBUGGING) $(OPTIMIZE) $(PROFILING) 

.SUFFIXES: .cpp


.cpp.o:
	$(GCC) $(CC_FLAGS) $<
 
SRC =  ArchivePareto.h ArchivePareto.cpp Distance.h Distance.cpp Global.h Global.cpp Group.h Group.cpp Lambda.h Lambda.cpp NEH.h NEH.cpp Problem.h Problem.cpp Subprob.h Subprob.cpp NTools.h NTools.cpp NMODEL.h NMODEL.cpp MOEAD.h MOEAD.cpp moead_main.cpp

OBJS = ArchivePareto.o Distance.o Global.o Group.o Lambda.o NEH.o Problem.o Subprob.o NTools.o NMODEL.o MOEAD.o moead_main.o



all: moeadeda

moeadeda  : $(OBJS)
	$(GCC) $(LN_FLAGS) $(OBJS) -o moeadeda

edit:
	emacs $(SRC) &

tags:
	etags *.h *.cpp

clean:
	rm -f $(OBJS) moeadeda










