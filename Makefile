PROG =	VM.exe

SRCS =	zone.f90 quietstart.f90 particules.f90 \
	maxwell.f90 main.f90

OBJS =	zone.o quietstart.o particules.o \
	maxwell.o main.o

LIBS =	-llapack

F90      = gfortran
F90FLAGS = -O3 -Wall
LDFLAGS  = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod fort.* *.gnu

debug:	
	make F90FLAGS=-g FFLAGS=-g

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90
.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90
.f90.ps:
	a2ps $*.f90 -o $*.ps

main.o : main.f90 particules.o maxwell.o zone.o 
particules.o : particules.f90 quietstart.o zone.o
quietstart.o : quietstart.f90 zone.o
maxwell.o : maxwell.f90 zone.o
zone.o : zone.f90
