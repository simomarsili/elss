ELSS = elss
FC = mpif90
FFLAGS = -O2 -fconvert=little-endian
# important: use -heap-arrays with the Intel ifort 
# ifort 
#FFLAGS = -fast -heap-arrays 
# gfortran 
#FFLAGS = -O3 -ffast-math
LDFLAGS = -s 
SRCS = main.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 errors.f90 averages.f90 likelihood.f90 mpi_wrapper.f90 map.f90 fasta.f90
OBJS = main.o units.o random.o mcmc.o kinds.o data.o command_line.o dump.o parser.o scoring.o constants.o errors.o averages.o likelihood.o mpi_wrapper.o map.o fasta.o

$(ELSS): $(OBJS) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90 .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

#mpi_wrapper.o:
#	$(MPI_FC) $(FFLAGS) -c mpi_wrapper.f90
#map.o:
#	$(MPI_FC) $(FFLAGS) -c map.f90
#main.o:
#	$(MPI_FC) $(FFLAGS) -c main.f90

kinds.o: kinds.f90
constants.o: constants.f90
fasta.o: fasta.f90 constants.o
parser.o: parser.f90
units.o: units.f90 kinds.o
random.o: random.f90 kinds.o
dump.o: dump.f90 kinds.o units.o random.o parser.o constants.o fasta.o
likelihood.o: likelihood.f90 kinds.o
errors.o: errors.f90 kinds.o constants.o
command_line.o: command_line.f90 kinds.o units.o constants.o
scoring.o: scoring.f90 kinds.o constants.o 
data.o: data.f90 kinds.o parser.o constants.o fasta.o units.o 
averages.o: averages.f90 kinds.o
mcmc.o: mcmc.f90 kinds.o dump.o constants.o averages.o 
mpi_wrapper.o: mpi_wrapper.f90 kinds.o #mpi
map.o: map.f90 units.o kinds.o mpi_wrapper.o mcmc.o dump.o data.o constants.o likelihood.o # mpi
main.o: main.f90 command_line.o mpi_wrapper.o random.o kinds.o units.o data.o map.o dump.o scoring.o constants.o errors.o mcmc.o #mpi

.PHONY: default clean realclean debug showtargets

# clean: remove useless files, but keep executables
clean:
	$(RM) core TAGS ?*[~#] *.o __* ...* *.mod

# realclean: remove all regenerable files, including executables
realclean: clean
	$(RM) $(ELSS) $(OBJS)

# debug: debug options (enable D comments lines, disable -O for ABSOFT)
debug:
	$(MAKE) "FFLAGS=-g " "OSTYPE=$(OSTYPE:-gnu=)"

# showtargets: list the most important targets of this makefile
showtargets:
	@ echo clean realclean debug fcheck fdepend 

distclean: realclean
