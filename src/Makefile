EXE = mcdca mcdca-eval mcdca-sample mcdca-learn
FC = mpif90
FFLAGS = -O2 -fconvert=little-endian
# important: use -heap-arrays with the Intel ifort 
# ifort 
#FFLAGS = -fast -heap-arrays 
# gfortran 
#FFLAGS = -O3 -ffast-math
LDFLAGS = -s 
SRCS = main.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 errors.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90
OBJS = main.o units.o random.o mcmc.o kinds.o data.o command_line.o dump.o parser.o scoring.o constants.o errors.o averages.o cost.o mpi_wrapper.o map.o fasta.o
SRCS-learn = learn.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 learn_command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 errors.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90
OBJS-learn = learn.o units.o random.o mcmc.o kinds.o data.o learn_command_line.o dump.o parser.o scoring.o constants.o errors.o averages.o cost.o mpi_wrapper.o map.o fasta.o
SRCS-eval = eval.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 eval_command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 errors.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90
OBJS-eval = eval.o units.o random.o mcmc.o kinds.o data.o eval_command_line.o dump.o parser.o scoring.o constants.o errors.o averages.o cost.o mpi_wrapper.o map.o fasta.o
SRCS-sample = sample.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 sample_command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 errors.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90
OBJS-sample = sample.o units.o random.o mcmc.o kinds.o data.o sample_command_line.o dump.o parser.o scoring.o constants.o errors.o averages.o cost.o mpi_wrapper.o map.o fasta.o

all: $(EXE)
mcdca: $(OBJS) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
mcdca-learn: $(OBJS-learn) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-learn) $(LIBS)
mcdca-eval: $(OBJS-eval) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-eval) $(LIBS)
mcdca-sample: $(OBJS-sample) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-sample) $(LIBS)

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
cost.o: cost.f90 kinds.o
errors.o: errors.f90 kinds.o constants.o
command_line.o: command_line.f90 kinds.o units.o constants.o
learn_command_line.o: learn_command_line.f90 kinds.o units.o constants.o
eval_command_line.o: eval_command_line.f90 kinds.o units.o constants.o
sample_command_line.o: sample_command_line.f90 kinds.o units.o constants.o
scoring.o: scoring.f90 kinds.o constants.o 
data.o: data.f90 kinds.o parser.o constants.o fasta.o units.o 
averages.o: averages.f90 kinds.o
mcmc.o: mcmc.f90 kinds.o dump.o constants.o averages.o 
mpi_wrapper.o: mpi_wrapper.f90 kinds.o #mpi
map.o: map.f90 units.o kinds.o mpi_wrapper.o mcmc.o dump.o data.o constants.o cost.o # mpi
main.o: main.f90 command_line.o mpi_wrapper.o random.o kinds.o units.o data.o map.o dump.o scoring.o constants.o errors.o mcmc.o #mpi
learn.o: learn.f90 learn_command_line.o mpi_wrapper.o random.o kinds.o units.o data.o map.o dump.o scoring.o constants.o errors.o mcmc.o #mpi
eval.o: eval.f90 eval_command_line.o random.o kinds.o units.o data.o map.o dump.o scoring.o constants.o errors.o mcmc.o 
sample.o: sample.f90 sample_command_line.o random.o kinds.o units.o dump.o constants.o errors.o mcmc.o

.PHONY: default clean realclean debug showtargets

# clean: remove useless files, but keep executables
clean:
	$(RM) core TAGS ?*[~#] *.o __* ...* *.mod

# realclean: remove all regenerable files, including executables
realclean: clean
	$(RM) $(EXE) $(OBJS) $(OBJS-eval) $(OBJS-sample)

# debug: debug options (enable D comments lines, disable -O for ABSOFT)
debug:
	$(MAKE) "FFLAGS=-g " "OSTYPE=$(OSTYPE:-gnu=)"

# showtargets: list the most important targets of this makefile
showtargets:
	@ echo clean realclean debug fcheck fdepend 

distclean: realclean
