#INSTALLDIR=/usr/local/bin
INSTALLDIR=~/.local/bin

EXE = elss-eval elss-sample elss-learn
TOOLS = pchk
FC = mpif90
FFLAGS = -O2 -fconvert=little-endian -std=f2008
#FFLAGS = -fast -heap-arrays # ifort (important: use -heap-arrays)
#FFLAGS = -O3 -ffast-math # gfortran
LDFLAGS = -s 
SRCS-learn = learn.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 learn_command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90 arguments.f90
OBJS-learn = learn.o units.o random.o mcmc.o kinds.o data.o learn_command_line.o dump.o parser.o scoring.o constants.o averages.o cost.o mpi_wrapper.o map.o fasta.o arguments.o
SRCS-eval = eval.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 eval_command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90 arguments.f90
OBJS-eval = eval.o units.o random.o mcmc.o kinds.o data.o eval_command_line.o dump.o parser.o scoring.o constants.o averages.o cost.o mpi_wrapper.o map.o fasta.o arguments.o
SRCS-sample = sample.f90 units.f90 random.f90 mcmc.f90 kinds.f90 data.f90 sample_command_line.f90 dump.f90 parser.f90 scoring.f90 constants.f90 averages.f90 cost.f90 mpi_wrapper.f90 map.f90 fasta.f90 arguments.f90
OBJS-sample = sample.o units.o random.o mcmc.o kinds.o data.o sample_command_line.o dump.o parser.o scoring.o constants.o averages.o cost.o mpi_wrapper.o map.o fasta.o arguments.o
SRCS-pchk = pchk.f90 kinds.f90 constants.f90 units.f90 arguments.f90
OBJS-pchk = getprm.o kinds.o constants.o units.o arguments.o

all: $(EXE) $(TOOLS)
elss-learn: $(OBJS-learn) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-learn) $(LIBS)
elss-eval: $(OBJS-eval) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-eval) $(LIBS)
elss-sample: $(OBJS-sample) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-sample) $(LIBS)
pchk: $(OBJS-pchk) 
	$(FC) $(LDFLAGS) -o $@ $(OBJS-pchk) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90 .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

%.o: %.mod

kinds.o: kinds.f90
mpi_wrapper.o: mpi_wrapper.f90 kinds.o #mpi
constants.o: constants.f90
arguments.o: arguments.f90 constants.o
fasta.o: fasta.f90 constants.o
parser.o: parser.f90
units.o: units.f90 kinds.o
random.o: random.f90 kinds.o
dump.o: dump.f90 kinds.o units.o random.o parser.o constants.o fasta.o
cost.o: cost.f90 kinds.o
learn_command_line.o: learn_command_line.f90 kinds.o units.o constants.o arguments.o
eval_command_line.o: eval_command_line.f90 kinds.o units.o constants.o arguments.o
sample_command_line.o: sample_command_line.f90 kinds.o units.o constants.o arguments.o
scoring.o: scoring.f90 kinds.o constants.o 
data.o: data.f90 kinds.o parser.o constants.o fasta.o units.o 
averages.o: averages.f90 kinds.o
mcmc.o: mcmc.f90 kinds.o dump.o constants.o averages.o 
map.o: map.f90 units.o kinds.o mpi_wrapper.o mcmc.o dump.o data.o constants.o cost.o # mpi
learn.o: learn.f90 learn_command_line.o mpi_wrapper.o random.o kinds.o units.o data.o map.o dump.o scoring.o constants.o mcmc.o #mpi
eval.o: eval.f90 eval_command_line.o random.o kinds.o units.o data.o map.o dump.o scoring.o constants.o mcmc.o 
sample.o: sample.f90 sample_command_line.o random.o kinds.o units.o dump.o constants.o mcmc.o
pchk.o = pchk.f90 kinds.o constants.o units.o arguments.o

.PHONY: default install uninstall clean realclean debug showtargets

install: $(EXE)
	mkdir -p $(INSTALLDIR)
	cp elss-* $(INSTALLDIR)

uninstall:
	rm -f $(INSTALLDIR)/elss-*

# clean: remove useless files, but keep executables
clean:
	$(RM) core TAGS ?*[~#] *.o __* ...* *.mod

# realclean: remove all regenerable files, including executables
realclean: clean
	$(RM) $(EXE) $(OBJS-learn) $(OBJS-eval) $(OBJS-sample)

# debug: debug options (enable D comments lines, disable -O for ABSOFT)
debug:
	$(MAKE) "FFLAGS=-g " "OSTYPE=$(OSTYPE:-gnu=)"

# showtargets: list the most important targets of this makefile
showtargets:
	@ echo clean realclean debug fcheck fdepend 

distclean: realclean
