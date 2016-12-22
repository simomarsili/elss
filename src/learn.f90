! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program learn
  use kinds
  use constants
  use fasta,       only : fasta_read
  use mpi_wrapper
  use dump,        only: read_rst,read_prm_unit,dump_energies,dump_rst
  use learn_command_line
  use units,       only: units_initialize,units_open
  use data,        only: data_read,data_average
  use random,      only: random_initialize,random_seq
  use map,         only: map_learn
  use scoring,     only: gauge,compute_scores,print_scores
  use mcmc,        only: mcmc_simulate,mcmc_compute_energy
  implicit none
  integer                 :: nvars        ! total number of variables
  integer                 :: nclasses     ! total number of classes
  integer                 :: nseqs     ! total number of samples
  real(kflt)              :: neff         ! effective number of seqs
  integer,    allocatable :: seq(:)    ! seq array
  integer,    allocatable :: seqs(:,:) ! data matrix
  integer,    allocatable :: seqs0(:,:)
  real(kflt), allocatable :: prm(:)       ! parameters array
  real(kflt), allocatable :: fmodel(:)    ! model frequencies
  real(kflt), allocatable :: fdata(:)     ! empirical frequencies
  real(kflt), allocatable :: scores(:,:)  ! final scores
  integer, allocatable    :: seqs_table(:,:)
  ! command line parameters
  integer                        :: udata        ! data unit
  character(len=string_size)     :: data_format  ! data format ('raw', 'table', 'protein')
  integer                        :: uwgt         ! ww unit
  real(kflt)                     :: wid          ! %id for weights calculation
  integer                        :: uprm         ! prm unit
  integer                        :: urst         ! rst unit
  integer                        :: useq         ! seq unit
  integer                        :: mc_nsweeps   ! number of MC sweeps per gradient estimate
  integer                        :: nupdate      ! stride for averages aupdate
  integer                        :: niter_agd    ! number of iter. (Nesterov alg.)
  integer                        :: niter_gd     ! number of iter. (gradient descent)
  real(kflt)                     :: lambda       ! Gaussian prior hyperparameter
  character(len=string_size)     :: mode         ! run mode: EVAL, SIM, LEARN
  integer                        :: rseed        ! random seed
  real(kflt)                     :: beta         ! temperature of the run
  character(len=long_string_size) :: prefix
  integer                         :: err
  character(long_string_size)     :: err_string
  integer                         :: dim1,dim2
  integer                         :: utrj,uene,ulog
  character(len=long_string_size) :: filename
  real(kflt)                      :: facc
  integer                         :: j
  real(kflt)                      :: etot,efields,ecouplings
  logical                         :: hot_start

  !================================================ set defaults

  nvars = 0
  nclasses = 0
  udata = 0
  data_format = ''
  uwgt = 0
  wid = -1
  uprm = 0
  urst = 0
  useq = 0
  mc_nsweeps = 0
  nupdate = 10
  niter_gd = 0
  niter_agd = 0
  lambda = 0.01_kflt
  mode = 'LEARN'
  rseed = 0
  beta = 1.0_kflt
  prefix = ''

  !================================================ init. unit identifiers

  call units_initialize()

  !================================================ read args

  call command_line_read(udata,data_format,uwgt,&
       wid,uprm,urst,useq,rseed,beta,mc_nsweeps,nupdate,niter_agd,&
       niter_gd,lambda,prefix,err)
  
  if (err /= 0) then
     stop
  end if

  !================================================ init. mpi

  call mpi_wrapper_initialize(err)

  !================================================ init. the random number generator

  call random_initialize(rseed,iproc)

  !================================================ read prms for the run

  if (uprm > 0) then
     call read_prm_unit(uprm,nvars,nclasses,&
          prm,data_format,err)
     if (err /= 0) then
        if(iproc == 0) write(0,*) 'ERROR ! cannot read from prm'
        call mpi_wrapper_finalize(err)
        stop
     end if
     close(uprm)
  end if

  !================================================ read restart file
  
  if (urst > 0) then
     call read_rst(urst,data_format,nvars,nclasses,iproc,nproc,seq,prm,err)
     if (err /= 0) then
        if(iproc == 0) write(0,*) 'ERROR ! cannot read from rst'
        call mpi_wrapper_finalize(err)
        stop
     end if
     close(urst)
  end if

  !================================================ read data

  call data_read(iproc,udata,data_format,uwgt,wid,&
       nvars,nclasses,nseqs,neff,seqs,err,err_string)

  if (err /= 0) then
     if(iproc == 0) write(0,*) 'ERROR ! cannot read from msa file'
     call mpi_wrapper_finalize(err)
     stop
  end if

  !================================================ allocate memory for the run and initialize

  dim1 = nvars * nclasses
  dim2 = nvars * (nvars - 1) * nclasses**2 / 2
  if (uprm == 0 .and. urst == 0) then
     allocate(prm(dim1+dim2),stat=err)
     prm = 0.0_kflt
  end if
  if (urst == 0) then
     ! allocate seq
     allocate(seq(nvars),stat=err)
     seq = 0
     call random_seq(nvars,nclasses,seq)
  end if
  ! allocate and fill up seqs_table
  allocate(seqs_table(nvars,nproc),stat=err)
  seqs_table = 0
  seqs_table(:,iproc+1) = seq
  CALL mpi_allgather(seq, nvars, MPI_INTEGER, seqs_table, nvars, MPI_INTEGER, MPI_COMM_WORLD, err)
  ! allocate model frequencies
  allocate(fmodel(dim1+dim2),stat=err)
  fmodel = 0.0_kflt
  
  close(udata)

  ! open output
  call units_open(trim(mode)//'.log','unknown',ulog,err)
  if(err /= 0) then
     if (iproc == 0) write(0,*) "error opening file ", trim(mode)//'.log'
     call mpi_wrapper_finalize(err)
     stop
  end if


  !================================================ print a header

  if ( iproc == 0 ) then
     write(ulog,'(a)')         '#'
     write(ulog,'(a)')         '#==========================================='
     write(ulog,'(a)')         '#_________________mcsg_v0.3.1_______________'
     write(ulog,'(a)')         '#'
     write(ulog,'(a)')         '#  mode  :    '//trim(mode)
     write(ulog,'(a)')         '#  format:    '//trim(data_format)
     if (urst > 0) write(ulog,'(a)')&
          '#  reading restart file           '
     if (uprm > 0) write(ulog,'(a)')&
          '#  reading parameter file         '
     write(ulog,'(a,1x,i8)')   '#  n. of procs                  = ', nproc
     write(ulog,'(a,1x,i8)')   '#  n. of variables              = ', nvars
     write(ulog,'(a,1x,i8)')   '#  n. of classes                = ', nclasses
     write(ulog,'(a,1x,i8)')   '#  n. of seqs                   = ', nseqs
     write(ulog,'(a,1x,i8)')   '#  stride (sweeps) for updates  = ', nupdate
     if (trim(mode) == 'SIM') &
          write(ulog,'(a,1x,i8)')   '#  n. of MC sweeps              = ', mc_nsweeps
     write(ulog,'(a,1x,f8.1)') '#  temperature factor           = ', 1.0_kflt / beta
     if (trim(mode) == 'LEARN') then
        if (wid > 0.0 .or. uwgt > 0.0) then
           if (wid > 0.0) then
              write(ulog,'(a,1x,f8.1)') &
                   '#  %id for reweighting          = ', wid
           else if (uwgt > 0.0) then
              write(ulog,'(a)') &
                   '#  reading weights from file...    '
           end if
           write(ulog,'(a,1x,f8.1)') '#  effective n. of seqs         = ', neff
        end if
        write(ulog,'(a,1x,f8.3)') '#  scaled prior (hyper-)parameter  = ', lambda
        write(ulog,'(a,1x,i8)')   '#  MC sweeps per gradient est.  = ', mc_nsweeps
        write(ulog,'(a,1x,i8)')   '#  n. of iterations (GD)        = ', niter_gd
        write(ulog,'(a,1x,i8)')   '#  n. of iterations (AGD)       = ', niter_agd
     end if
     write(ulog,'(a)')         '#'
     write(ulog,'(a)')         '#==========================================='
     write(ulog,'(a)')         '#'
  end if


  !================================================ allocate memory for map algorithm
  
  allocate(fdata(dim1+dim2),stat=err)
  allocate(scores(nvars,nvars),stat=err)
  fdata = 0.0_kflt
  scores = 0.0_kflt
  
  !================================================ compute averages from data
  
  call data_average(nvars,nclasses,nseqs,neff,seqs,fdata(1:dim1),fdata(dim1+1:dim1+dim2))
  
  !================================================ maximum a posteriori estimate of parameters
  
  ! inv. temperature for MAP estimation is set to 1
  call map_learn(nvars,nclasses,niter_agd,niter_gd,lambda,mc_nsweeps,&
       1.0_kflt,nupdate,data_format,ulog,fdata,seq,seqs_table,prm,fmodel)
  
  !================================================ compute and print scores
  
  !     call gauge(nvars,nclasses,prm(1:dim1),prm(dim1+1:dim1+dim2))
  
  !     call compute_scores(nvars,nclasses,prm(dim1+1:dim1+dim2),scores)
  
  !     if ( iproc == 0 ) call print_scores(nvars,scores)
  
  
  !================================================ finalize mpi
  
  call mpi_wrapper_finalize(err)

end program learn
