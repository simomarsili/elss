! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program learn
  use kinds
  use constants
  use mpi_wrapper
  use dump,        only: read_chk,dump_energies
  use learn_command_line
  use units,       only: units_initialize,units_open
  use data,        only: data_read,data_average
  use random,      only: random_initialize,random_data
  use map,         only: map_learn
  use scoring,     only: gauge,compute_scores,print_scores
  use mcmc,        only: mcmc_simulate,mcmc_compute_energy
  implicit none
  integer                 :: nvars        ! total number of variables
  integer                 :: nclasses     ! total number of classes
  integer                 :: ndata     ! total number of samples
  real(kflt)              :: neff         ! effective number of data
  integer,    allocatable :: seq(:)    ! nvars seq array
  integer,    allocatable :: seqs(:,:)    ! (nvars, nreplicas) seqs array
  integer,    allocatable :: data(:,:) ! data matrix
  real(kflt), allocatable :: prm(:)       ! parameters array
  real(kflt), allocatable :: fmodel(:)    ! model frequencies
  real(kflt), allocatable :: fdata(:)     ! empirical frequencies
  real(kflt), allocatable :: scores(:,:)  ! final scores
  integer, allocatable    :: seqs_table(:, :, :) ! (nvars, nreplicas, nproc)
  character(len=string_size) :: data_type    ! data tytpe ('unk', 'bio', 'protein', 'nuc_acid')
  ! command line parameters
  integer                        :: udata        ! data unit
  integer                        :: uwgt         ! ww unit
  real(kflt)                     :: wid          ! %id for weights calculation
  integer                        :: uchk         ! chk unit
  integer                        :: nreplicas   ! number of replicas in the swarm
  integer                        :: mc_nsweeps   ! number of MC sweeps per gradient estimate
  integer                        :: nupdate      ! stride for averages aupdate
  integer                        :: niter        ! number of iter.
  real(kflt)                     :: lambda       ! Gaussian prior hyperparameter
  integer                        :: rseed        ! random seed
  real(kflt)                     :: beta         ! temperature of the run
  real(kflt)                     :: rate         ! learning rate
  character(len=string_size)     :: algorithm
  character(len=long_string_size) :: prefix
  integer                         :: err
  integer                         :: nv,nc,dim1,dim2
  integer                         :: ulog
  character(len=long_string_size) :: filename
  integer, allocatable            :: chk_data(:,:)

  ! set default values
  nvars = 0
  nclasses = 0
  udata = 0
  data_type = 'unk'
  uwgt = 0
  wid = -1
  neff =  0.0_kflt
  uchk = 0
  nreplicas = 1
  mc_nsweeps = 1000
  nupdate = 10
  algorithm = 'adam'
  rate = 0.01_kflt
  niter = 2000
  lambda = 0.01_kflt
  rseed = 0
  beta = 1.0_kflt
  prefix = ''

  ! init. unit identifiers
  call units_initialize()

  ! read command line args
  call command_line_read(udata, data_type, uwgt, wid, uchk, rseed, beta, &
       nreplicas, mc_nsweeps, nupdate, algorithm, rate, niter, lambda, &
       prefix, err)
  if (err /= 0) then
     call mpi_wrapper_finalize(err)
     stop
  end if
  
  ! init mpi
  call mpi_wrapper_initialize(err) 

  ! init. the random number generator
  call random_initialize(rseed,iproc) 

  ! read checkpoint file
  if (uchk > 0) then
     call read_chk(uchk, nvars, nclasses, data_type, chk_data, prm, err)
     if (err /= 0) then
        if(iproc == 0) &
             write(0,*) 'ERROR ! cannot read from chk file'
        call mpi_wrapper_finalize(err)
        stop
     end if
     allocate(seq(nvars), stat=err)
     call random_data(nclasses, seq)
     if (allocated(chk_data)) then
        if (nproc <= size(chk_data(1, :))) then
           seq = chk_data(:, iproc + 1)
        end if
     end if
     close(uchk)
  end if

  ! read data
  call data_read(iproc,udata,uwgt,wid,& 
       nv,nc,data_type,ndata,neff,data,err)
  close(udata)
  if (err /= 0) then
     call mpi_wrapper_finalize(err)
     stop
  end if
  if (nvars > 0) then
     if (nv /= nvars .or. nc > nclasses) then
        if(iproc == 0) write(0,*) 'ERROR ! data are not consistent with chk file'
        call mpi_wrapper_finalize(err)
        stop
     end if
  else
     nvars = nv
     nclasses = nc
  end if

  ! allocate memory for the run and initialize chains
  dim1 = nvars * nclasses
  dim2 = nvars * (nvars - 1) * nclasses**2 / 2
  if (uchk == 0) then
     allocate(prm(dim1 + dim2), stat=err)
     prm = 0.0_kflt
     allocate(seqs(nvars, nreplicas), seq(nvars), stat=err)
     call random_data(nclasses, seqs)
     seq = seqs(:, 1)
  end if
  
  
  ! allocate and fill up seqs_table
  allocate(seqs_table(nvars, nreplicas, nproc), stat=err)
  seqs_table = 0
  seqs_table(:, :, iproc+1) = seqs
  CALL mpi_allgather(seqs, nvars * nreplicas, MPI_INTEGER, seqs_table, nvars * nreplicas, MPI_INTEGER, MPI_COMM_WORLD, err)
  ! allocate model frequencies
  allocate(fmodel(dim1 + dim2), stat=err)
  fmodel = 0.0_kflt

  close(udata)

  if (iproc == 0) then
     
     ! open log file
     if (len_trim(prefix) == 0) then
        filename = 'log'
     else
        filename = trim(prefix)//'.log'        
     end if
     call units_open(filename,'unknown',ulog,err)
     if(err /= 0) then
        write(0,*) "error opening file ", trim(filename)
        call mpi_wrapper_finalize(err)
        stop
     end if

     ! and print a header
     write(ulog, 101) adjustr(trim(data_type)), uchk, nproc, nvars, &
          nclasses, ndata, nupdate, beta, uwgt, wid, neff, lambda, &
          nreplicas, mc_nsweeps, niter, adjustr(trim(algorithm))
  end if
  
101 format(&
         '# elss-learn (elss v0.4)              '/& 
         '#                                     '/&
         '# data type:                          ',    a12  /&
         '# chk file unit:                      ',    i12  /&
         '# n. procs:                           ',    i12  /&
         '# n. features:                        ',    i12  /&
         '# n. classes:                         ',    i12  /&
         '# n. samples:                         ',    i12  /&
         '# n. sweeps (per update):             ',    i12  /&
         '# temperature factor:                 ',  f12.3  /&
         '# weights file unit:                  ',    i12  /&
         '# %id threshold:                      ',  f12.3  /&
         '# neff:                               ',  f12.3  /&
         '# lambda:                             ',  f12.3  /&
         '# n. chains (per process):            ',    i12  /&
         '# n. sweeps (per iteration, process): ',    i12  /&
         '# n. iterations:                      ',    i12  /&
         '# optimization algorithm:             ',    a12  /&
         '# learning rate:                      ',  f12.3  /)

  !================================================ allocate memory for map algorithm
  
  allocate(fdata(dim1+dim2),stat=err)
  allocate(scores(nvars,nvars),stat=err)
  fdata = 0.0_kflt
  scores = 0.0_kflt
  
  !================================================ compute averages from data
  
  call data_average(nvars,nclasses,ndata,neff,data,fdata(1:dim1),fdata(dim1+1:dim1+dim2))
  
  !================================================ maximum a posteriori estimate of parameters
  
  ! inv. temperature for MAP estimation is set to 1
  call map_learn(algorithm,rate,nvars,nclasses,niter,lambda,mc_nsweeps,&
       1.0_kflt,nupdate,data_type,ulog,fdata,prefix,seq,seqs_table,prm,fmodel)

  !================================================ compute and print scores
  
  !     call gauge(nvars,nclasses,prm(1:dim1),prm(dim1+1:dim1+dim2))
  
  !     call compute_scores(nvars,nclasses,prm(dim1+1:dim1+dim2),scores)
  
  !     if ( iproc == 0 ) call print_scores(nvars,scores)
  
  
  !================================================ finalize mpi
  
  call mpi_wrapper_finalize(err)

end program learn
