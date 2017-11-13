! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program sample
  use kinds
  use constants
  use fasta,       only : fasta_read
  use dump,        only: read_chk
  use sample_command_line
  use units,       only: units_initialize,units_open
  use random,      only: random_initialize,random_data
  use mcmc,        only: mcmc_simulate
  implicit none
  integer                 :: nvars        ! total number of variables
  integer                 :: nclasses     ! total number of classes
  integer,    allocatable :: seq(:)    ! seq array
  integer,    allocatable :: seqs0(:,:)
  real(kflt), allocatable :: prm(:)       ! parameters array
  real(kflt), allocatable :: fmodel(:)    ! model frequencies
  character(len=string_size) :: data_type    ! data tytpe ('unknown', 'protein', 'nuc_acid')
  ! command line parameters
  integer                    :: uchk         ! chk unit
  integer                    :: istart       ! index of starting point in chk file
  integer                    :: mc_nsweeps   ! number of MC sweeps per gradient estimate
  integer                    :: nupdate      ! stride for averages aupdate
  integer                    :: rseed        ! random seed
  real(kflt)                 :: beta         ! temperature of the run
  character(len=long_string_size) :: prefix
  integer                         :: iproc=0,nproc=1
  integer                         :: err
  integer                         :: dim1,dim2
  integer                         :: utrj,ulog
  character(len=long_string_size) :: filename
  real(kflt)                      :: facc
  logical                         :: hot_start

  !================================================ set defaults

  nvars = 0
  nclasses = 0
  data_type = 'unk'
  uchk = 0
  istart = 0
  mc_nsweeps = 1000
  nupdate = 10
  rseed = 0
  beta = 1.0_kflt
  prefix = ''

  !================================================ init. unit identifiers

  call units_initialize()

  !================================================ read args

  call command_line_read(uchk,rseed,beta,istart,mc_nsweeps,nupdate,prefix,err)

  if (err /= 0) stop

  !================================================ init. the random number generator

  call random_initialize(rseed,iproc)

  !================================================ read checkpoint file
  
  if (uchk > 0) then
     call read_chk(uchk,data_type,nvars,nclasses,istart,seq,prm,err)
     if (err /= 0) then
        write(0,*) 'ERROR ! cannot read from chk'
        stop
     end if
     close(uchk)
  end if

  !================================================ allocate memory for the run and initialize

  dim1 = nvars * nclasses
  dim2 = nvars * (nvars - 1) * nclasses**2 / 2
  if (uchk == 0) then
     ! allocate seq
     allocate(seq(nvars),stat=err)
     seq = 0
     call random_data(nclasses, seq)
  end if
  allocate(fmodel(dim1+dim2),stat=err)
  fmodel = 0.0_kflt

  if (iproc == 0) then
     
     ! open log file
     if (len_trim(prefix) == 0) then
        filename = 'log'
     else
        filename = trim(prefix)//'.log'        
     end if
     call units_open(trim(filename),'unknown',ulog,err)
     if(err /= 0) then
        write(0,*) "error opening file ", trim(filename)
        stop
     end if

     ! open trj file
     if (len_trim(prefix) == 0) then
        filename = 'trj'
     else
        filename = trim(prefix)//'.trj'        
     end if
     call units_open(filename,'unknown',utrj,err)
     if(err /= 0) then
        write(0,*) "error opening file ", trim(filename)
        stop
     end if

     ! print a header
     write(ulog, 101) adjustr(trim(data_type)), uchk, nvars,&
          nclasses, istart, mc_nsweeps, nupdate, beta
  end if

101 format(&
         '# elss-sample (elss v0.3.2)          '/& 
         '#                                   '/&
         '# data type:              ',    a12  /&
         '# chk file unit:          ',    i12  /&
         '# n. features:            ',    i12  /&
         '# n. classes:             ',    i12  /&
         '# index of starting point:',    i12  /&
         '# n. sweeps:              ',    i12  /&
         '# n. sweeps per update:   ',    i12  /&
         '# temperature factor:     ',  f12.3  /)

  !================================================ run a simulation
  
  hot_start = .false.
  call mcmc_simulate(nvars,nclasses,seq,&
       prm(1:dim1),prm(dim1+1:dim1+dim2),data_type,&
       fmodel(1:dim1),fmodel(dim1+1:dim1+dim2),&
       beta,mc_nsweeps,hot_start,nupdate,utrj,facc)

end program sample
