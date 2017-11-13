! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program eval
  use kinds
  use constants
  use dump,        only: read_chk,dump_energies
  use eval_command_line
  use units,       only: units_initialize,units_open
  use data,        only: data_read,data_average
  use random,      only: random_initialize
  use mcmc,        only: mcmc_compute_energy
  implicit none
  integer                 :: nvars        ! total number of variables
  integer                 :: nclasses     ! total number of classes
  integer                 :: ndata     ! total number of samples
  real(kflt)              :: neff         ! effective number of data
  integer,    allocatable :: seq(:)    ! seq array
  integer,    allocatable :: data(:,:) ! data matrix
  real(kflt), allocatable :: prm(:)       ! parameters array
  character(len=string_size) :: data_type    ! data tytpe ('unk', 'bio', 'protein', 'nuc_acid')
  ! command line parameters
  integer                    :: udata        ! data unit
  integer                    :: uwgt         ! ww unit
  real(kflt)                 :: wid          ! %id for weights calculation
  integer                    :: uchk         ! chk unit
  integer                    :: rseed        ! random seed
  character(len=long_string_size) :: prefix  ! prefix for output file
  integer                         :: iproc=0,nproc=1
  integer                         :: err
  integer                         :: dim1,dim2
  integer                         :: uene,ulog
  character(len=long_string_size) :: filename
  integer                         :: j
  real(kflt)                      :: etot,efields,ecouplings

  !================================================ set defaults

  nvars = 0
  nclasses = 0
  udata = 0
  data_type = 'unk'
  uwgt = 0
  wid = -1
  uchk = 0
  rseed = 0
  prefix = ''

  !================================================ init. unit identifiers

  call units_initialize()

  !================================================ read args

  call command_line_read(udata,data_type,uchk,prefix,err)
  if (err /= 0) stop

  !================================================ init. the random number generator

  call random_initialize(rseed,iproc)

  !================================================ read checkpoint file

  if (uchk > 0) then
     call read_chk(uchk,data_type,nvars,nclasses,iproc,seq,prm,err)
     if (err /= 0) then
        stop
     end if
     close(uchk)
  end if

  !================================================ read data

  call data_read(iproc,udata,uwgt,wid,&
       nvars,nclasses,data_type,ndata,neff,data,err)
  
  if (err /= 0) then
     stop
  end if
  
  !================================================ allocate memory for the run and initialize

  dim1 = nvars * nclasses
  dim2 = nvars * (nvars - 1) * nclasses**2 / 2

  close(udata)

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
        filename = 'ene'
     else
        filename = trim(prefix)//'.ene'        
     end if
     call units_open(filename,'unknown',uene,err)
     if(err /= 0) then
        write(0,*) "error opening file ", trim(filename)
        stop
     end if
     
     ! print a header
     write(ulog, 101) adjustr(trim(data_type)), uchk, nvars,&
          nclasses, ndata
  end if

101 format(&
         '# elss-eval (elss v0.3.2)           '/& 
         '#                                   '/&
         '# data type:              ',    a12  /&
         '# chk file unit:          ',    i12  /&
         '# n. features:            ',    i12  /&
         '# n. classes:             ',    i12  /&
         '# n. samples:             ',    i12  /)
         
  do j = 1,ndata
     call mcmc_compute_energy(nvars,nclasses,data(:,j),prm(1:dim1),prm(dim1+1:dim1+dim2),efields,ecouplings,etot)
     call dump_energies(uene,etot,efields,ecouplings)
  end do

end program eval
