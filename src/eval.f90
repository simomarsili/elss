! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program eval
  use kinds
  use constants
  use dump,        only: read_chk,read_prm_unit,dump_energies
  use eval_command_line
  use units,       only: units_initialize,units_open
  use data,        only: data_read,data_average
  use random,      only: random_initialize,random_seq
  use mcmc,        only: mcmc_compute_energy
  implicit none
  integer                 :: nvars        ! total number of variables
  integer                 :: nclasses     ! total number of classes
  integer                 :: nseqs     ! total number of samples
  real(kflt)              :: neff         ! effective number of seqs
  integer,    allocatable :: seq(:)    ! seq array
  integer,    allocatable :: seqs(:,:) ! data matrix
  real(kflt), allocatable :: prm(:)       ! parameters array
  character(len=string_size) :: data_type    ! data tytpe ('unk', 'bio', 'protein', 'nuc_acid')
  ! command line parameters
  integer                    :: udata        ! data unit
  integer                    :: uwgt         ! ww unit
  real(kflt)                 :: wid          ! %id for weights calculation
  integer                    :: uprm         ! prm unit
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
  uprm = 0
  uchk = 0
  rseed = 0
  prefix = 'eval'

  !================================================ init. unit identifiers

  call units_initialize()

  !================================================ read args

  call command_line_read(udata,data_type,uprm,uchk,prefix,err)
  if (err /= 0) stop

  !================================================ init. the random number generator

  call random_initialize(rseed,iproc)

  !================================================ read prms for the run

  if (uprm > 0) then
     call read_prm_unit(uprm,nvars,nclasses,&
          prm,data_type,err)
     if (err /= 0) then
        stop
     end if
     close(uprm)
  end if

  !================================================ read checkpoint file

  if (uchk > 0) then
     call read_chk(uchk,data_type,nvars,nclasses,iproc,nproc,seq,prm,err)
     if (err /= 0) then
        stop
     end if
     close(uchk)
  end if

  !================================================ read data

  call data_read(iproc,udata,uwgt,wid,&
       nvars,nclasses,data_type,nseqs,neff,seqs,err)
  
  if (err /= 0) then
     stop
  end if
  
  !================================================ allocate memory for the run and initialize

  dim1 = nvars * nclasses
  dim2 = nvars * (nvars - 1) * nclasses**2 / 2
  if (uprm == 0 .and. uchk == 0) then
     
     allocate(seq(nvars),stat=err)
     allocate(prm(dim1+dim2),stat=err)
     seq = 0
     prm = 0.0_kflt
     
     !================================================ initialize system configuration
     
     call random_seq(nvars,nclasses,seq)
     
  end if
  close(udata)

  if (iproc == 0) then
     
     ! open log file
     filename = trim(prefix)//'.log'
     call units_open(trim(filename),'unknown',ulog,err)
     if(err /= 0) then
        write(0,*) "error opening file ", trim(filename)//'.log'
        stop
     end if

     ! print a header
     write(ulog, 101) adjustr(trim(data_type)), uchk, uprm, nvars,&
          nclasses, nseqs
  end if

101 format(&
         '# elss-eval (elss v0.2.1)           '/& 
         '#                                   '/&
         '# data type:              ',    a12  /&
         '# chk file unit:          ',    i12  /&
         '# prm file unit:          ',    i12  /&
         '# n. features:            ',    i12  /&
         '# n. classes:             ',    i12  /&
         '# n. samples:             ',    i12  /)
         
  write(filename,*) iproc
  call units_open(trim(adjustl(filename))//'.ene','unknown',uene,err)
  if(err /= 0) then
     if (iproc == 0) write(0,*) "error opening file ", trim(adjustl(filename))//'.ene'
     stop
  end if
  
  do j = 1,nseqs
     call mcmc_compute_energy(nvars,nclasses,seqs(:,j),prm(1:dim1),prm(dim1+1:dim1+dim2),efields,ecouplings,etot)
     call dump_energies(uene,etot,efields,ecouplings)
  end do

end program eval
