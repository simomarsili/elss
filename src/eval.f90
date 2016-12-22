! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program eval
  use kinds
  use constants
  use fasta,       only : fasta_read
  use dump,        only: read_rst,read_prm_unit,dump_energies,dump_rst
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
  ! command line parameters
  integer                    :: udata        ! data unit
  character(len=string_size) :: data_format  ! data format ('raw', 'table', 'protein')
  integer                    :: uwgt         ! ww unit
  real(kflt)                 :: wid          ! %id for weights calculation
  integer                    :: uprm         ! prm unit
  integer                    :: urst         ! rst unit
  character(len=string_size) :: mode         ! run mode: EVAL, SIM, LEARN
  integer                    :: rseed        ! random seed
  character(len=long_string_size) :: prefix  ! prefix for output file
  integer                         :: iproc=0,nproc=1
  integer                         :: err
  character(long_string_size)     :: err_string
  integer                         :: dim1,dim2
  integer                         :: uene,ulog
  character(len=long_string_size) :: filename
  integer                         :: j
  real(kflt)                      :: etot,efields,ecouplings

  !================================================ set defaults

  nvars = 0
  nclasses = 0
  udata = 0
  data_format = ''
  uwgt = 0
  wid = -1
  uprm = 0
  urst = 0
  mode = 'EVAL'
  rseed = 0
  prefix = ''

  !================================================ init. unit identifiers

  call units_initialize()

  !================================================ read args

  call command_line_read(udata,data_format,uprm,urst,prefix,err)
  if (prefix=='') prefix = trim(mode)
  if (err /= 0) stop

  !================================================ init. the random number generator

  call random_initialize(rseed,iproc)

  !================================================ read prms for the run

  if (uprm > 0) then
     call read_prm_unit(uprm,nvars,nclasses,&
          prm,data_format,err)
     if (err /= 0) then
        stop
     end if
     close(uprm)
  end if

  !================================================ read restart file

  if (urst > 0) then
     call read_rst(urst,data_format,nvars,nclasses,iproc,nproc,seq,prm,err)
     if (err /= 0) then
        stop
     end if
     close(urst)
  end if

  !================================================ read data

  if (udata > 0) then ! TODO: udata should always be > 0 for eval; check in command line and remove if then

     call data_read(iproc,udata,data_format,uwgt,wid,&
          nvars,nclasses,nseqs,neff,seqs,err,err_string)

     if (err /= 0) then
        stop
     end if

     !================================================ allocate memory for the run and initialize

     dim1 = nvars * nclasses
     dim2 = nvars * (nvars - 1) * nclasses**2 / 2
     if (uprm == 0 .and. urst == 0) then

        allocate(seq(nvars),stat=err)
        allocate(prm(dim1+dim2),stat=err)
        seq = 0
        prm = 0.0_kflt

        !================================================ initialize system configuration

        call random_seq(nvars,nclasses,seq)

     end if
     close(udata)
  end if

  ! open log file
  filename = trim(prefix)//'.log'
  call units_open(trim(filename),'unknown',ulog,err)
  if(err /= 0) then
     write(0,*) "error opening file ", trim(filename)//'.log'
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
     write(ulog,'(a,1x,i8)')   '#  n. of variables              = ', nvars
     write(ulog,'(a,1x,i8)')   '#  n. of classes                = ', nclasses
     write(ulog,'(a,1x,i8)')   '#  n. of seqs                   = ', nseqs
     write(ulog,'(a)')         '#'
     write(ulog,'(a)')         '#==========================================='
     write(ulog,'(a)')         '#'
  end if

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
