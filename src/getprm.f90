! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program getprm
  use kinds
  use constants
  use units, only: units_open_unf
  use arguments, only: read_opt,read_opt_arg
  
  integer                    :: unt=11
  character(len=string_size) :: data_type
  character(len=string_size) :: data_format
  integer                    :: nvars,nclasses
  integer,       allocatable :: seq(:,:)
  real(kflt),    allocatable :: prm(:)
  integer                    :: p,np,err
  integer                    :: nv,nc
  character(len=long_string_size) :: cmd
  integer                         :: nargs
  character(len=long_string_size) :: arg
  integer                         :: iarg
  logical                         :: file_exists
  character(len=long_string_size) :: syntax = 'syntax'
  character(len=long_string_size) :: rst_file
  integer :: iv,jv,is,js,k
  character(len=1) :: protein_alphabet(21) = &
       ['A','C','D','E','F','G','H','I','K',&
       'L','M','N','P','Q','R','S','T','V','W','Y','-']
  character(len=1) :: nuc_acid_alphabet(6) = &
       ['A','C','G','T','U','-']
  
  ! read arguments
  call get_command(cmd)
  nargs = command_argument_count()
  
  ! call with no args
  if (nargs == 0) stop

  iarg = 0
  args_loop: do while(iarg < nargs)
     call read_opt(iarg,nargs,arg,err)
     select case(trim(arg))
     case('-h','--help')
        ! print help and exit
        write(0,*) trim(syntax)
        stop
     case('-i')
        ! prm file
        call read_opt_arg(iarg,nargs,rst_file,err)
        if (err == 1) then
           write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
           stop
        end if
     case default
        write(0,*) 'ERROR ! unkown argument '//trim(arg)
        stop
     end select
  end do args_loop

  call units_open_unf(rst_file,'old',unt,err)
  
  ! read binary file
  read(unt) data_type
  read(unt) data_format
  read(unt) nvars
  read(unt) nclasses
  read(unt) np
  allocate(prm(nvars*nclasses + nvars*(nvars - 1)*nclasses**2/2),stat=err)
  if (np > 0) then
     ! seqs found in restart
     allocate(seq(nvars,np),stat=err)
     seq = 0
     do p = 1,np
        read(unt) seq(:,p)
     end do
  end if
  prm = 0.0_kflt
  read(unt) prm

  ! dump file
  write(*,'(a)') '# <data_type> <data_format> <nvars> <nclasses> <nseq>'
  write(*,'(a,1x,a,3(1x,i4))') trim(data_type), trim(data_format), nvars, nclasses, np
  if (np > 0) then
     do p = 1,np
        write(*,'(a,1000(1x,i2))') '# ',seq(:,p)
     end do
  end if
  k = 0
  do iv = 1,nvars
     do is = 1,nclasses
        k = k + 1
        write(*,'(i3,1x,i2,1x,f8.4)') iv,is,prm(k)
     end do
  end do
  do jv = 1,nvars-1
     do iv = jv+1,nvars
        do js = 1,nclasses
           do is = 1,nclasses
              k = k + 1
              ! order is inverted for printing 
              write(*,'(2(i3,1x),2(i2,1x),f8.4)') jv,iv,js,is,prm(k)
           end do
        end do
     end do
  end do

end program getprm
