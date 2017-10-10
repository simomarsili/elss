! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program pchk
  use kinds
  use constants
  use units, only: units_open_unf
  use arguments, only: read_opt,read_opt_arg
  
  integer                    :: unt=11
  character(len=string_size) :: data_type
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
  character(len=long_string_size) :: chk_file
  integer :: iv,jv,is,js,k
  integer :: n_digits = -1
  character(len=1) :: protein_alphabet(21) = &
       ['A','C','D','E','F','G','H','I','K',&
       'L','M','N','P','Q','R','S','T','V','W','Y','-']
  character(len=1) :: nuc_acid_alphabet(6) = &
       ['A','C','G','T','U','-']
  character(len=long_string_size) :: frmt
  
  ! read arguments
  call get_command(cmd)
  nargs = command_argument_count()
  
  ! call with no args
  if (nargs == 0) stop

  iarg = 0
  args_loop: do while(iarg < nargs)
     iarg = iarg + 1
     call read_opt(iarg,nargs,arg,err)
     select case(trim(arg))
     case('-h','--help')
        ! print help and exit
        write(0,*) trim(syntax)
        stop
     case('-i')
        ! prm file
        iarg = iarg + 1
        call read_opt_arg(iarg,nargs,chk_file,err)
        if (err == 2) then
           write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
           stop
        end if
     case('-d', '--digits')
        ! the number of digits to the right of the decimal point
        ! default: print with full precision
        iarg = iarg + 1
        call read_opt_arg(iarg,nargs,n_digits,err)
        if (err == 2) then
           write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <n>'
           stop
        end if
     case default
        write(0,*) 'ERROR ! unkown argument '//trim(arg)
        stop
     end select
  end do args_loop

  write(frmt,*) 

  call units_open_unf(chk_file,'old',unt,err)
  
  ! read binary file
  read(unt) data_type
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
  write(*,'(a)') '# <data_type> <nvars> <nclasses> <nseq>'
  write(*,'(a,3(1x,i4))') trim(data_type), nvars, nclasses, np
  if (np > 0) then
     do p = 1,np
        write(*,'(a,1000(1x,i2))') '# ',seq(:,p)
     end do
  end if
  k = 0
  write(frmt,*) '(i3,1x,1000f', n_digits + 4, '.', n_digits, ')'
  do iv = 1,nvars
     !write(*,'(i3,1x,1000f8.4)') iv,prm(k+1:k+nclasses)
     !write(*,*) iv,prm(k+1:k+nclasses)
     if (n_digits >= 0) then
        write(*,frmt) iv,prm(k+1:k+nclasses)
     else
        write(*,*) iv,prm(k+1:k+nclasses)
     end if
     k = k + nclasses
  end do
  write(frmt,*) '(2(i3,1x),1000f', n_digits + 4, '.', n_digits, ')'
  do iv = 1,nvars-1
     do jv = iv+1,nvars
        ! order is inverted for printing 
        !write(*,'(2(i3,1x),1000f8.4)') jv,iv,prm(k+1:k+nclasses**2)
        ! (iv, jv) pair of features; inner loop over the states of jv
        if (n_digits >= 0) then
           write(*,frmt) iv,jv,prm(k+1:k+nclasses**2)
        else
           write(*,*) iv,jv,prm(k+1:k+nclasses**2)
        end if
        k = k + nclasses**2
     end do
  end do

end program pchk
