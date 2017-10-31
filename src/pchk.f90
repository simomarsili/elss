! Copyright (C) 2015-2017, Simone Marsili
! All rights reserved.
! License: BSD 3 clause

program pchk
  use kinds
  use constants
  use units, only: units_initialize, units_open, units_open_unf
  use arguments, only: read_opt, read_opt_arg
  use parser, only: parser_nfields,remove_comments, parse_next_line
  implicit none

  character(len=string_size) :: data_type
  integer                    :: nvars=0, nclasses=0
  integer,       allocatable :: seqs(:,:)
  real(kflt),    allocatable :: prm(:)
  real(kflt), allocatable :: array1(:), array2(:)
  integer                    :: id,ndata,err
  integer :: unt
  character(len=long_string_size) :: cmd
  character(len=long_string_size) :: arg
  character(len=long_string_size) :: source
  integer                         :: iarg, nargs
  integer :: nfields
  integer :: iv,jv,k,a,b
  integer :: n_digits = -1
  character(len=long_string_size) :: frmt
  character(len=long_string_size) :: line, parsed_line
  integer :: flag=0 ! 1: u -> f; 2: f-> u
  integer :: nerrs=0

  call units_initialize()

  ! read arguments
  call get_command(cmd)
  nargs = command_argument_count()

  ! call with no args
  if (nargs == 0) then
     write(0,100)
     stop
  end if

  source = ''
  iarg = 0
  args_loop: do while(iarg < nargs)
     iarg = iarg + 1
     call read_opt(iarg,nargs,arg,err)
     select case(trim(arg))
     case('-h','--help')
        ! print help and exit
        write(0,100)
        stop
     case('-u')
        ! unformatted file
        if (flag == 0) then
           flag = 1
        else
           write(0,*) "ERROR: either -u <chk_file > or -f <formatted_file> option"
           stop
        end if
        iarg = iarg + 1
        call read_opt_arg(iarg,nargs,source,err)
        if (err > 0) then
           write(0,*) 'ERROR ! option "-u": missing argument'
           stop
        end if
     case('-f')
        ! formatted file
        if (flag == 0) then
           flag = 2
        else
           write(0,*) "ERROR: either -u <chk_file > or -f <formatted_file> option"
           stop
        end if
        iarg = iarg + 1
        call read_opt_arg(iarg,nargs,source,err)
        if (err > 0) then
           write(0,*) 'ERROR ! option "-f": missing argument'
           stop
        end if
     case('-d', '--digits')
        ! the number of digits to the right of the decimal point
        ! default: print with full precision
        iarg = iarg + 1
        call read_opt_arg(iarg,nargs,n_digits,err)
        if (err == 2) then
           write(0,*) 'ERROR ! option "-d": missing argument'
           stop
        end if
     case default
        write(0,*) 'ERROR ! unkown argument '//trim(arg)
        write(0,100)
        stop
     end select
  end do args_loop

  select case(flag)
  case(1)
     ! open unformatted file
     call units_open_unf(source,'old',unt,err)
     ! read unformatted file
     read(unt) nvars
     read(unt) nclasses
     read(unt) data_type
     read(unt) ndata
     allocate(prm(nvars*nclasses + nvars*(nvars - 1)*nclasses**2/2),stat=err)
     if (ndata > 0) then
        allocate(seqs(nvars,ndata),stat=err)
        seqs = 0
        do id = 1,ndata
           read(unt) seqs(:,id)
        end do
     end if
     prm = 0.0_kflt
     read(unt) prm
     close(unt)
     ! dump formatted file to standard output
     source = trim(source)//'.txt'
     call units_open(source,'unknown',unt,err)
     write(unt,'(i10, 1x, a)') nvars, "# n. of variables"
     write(unt,'(i10, 1x, a)') nclasses, "# n. of classes"
     write(unt,'(a10, 1x, a)') trim(data_type), "# n. output data format"
     write(unt,'(i10, 1x, a)') ndata, "# n. samples in checkpoint"
     if (ndata > 0) then
        write(unt, '(a)') "# samples start here"
        do id = 1,ndata
           write(unt,'(1000(i2, 1x))') seqs(:,id)
        end do
     end if
     write(unt, '(a)') "# parms start here"
     k = 0
     write(frmt,*) '(i3,1x,1000f', n_digits + 4, '.', n_digits, ')'
     do iv = 1,nvars
        if (n_digits >= 0) then
           write(unt,frmt) iv,prm(k+1:k+nclasses)
        else
           write(unt,*) iv,prm(k+1:k+nclasses)
        end if
        k = k + nclasses
     end do
     write(frmt,*) '(2(i3,1x),1000f', n_digits + 4, '.', n_digits, ')'
     do jv = 1,nvars-1
        do iv = jv+1,nvars
           if (n_digits >= 0) then
              write(unt,frmt) jv, iv, prm(k+1:k+nclasses**2)
           else
              write(unt,*) jv, iv, prm(k+1:k+nclasses**2)
           end if
           k = k + nclasses**2
        end do
     end do
     close(unt)
  case (2)
     ! open formatted file
     call units_open(source,'old',unt,err)

     ! read nvars
     call parse_next_line(unt, line, parsed_line, nfields, err)
     if (err < 0) then
        write(0,*) "ERROR: invalid input file"
        stop
     end if
     read(parsed_line, *, iostat=err) nvars
     if (err < 0) then
        write(0,*) "check n. of features"
        stop
     end if

     ! read nclasses
     call parse_next_line(unt, line, parsed_line, nfields, err)
     if (err < 0) then
        write(0,*) "ERROR: invalid input file"
        stop
     end if
     read(parsed_line, *, iostat=err) nclasses
     if (err < 0) then
        write(0,*) "check n. of classes"
        stop
     end if

     ! read data type
     call parse_next_line(unt, line, parsed_line, nfields, err)
     if (err < 0) then
        write(0,*) "ERROR: invalid input file"
        stop
     end if
     read(parsed_line, *, iostat=err) data_type
     if (err < 0) then
        write(0,*) "check output data type format"
        stop
     end if

     ! read ndata
     call parse_next_line(unt, line, parsed_line, nfields, err)
     if (err < 0) then
        write(0,*) "ERROR: invalid input file"
        stop
     end if
     read(parsed_line, *, iostat=err) ndata
     if (err < 0) then
        write(0,*) "check n. of samples in .chk file"
        stop
     end if

     ! if the first lines are ok, allocate
     allocate(prm(nvars*nclasses + nvars*(nvars - 1)*nclasses**2/2),stat=err)
     prm = 0.0_kflt
     allocate(array1(nclasses), array2(nclasses**2), stat=err)

     ! read samples
     if (ndata > 0) then
        allocate(seqs(nvars,ndata),stat=err)
        id = 1
        do
           if (id > ndata) exit
           ! read sequence line
           read(unt,'(a)',iostat=err) line
           if (err < 0) exit
           call remove_comments(line)
           if (len_trim(line) == 0) cycle
           call parser_nfields(line, parsed_line, nfields)
           read(parsed_line,*) seqs(:,id)
           id = id + 1
        end do
     end if
     do
        call parse_next_line(unt, line, parsed_line, nfields, err)
        if (err < 0) then
           exit
        end if
        if (nfields == nclasses + 1) then
           read(parsed_line, *) iv
           k = (iv - 1) * nclasses
           read(parsed_line, *) iv, prm(k + 1 : k + nclasses)
        else if (nfields == nclasses**2 + 2) then
           read(parsed_line, *) jv, iv
           if (jv > iv) then
              write(0,*) "ERROR: when reading couplings line p q {couplings}, it must be q > p"
              stop
           end if
           k = (2 * nvars - jv) * (jv - 1) / 2 + iv - jv
           k = nvars * nclasses + (k - 1) * nclasses**2
           read(parsed_line, *) jv, iv, prm(k + 1 : k + nclasses**2)
        else
           write(0,*) "invalid checkpoint file"
           stop
        end if
     end do
     close(unt)
     ! dump an unformatted file <source>.chk
     source = trim(source)//'.chk'
     call units_open_unf(source,'unknown',unt,err)
     write(unt) nvars
     write(unt) nclasses
     write(unt) data_type
     write(unt) ndata
     do id = 1,ndata
        write(unt) seqs(:,id)
     end do
     write(unt) prm
     close(unt)

  end select

100 format(&
       'elss-pchk (elss v0.3.1)                                               '/&
       '                                                                    '/&
       'Usage:                                                              '/&
       '    elss-pchk [-d <digits>] (-u <chk_file> | -f <formatted_file>)   '/&
       '                                                                    '/&
       'Description:                                                        '/&
       '    Convert a checkpoint file to a human-readable plain text file   '/&
       '    and viceversa.                                                  '/&
       '                                                                    '/&
       'Required:                                                           '/&
       '-u <chk_file>                                                       '/&
       '    Print the content of <chk_file> to text file <chk_file>.txt     '/&
       '-f <formatted_file>                                                 '/&
       '    Read <formatted_file> and generate a checkpoint file named      '/&
       '    <formatted_file>.chk                                            '/&
       '                                                                    '/&
       'Options:                                                            '/&
       '-h, --help                                                          '/&
       '    Display this help and exit.                                     '/&
       '                                                                    '/&
       '-d, --digits <number_of_digits>, integer                            '/&
       '    For formatted floating points, the number of digits             '/&
       '    to the right of the decimal point.                              '/&
       '                                                                    '/&
       'Examples:                                                           '/&
       '                                                                    '/)

end program pchk
