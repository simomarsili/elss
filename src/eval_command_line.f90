! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module eval_command_line
  use kinds
  use constants
  implicit none
  private
  public :: command_line_read

contains
  
  subroutine command_line_read(udata,data_type,uchk,prefix,error_code)
    use units, only: units_open,units_open_unf
    use arguments, only: read_opt,read_opt_arg
    integer,                    intent(inout) :: udata
    character(len=*),           intent(inout) :: data_type
    integer,                    intent(inout) :: uchk
    character(len=*),           intent(out) :: prefix
    integer,                    intent(out)   :: error_code
    integer                         :: err
    character(len=long_string_size) :: data_file
    character(len=long_string_size) :: chk_file
    character(len=long_string_size) :: cmd
    integer                         :: nargs
    character(len=long_string_size) :: arg
    integer                         :: iarg
    logical                         :: file_exists

    error_code = 0

    call get_command(cmd)
    nargs = command_argument_count()

    if (nargs == 0) then
       ! no args: print syntax and stop
       write(0,100)
       error_code = 1
       stop
    end if

    ! default 
    data_file = ''
    chk_file = ''

    iarg = 0
    args_loop: do while(iarg < nargs)
       iarg = iarg + 1
       call read_opt(iarg,nargs,arg,err)
       select case(trim(arg))
       case('-h','--help')
          ! help: print syntax and exit
          write(0,100)
          error_code = 1
          return
       case('-i','--int','--fasta')
          ! input file
          select case(trim(arg))
          case('-i','--int')
             data_type = 'int'
          case('--fasta')
             data_type = 'bio'
          end select
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,data_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('-c','--chk')
          ! chk file
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,chk_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('--prefix')
          ! prefix
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,prefix,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <prefix>'
             error_code = 1
             return
          end if
       case default
          write(0,*) 'ERROR ! invalid option '//trim(arg)
          error_code = 1
          return
       end select
    end do args_loop

    if ( chk_file /= "" ) then
       inquire( file = chk_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(0,*) 'ERROR ! cannot access '//trim(chk_file)
          error_code = 1
          return
       end if
       call units_open_unf(chk_file,'old',uchk,err)
       if( err /= 0 ) then
          write(0,*) 'ERROR ! error opening file '//trim(chk_file)
          error_code = 1
          return
       end if
    else
       write(0,*) 'ERROR ! chk file is needed'
       error_code = 1
       return
    end if

    if( data_file == '' ) then
       write(0,*) 'ERROR ! missing datafile'
       error_code = 1
       return
    else
       inquire( file = data_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(0,*) 'ERROR ! cannot access '//trim(data_file)
          error_code = 1
          return
       end if
       ! open data file
       call units_open(data_file,'old',udata,err)
       if( err /= 0 ) then
          write(0,*) 'ERROR ! error opening file '//trim(data_file)
          error_code = 1
          return
       end if
    end if

100 format(&
         'elss-eval (elss v0.4)                                             '/&
         '                                                                  '/&
         'Usage:                                                            '/&
         '    elss-eval [options] -c <chk_file> (-i <data_file>|--fasta <MSA_file>)'/&
         '                                                                  '/&
         'Description:                                                      '/&
         '    Given a set of parameters, compute the (non-normalized)       '/&
         '    probability of each sample in a data file. The output         '/&
         '    first column contains the log(probability) of the samples     '/&
         '    (upan additive constant that cancels out in probability ratios).'/&
         '                                                                  '/&
         'Required:                                                         '/&
         '-c, --chk <chk_file>                                              '/&
         '    Read the parameters of the model from checkpoint file.        '/&
         '                                                                  '/&
         '-i, --int <data_file>                                             '/&
         '    Path to data file.                                            '/&
         '                                                                  '/&
         '--fasta <data_file>                                               '/&
         '    Path to FASTA multiple sequence alignment.                    '/&
         '                                                                  '/&
         'Options:                                                          '/&
         '-h, --help                                                        '/&
         '    Display this help and exit.                                   '/&
         '                                                                  '/&
         '--prefix <str>                                                    '/&         
         '    Prefix of output files.                                       '/&
         '                                                                  '/&
         'Examples:                                                         '/&
         '                                                                  '/)



  end subroutine command_line_read

end module eval_command_line
