! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module sample_command_line
  use kinds
  use constants
  implicit none
  private
  public :: command_line_read

contains

  subroutine command_line_read(uchk,rseed,beta,istart,mc_nsweeps,nupdate,&
                               prefix,error_code)
    use units, only: units_open,units_open_unf
    use arguments, only: read_opt,read_opt_arg
    integer,                    intent(inout) :: uchk
    integer,                    intent(inout) :: rseed
    real(kflt),                 intent(inout) :: beta
    integer,                    intent(inout) :: istart
    integer,                    intent(inout) :: mc_nsweeps
    integer,                    intent(inout) :: nupdate
    character(len=*),           intent(out) :: prefix
    integer,                    intent(out)   :: error_code
    integer                         :: err
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
       return
    end if

    ! default 
    chk_file = ''

    iarg = 0
    args_loop: do while(iarg < nargs)
       iarg = iarg + 1
       call read_opt(iarg,nargs,arg,err)
       select case(trim(arg))
       case('-h','--help')
          ! print help and exit
          write(0,100)
          error_code = 1
          return
       case('-c','--chk')
          ! chk file
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,chk_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('-t','--temp')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,beta,err)
          if ( beta < 0.0) then
             write(0,*) 'ERROR ! check beta'
             error_code = 1
             return
          end if
          beta = 1.0_kflt / beta 
       case('-n','--nsweeps')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,mc_nsweeps,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nsweeps'
             error_code = 1
             return
          end if
       case('--seed')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,rseed,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check seed'
             error_code = 1
             return
          end if
       case('-u','--nupdate')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,nupdate,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nupdate'
             error_code = 1
             return
          end if
       case('--start')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,istart,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check --start argument'
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
       write(0,*) 'ERROR ! a chk file is needed'
       error_code = 1
       return
    end if

    if (istart < 0) then
       write(0,*) "ERROR ! index for starting point cant be < 0"
       error_code = 1
       return
    end if

100 format(&
         'elss-sample (elss v0.3)                                                      '/& 
         '                                                                               '/&
         'Usage:                                                                         '/&
         '    elss-sample [options] -c <chk_file>                                        '/&
         '                                                                               '/&
         'Description:                                                                   '/&
         '    Simulate a trajectory in "sequence" space using the parameters from        '/&
         '    checkpoint file <chk_file>. The unformatted checkpoint file can be         '/&
         '    generated using the elss-pchk tool from a set of user-defined parameters,  '/&
         '    or directly from data using the elss-learn executable.                     '/&
         '                                                                               '/&
         'Required:                                                                      '/&
         '-c, --chk <chk_file>                                                           '/&
         '    Read the parameters of the model from checkpoint file <chk_file>.          '/&
         '                                                                               '/&
         'Options:                                                                       '/&
         '-h, --help                                                                     '/&
         '    Display this help and exit.                                                '/&
         '                                                                               '/&
         '--prefix <str>                                                                 '/&         
         '    Prefix of output files.                                                    '/&
         '                                                                               '/&
         '--start <i>, integer                                                           '/&
         '    Use the <i>-th sample in <chk_file> as starting point for the MC chain.   '/&
         '    (indexing starts from zero)                                                '/&
         '    [default: 0]                                                               '/&
         '                                                                               '/&
         '-n, --nsweeps <n>, integer                                                     '/&
         '    Simulate <n> MC sweeps.                                                    '/&
         '    [default: 1000]                                                            '/&
         '                                                                               '/&
         '-u, --nupdate <n>, integer                                                     '/&
         '    Update model averages every <n> sweeps.                                    '/&
         '    [default: 10]                                                              '/&
         '                                                                               '/&
         '--seed <seed>, integer                                                         '/&
         '    Seed for the initialization of the pseudo-random number generator.         '/&
         '    If == 0, the seed is generated from system clock.                          '/&
         '    [default: 0]                                                               '/&
         '                                                                               '/&
         'Examples:                                                                      '/&
         '                                                                               '/)

  end subroutine command_line_read

end module sample_command_line
