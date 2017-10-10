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

  subroutine command_line_read(uprm,uchk,useq,rseed,beta,mc_nsweeps,nupdate,&
                               prefix,error_code)
    use units, only: units_open,units_open_unf
    use arguments, only: read_opt,read_opt_arg
    integer,                    intent(inout) :: uprm
    integer,                    intent(inout) :: uchk
    integer,                    intent(inout) :: useq
    integer,                    intent(inout) :: rseed
    real(kflt),                 intent(inout) :: beta
    integer,                    intent(inout) :: mc_nsweeps
    integer,                    intent(inout) :: nupdate
    character(len=*),           intent(out) :: prefix
    integer,                    intent(out)   :: error_code
    integer                         :: err
    character(len=long_string_size) :: prm_file
    character(len=long_string_size) :: chk_file
    character(len=long_string_size) :: seq_file
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
    prm_file = ''
    chk_file = ''
    seq_file = ''

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
       case('-p','--prm')
          ! prm file
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,prm_file,err)
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
       case('-s','--seq')
          ! seq file
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,seq_file,err)
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

    if (prm_file /= "" .and. chk_file /= "") then
       write(0,*) 'ERROR ! either a chk or a prm file'
       error_code = 1
       return
    end if

    if (prm_file == "" .and. chk_file == "") then
       write(0,*) 'ERROR ! either a chk or a prm file'
       error_code = 1
       return
    end if

    if ( prm_file /= "" ) then
       inquire( file = prm_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(0,*) 'ERROR ! cannot access '//trim(prm_file)
          error_code = 1
          return
       end if
       call units_open(prm_file,'old',uprm,err)
       if( err /= 0 ) then
          write(0,*) 'ERROR ! error opening file '//trim(prm_file)
          error_code = 1
          return
       end if
    end if

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
    end if

    if ( seq_file /= "" ) then
       inquire( file = seq_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(0,*) 'ERROR ! cannot access '//trim(seq_file)
          error_code = 1
          return
       end if
       call units_open(seq_file,'old',useq,err)
       if( err /= 0 ) then
          write(0,*) 'ERROR ! error opening file '//trim(seq_file)
          error_code = 1
          return
       end if
    end if

100 format(&
         'elss-sample (elss v0.2.1)                                                      '/& 
         '                                                                               '/&
         'Usage:                                                                         '/&
         '    elss-sample [options] (-c <chk_file>|-p <prm_file>)                        '/&
         '                                                                               '/&
         'Description:                                                                   '/&
         '    Simulate a trajectory in sequence space, using a model of pairwise         '/&
         '    interacting variables. The parameters of the model can either be read from '/&
         '    a file of parameters or from a checkpoint file generated by elss-learn.    '/&
         '                                                                               '/&
         'Required:                                                                      '/&
         '-p, --prm <prm_file>                                                           '/&
         '    Read the parameters of the model from file <prm_file>.                     '/&
         '                                                                               '/&
         '-c, --chk <chk_file>                                                           '/&
         '    Read the parameters of the model and the starting point for the simulation '/&
         '    from checkpoint file <chk_file>.                                           '/&
         '                                                                               '/&
         'Options:                                                                       '/&
         '-h, --help                                                                     '/&
         '    Display this help and exit.                                                '/&
         '                                                                               '/&
         '--prefix <str>                                                                 '/&         
         '    Prefix of output files.                                                    '/&
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
         '    Seed for initialization of pseudo-random number generator.                 '/&
         '    If == 0, the seed is computed from the system clock.                       '/&
         '    [default: 0]                                                               '/&
         '                                                                               '/&
         '-s, --seq <file>, string                                                       '/&
         '    Read the starting sequence from file <file>.                               '/&
         '    Overrides the starting sequence read in checkpoint file.                   '/&
         '    [default: generate random starting sequence]                               '/&
         '                                                                               '/&
         'Examples:                                                                      '/&
         '                                                                               '/)

  end subroutine command_line_read

end module sample_command_line
