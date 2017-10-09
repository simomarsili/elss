! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module learn_command_line
  use kinds
  use constants
  implicit none
  private
  public :: command_line_read

contains

  subroutine command_line_read(udata, data_type, uwgt, wid, uprm, uchk, rseed, &
       beta, mc_nsweeps, nupdate, niter_agd, niter_gd, lambda,prefix, &
       error_code)
    use units, only: units_open, units_open_unf
    use arguments, only: read_opt, read_opt_arg
    integer,          intent(inout) :: udata
    character(len=*), intent(inout) :: data_type
    integer,          intent(inout) :: uwgt
    real(kflt),       intent(inout) :: wid
    integer,          intent(inout) :: uprm
    integer,          intent(inout) :: uchk
    integer,          intent(inout) :: rseed
    real(kflt),       intent(inout) :: beta
    integer,          intent(inout) :: mc_nsweeps
    integer,          intent(inout) :: nupdate
    integer,          intent(inout) :: niter_agd
    integer,          intent(inout) :: niter_gd
    real(kflt),       intent(inout) :: lambda
    character(len=*), intent(out) :: prefix
    integer,          intent(out) :: error_code
    integer                         :: err
    character(len=long_string_size) :: data_file
    character(len=long_string_size) :: ww_file
    character(len=long_string_size) :: prm_file
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
       write(0,100)! no args: print syntax and stop
       error_code = 1
       return
    end if

    ! initialize filenames
    data_file = ''
    ww_file = ''
    prm_file = ''
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
       case('-w','--weights')
          ! ww file
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,ww_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('--wid')
          ! identity threshold for reweighting
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,wid,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check wid threshold'
             error_code = 1
             return
          end if
          if ( wid < 0.0 .or. wid > 100.0 ) then
             write(0,*) 'ERROR ! wid must be > 0 and < 100'
             error_code = 1
             return
          end if
       case('-t','--temp')
          ! temperature for the run
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,beta,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check temperature'
             error_code = 1
             return
          end if
          if ( beta <= 0.0 ) then
             write(0,*) 'ERROR ! temperature must be > 0'
             error_code = 1
             return
          end if
          beta = 1.0_kflt / beta 
       case('-n','--nsweeps')
          ! num. of sweeps at each gradient evaluation
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,mc_nsweeps,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nsweeps'
             error_code = 1
             return
          end if
          if ( mc_nsweeps <= 0 ) then
             write(0,*) 'ERROR ! nsweeps must be > 0'
             error_code = 1
             return
          end if
       case('--random_seed')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,rseed,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check random_seed'
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
       case('--learn','--learn-agd')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,niter_agd,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nupdate'
             error_code = 1
             return
          end if
       case('--learn-gd')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,niter_gd,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nupdate'
             error_code = 1
             return
          end if
       case('-l','--lambda')
          ! regularization parameter
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,lambda,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check lambda'
             error_code = 1
             return
          end if
          if ( lambda <= 0.0 ) then
             write(0,*) 'ERROR ! lambda must be > 0'
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

    if ( ww_file /= "" ) then
       inquire( file = ww_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(0,*) 'ERROR ! cannot access '//trim(ww_file)
          error_code = 1
          return
       end if
       call units_open(ww_file,'old',uwgt,err)
       if( err /= 0 ) then
          write(0,*) 'ERROR ! error opening file '//trim(ww_file)
          error_code = 1
          return
       end if
    end if
    
    if ( uwgt > 0 .and. wid >= 0.0_kflt ) then
       write(0,*) 'ERROR ! computing and reading weights at the same time'
       error_code = 1
       return
    end if

    if( data_file == '' ) then

       write(0,*) 'ERROR ! a data file is needed'
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
         'elss-learn (elss v0.2.1)                                                       '/& 
         '                                                                               '/&
         'Usage:                                                                         '/&
         '    elss-learn [options] -i <data_file>                                        '/&
         '    elss-learn [options] --fasta <fasta_file>                                  '/&
         '    (preceded by mpiexec -n <n_proc> for parallel execution)                   '/&
         '                                                                               '/&
         'Description:                                                                   '/&
         '    Read a data matrix and fit a model of pairwise interacting categorical     '/&
         '    variables. Valid input formats are the FASTA format for biological multiple'/&
         '    sequence alignments, or plain space/tab separated values.                  '/&
         '    The fitted parameters are dumped in a `chk` file that can be an input to   '/&
         '    - elss-sample, to generate artificial samples according to the fitted model'/&
         '    - elss-eval, to measure the relative probability of new samples            '/&
         '      according to the model).                                                 '/&
         '                                                                               '/&
         'Required:                                                                      '/&
         '-i, --int <data_file>                                                          '/&
         '    Path to data file.                                                         '/&
         '                                                                               '/&
         '--fasta <data_file>                                                            '/&
         '    Path to FASTA multiple sequence alignment.                                 '/&
         '                                                                               '/&
         'Options:                                                                       '/&
         '-h, --help                                                                     '/&
         '    Display this help and exit.                                                '/&
         '                                                                               '/&
         '--prefix <str>                                                                 '/&         
         '    Prefix of output files.                                                    '/&
         '                                                                               '/&
         '-p, --prm <prm_file>                                                           '/&
         '    Read the initial guess for parameters from file <prm_file>.                '/&
         '                                                                               '/&
         '-c, --chk <chk_file>                                                           '/&
         '    Restart a previous calculation from checkpoint file <chk_file>             '/&
         '                                                                               '/&
         '-n, --nsweeps <n>, integer                                                     '/&
         '    Number of MC sweeps per gradient estimate.                                 '/&
         '    [default: 1000]                                                            '/&
         '                                                                               '/&
         '-u, --nupdate <n>, integer                                                     '/&
         '    During gradient estimate: update model averages every <n> sweeps.          '/&
         '    [default: 1000]                                                            '/&
         '                                                                               '/&
         '--wid <%id>, float                                                             '/&
         '    Compute data averages reweighting each sample according to its similarity  '/&
         '    to the rest of data, using a threshold <%id> to define data/sequence       '/&
         '    similarity.                                                                '/&
         '    [default: no reweighting]                                                  '/&
         '                                                                               '/&
         '-w, --weights <file>, path_to_file                                             '/&
         '    Read weights for data from file <file>.                                    '/&
         '                                                                               '/&
         '--learn, --learn-agd <n>, integer                                              '/&
         '    Number of accelerated gradient descent steps.                              '/&
         '    [default: 2000]                                                            '/&
         '                                                                               '/&
         '--learn, --learn-gd <n>, integer                                               '/&
         '    Number of gradient descent steps.                                          '/&
         '    [default: 0]                                                               '/&
         '                                                                               '/&
         '-l, --lambda <regularization_strength>, float                                  '/&
         '    Parameter controlling the strength of l2 regularization.                   '/&
         '    [default: 0.01]                                                            '/&
         '                                                                               '/&
         '--random_seed <seed>, integer                                                  '/&
         '    Seed for initialization of pseudo-random number generator.                 '/&
         '    If == 0, the seed is computed from the system clock.                       '/&
         '    [default: 0]                                                               '/)

  end subroutine command_line_read

end module learn_command_line
