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

  subroutine command_line_read(udata, data_type, uwgt, wid, uchk, rseed, &
       beta, mc_nsweeps, nupdate, algorithm, rate, niter, lambda,prefix, &
       error_code)
    use units, only: units_open, units_open_unf
    use arguments, only: read_opt, read_opt_arg
    integer,          intent(inout) :: udata
    character(len=*), intent(inout) :: data_type
    integer,          intent(inout) :: uwgt
    real(kflt),       intent(inout) :: wid
    integer,          intent(inout) :: uchk
    integer,          intent(inout) :: rseed
    real(kflt),       intent(inout) :: beta
    integer,          intent(inout) :: mc_nsweeps
    integer,          intent(inout) :: nupdate
    character(len=string_size), intent(inout) :: algorithm
    real(kflt),       intent(inout) :: rate
    integer,          intent(inout) :: niter
    real(kflt),       intent(inout) :: lambda
    character(len=*), intent(out) :: prefix
    integer,          intent(out) :: error_code
    integer                         :: err
    character(len=long_string_size) :: data_file
    character(len=long_string_size) :: ww_file
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
       case('-r','--rate')
          ! learning rate
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,rate,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check learning rate'
             error_code = 1
             return
          end if
          if ( rate <= 0.0 ) then
             write(0,*) 'ERROR ! learning rate must be > 0'
             error_code = 1
             return
          end if
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
       case('--algorithm')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,algorithm,err)
          if (err == 2) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <algorithm>'
             error_code = 1
             return
          end if
       case('--niter')
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,niter,err)
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
         'elss-learn (elss v0.5)                                            '/& 
         '                                                                  '/&
         'Usage:                                                            '/&
         '    elss-learn [options] (-i <file> | --fasta <file>)             '/&
         '    (prepend "mpiexec -n <n_proc>" for multiple parallel chains)  '/&
         '                                                                  '/&
         'Description:                                                      '/&
         '    Fit a model of pairwise-interacting categorical variables to  '/&
         '    input data. Valid formats are the FASTA format for biological '/&
         '    multiple sequence alignments, or plain space/tab separated    '/&
         '    integers. The fitted parameters are dumped to an unformatted  '/&
         '    checkpoint file `chk` that can be used as input to:           '/&
         '    - the elss-sample executable,                                 '/&
         '    to generate artificial data sampled from the fitted model     '/&
         '    - the elss-eval executable,                                   '/&
         '    to measure the estimated relative probability of new samples. '/&
         '    A matrix of scores measuring the coupling strength between    '/&
         '    different variables is printed to the scores file `scr`.      '/&
         '                                                                  '/&
         'Required:                                                         '/&
         '-i, --int <data_file>                                             '/&
         '    Path to data file (space/tab separated integers), or:         '/&
         '--fasta <data_file>                                               '/&
         '    Path to a multiple sequence alignment in FASTA format.        '/&
         '                                                                  '/&
         'Options:                                                          '/&
         '-h, --help                                                        '/&
         '    Display this help and exit.                                   '/&
         '                                                                  '/&
         '--prefix <str>                                                    '/&         
         '    Prefix of output files.                                       '/&
         '                                                                  '/&
         '-c, --chk <chk_file>                                              '/&
         '    Restart a previous calculation from checkpoint file.          '/&
         '                                                                  '/&
         '-n, --nsweeps <n>, integer                                        '/&
         '    Number of MC sweeps per gradient estimate.                    '/&
         '    [default: 1000]                                               '/&
         '                                                                  '/&
         '-u, --nupdate <n>, integer                                        '/&
         '    Gradient estimate: update model averages every <n> sweeps.    '/&
         '    [default: 10]                                                 '/&
         '                                                                  '/&
         '--wid <%id>, float                                                '/&
         '    Take averages reweighting for similarity of each sample to the'/&
         '    rest of data, using a threshold <%id> to define similarity.   '/&
         '    [default: no reweighting]                                     '/&
         '                                                                  '/&
         '-w, --weights <file>, path_to_file                                '/&
         '    Read weights for data from file <file>.                       '/&
         '                                                                  '/&
         '--algorithm <algorithm_name>, string                              '/&
         '    Options are: gd, momentum, nag, adam.                         '/&
         '    [default: adam]                                               '/&
         '                                                                  '/&
         '-r, --rate <learning_rate>, float                                 '/&
         '    Learning rate for the run.                                    '/&
         '    [default: 0.01]                                               '/&
         '                                                                  '/&
         '--niter <n_iter>, integer                                         '/&
         '    Number of iterations in the minimization process.             '/&
         '    [default: 2000]                                               '/&
         '                                                                  '/&
         '-l, --lambda <regularization_strength>, float                     '/&
         '    Parameter controlling the strength of l2 regularization.      '/&
         '    [default: 0.01]                                               '/&
         '                                                                  '/&
         '--seed <seed>, integer                                            '/&
         '    Seed for the pseudo-random number generator.                  '/&
         '    If == 0, the seed is generated from system clock.             '/&
         '    [default: 0]                                                  '/&
         '                                                                  '/&
         'Examples:                                                         '/&
         '                                                                  '/)

  end subroutine command_line_read

end module learn_command_line
