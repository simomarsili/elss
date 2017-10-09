! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module learn_command_line
  use kinds
  use constants
  implicit none
  private
  public :: command_line_read
  character(len=1), parameter                  :: nl=achar(10)
  character(len=long_string_size)              :: syntax = nl//& 
       '                                        elss-learn (v0.3.1)                                    '//nl//&
       '                                     =========================                                 '//nl//&
       nl//&
       nl//&
       'Option                         Description                                     (Default Value)  '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       ' --fasta <path_to_file>        MSA file (FASTA format)                         (None)           '//nl//&
       nl//&
       nl//&
       ' (-h|--help)                   print this help message                         (None)           '//nl//&
       nl//&
       ' --prefix <string>             prefix for output files                         (None)           '//nl//&
       nl//&
       ' (-p|--prm) <path_to_file>     parameters file                                 (None)           '//nl//&   
       '        OR'//nl//&   
       ' (-r|--rst) <path_to_file>     restart file                                    (None)           '//nl//&   
       nl//&
       ' (-n|--nsweeps) <int>          num. of MC sweeps per gradient estimate         (1000)           '//nl//&
       nl//&
       ' (-u|--nupdate) <int>          stride (as num. of sweeps) for averages updates (10)             '//nl//&
       nl//&
       ' --fasta <path_to_file>        MSA file (FASTA format)                         (None)           '//nl//&
       nl//&
       ' (-w|--weights) <path_to_file> weights file                                    (None)           '//nl//&
       '        OR'//nl//&   
       ' --wid <float>                 %id threshold for weights calculation           (-1)             '//nl//&
       nl//&
       ' --learn-gd <int>              num. of gradient descent steps                  (0)              '//nl//&
       nl//&
       ' (--learn|--learn-agd) <int>   num. of accelerated gradient descent steps      (2000)           '//nl//&
       nl//&
       ' (-l|--lambda) <float>         (scaled) regularization parameter               (0.01)           '//nl//&
       nl//&
       ' --random_seed <int>           init. the random seed (if 0, use system_clock)  (0)              '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       nl//&
       nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       ' For more information and examples, please check the project github repository:                 '//nl//&
       ' https://github.com/simomarsili/elss                                                            '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       '                                                                                                    '

contains

  subroutine command_line_read(udata,data_type,uwgt,wid,uprm,urst,&
                               rseed,beta,mc_nsweeps,nupdate,niter_agd,niter_gd,&
                               lambda,prefix,error_code)
    use units, only: units_open,units_open_unf
    use arguments, only: read_opt,read_opt_arg
    integer,                    intent(inout) :: udata
    character(len=*),           intent(inout) :: data_type
    integer,                    intent(inout) :: uwgt
    real(kflt),                 intent(inout) :: wid
    integer,                    intent(inout) :: uprm
    integer,                    intent(inout) :: urst
    integer,                    intent(inout) :: rseed
    real(kflt),                 intent(inout) :: beta
    integer,                    intent(inout) :: mc_nsweeps
    integer,                    intent(inout) :: nupdate
    integer,                    intent(inout) :: niter_agd
    integer,                    intent(inout) :: niter_gd
    real(kflt),                 intent(inout) :: lambda
    character(len=*),           intent(out) :: prefix
    integer,                    intent(out) :: error_code
    integer                         :: err
    character(len=long_string_size) :: data_file
    character(len=long_string_size) :: ww_file
    character(len=long_string_size) :: prm_file
    character(len=long_string_size) :: rst_file
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
       write(0,*) trim(syntax)
       error_code = 1
       return
    end if

    ! initialize filenames
    data_file = ''
    ww_file = ''
    prm_file = ''
    rst_file = ''
    
    iarg = 0
    args_loop: do while(iarg < nargs)
       iarg = iarg + 1
       call read_opt(iarg,nargs,arg,err)
       select case(trim(arg))
       case('-h','--help')
          ! print help and exit
          write(0,*) trim(syntax)
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
       case('-r','--rst')
          ! rst file
          iarg = iarg + 1
          call read_opt_arg(iarg,nargs,rst_file,err)
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

    if (prm_file /= "" .and. rst_file /= "") then
       write(0,*) 'ERROR ! either a rst or a prm file'
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

    if ( rst_file /= "" ) then
       inquire( file = rst_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(0,*) 'ERROR ! cannot access '//trim(rst_file)
          error_code = 1
          return
       end if
       call units_open_unf(rst_file,'old',urst,err)
       if( err /= 0 ) then
          write(0,*) 'ERROR ! error opening file '//trim(rst_file)
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

  end subroutine command_line_read

end module learn_command_line
