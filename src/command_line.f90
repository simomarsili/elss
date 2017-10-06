! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause
! License: BSD 3 clause

module command_line
  use kinds
  use constants
  implicit none
  private
  public :: command_line_read

contains

  subroutine command_line_read(udata,data_format,uwgt,wid,uprm,urst,useq,&
                               rseed,beta,mc_nsweeps,nupdate,niter_agd,niter_gd,&
                               lambda,mode,error_code,error_string)
    use units, only: units_open,units_open_unf
    integer,                    intent(inout) :: udata
    character(len=*),           intent(inout) :: data_format
    integer,                    intent(inout) :: uwgt
    real(kflt),                 intent(inout) :: wid
    integer,                    intent(inout) :: uprm
    integer,                    intent(inout) :: urst
    integer,                    intent(inout) :: useq
    integer,                    intent(inout) :: rseed
    real(kflt),                 intent(inout) :: beta
    integer,                    intent(inout) :: mc_nsweeps
    integer,                    intent(inout) :: nupdate
    integer,                    intent(inout) :: niter_agd
    integer,                    intent(inout) :: niter_gd
    real(kflt),                 intent(inout) :: lambda
    character(len=string_size), intent(inout) :: mode
    integer,                    intent(out) :: error_code
    character(len=*),           intent(out) :: error_string
    integer                         :: err
    character(len=long_string_size) :: data_file
    character(len=long_string_size) :: ww_file
    character(len=long_string_size) :: prm_file
    character(len=long_string_size) :: rst_file
    character(len=long_string_size) :: seq_file
    character(len=long_string_size) :: cmd
    integer                         :: nargs
    character(len=long_string_size) :: arg
    integer                         :: iarg
    logical                         :: file_exists

    error_code = 0
    error_string = ''

    call get_command(cmd)
    nargs = command_argument_count()

    ! call with no args
    if (nargs == 0) then
       error_code = -1
       return
    end if

    ! default 
    data_file = ''
    ww_file = ''
    prm_file = ''
    rst_file = ''
    seq_file = ''

    iarg = 1
    args_loop: do while(iarg <= nargs)
       call get_command_argument(iarg,arg)
       select case(trim(arg))
       case('-h','--help')
          ! print help and exit
          iarg = iarg + 1
          error_code = -1
          return
       case('-i','--raw','--fasta')
          ! input file
          if ( data_file /= "" ) then
             error_code = 16
             return
          end if
          select case(trim(arg))
          case('-i','--raw')
             data_format = 'raw'
          case('--fasta')
             data_format = 'FASTA'
          end select
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          data_file = arg
          if( data_file(1:1) == '-' ) then
             error_code = 4
             return
          end if
       case('-p','--prm')
          ! prm file
          if (prm_file /= "") then
             error_code = 26
             return
          end if
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          prm_file = arg
          if( prm_file(1:1) == '-' ) then
             error_code = 27
             return
          end if
       case('-r','--rst')
          ! rst file
          if (rst_file /= "") then
             error_code = 18
             return
          end if
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          rst_file = arg
          if( rst_file(1:1) == '-' ) then
             error_code = 6
             return
          end if
       case('-s','--seq')
          ! seq file
          if (seq_file /= "") then
             error_code = 40
             return
          end if
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          seq_file = arg
          if(seq_file(1:1) == '-') then
             error_code = 41
             return
          end if
       case('-w','--weights')
          ! input file
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          ww_file = arg
          if( ww_file(1:1) == '-' ) then
             error_code = 5
             return
          end if
       case('--wid')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) wid ! wid threshold for reweighting
          if ( wid < 0.0 .or. wid > 100.0 ) then
             error_code = 13
             write(error_string,*) trim(arg)
             return
          end if
       case('-t','--temp')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) beta ! temperature for the run 
          if ( beta < 0.0) then
             error_code = 45
             write(error_string,*) trim(arg)
             return
          end if
          beta = 1.0_kflt / beta 
       case('-n','--nsweeps')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) mc_nsweeps
          if ( err/= 0 ) then
             error_code = 10
             write(error_string,*) trim(arg)
             return
          end if
       case('--random_seed')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) rseed
          if ( err/= 0 ) then
             error_code = 44
             write(error_string,*) trim(arg)
             return
          end if
       case('-u','--nupdate')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) nupdate
          if ( err/= 0 .or. nupdate < 1) then
             error_code = 15
             write(error_string,*) trim(arg)
             return
          end if
       case('--learn','--learn-agd')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) niter_agd
          if ( err/= 0 ) then
             error_code = 11
             write(error_string,*) trim(arg)
             return
          end if
       case('--learn-gd')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) niter_gd
          if ( err/= 0 ) then
             error_code = 12
             write(error_string,*) trim(arg)
             return
          end if
       case('-l','--lambda')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) lambda ! 
          if ( lambda <= 0.0) then
             error_code = 38
             return
          end if
       case default
          error_code = 2
          write(error_string,*) trim(arg)
          return
       end select
       iarg = iarg + 1
    end do args_loop

    if (prm_file /= "" .and. rst_file /= "") then
       ! cannot read both rst and prm 
       error_code = 35
       return
    end if

    if ( prm_file /= "" ) then
       inquire( file = prm_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(error_string,*) trim(prm_file)
          error_code = 28
          return
       end if
       call units_open(prm_file,'old',uprm,err)
       if( err /= 0 ) then
          write(error_string,*) trim(prm_file)
          error_code = 19
          return
       end if
    end if

    if ( rst_file /= "" ) then
       inquire( file = rst_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(error_string,*) trim(rst_file)
          error_code = 9
          return
       end if
       call units_open_unf(rst_file,'old',urst,err)
       if( err /= 0 ) then
          write(error_string,*) trim(rst_file)
          error_code = 19
          return
       end if
    end if

    if ( seq_file /= "" ) then
       inquire( file = seq_file, exist = file_exists )
       if ( .not. file_exists ) then
          write(error_string,*) trim(seq_file)
          error_code = 42
          return
       end if
       call units_open(seq_file,'old',useq,err)
       if( err /= 0 ) then
          write(error_string,*) trim(seq_file)
          error_code = 19
          return
       end if
    end if


    !----------- run mode decision tree
    !
    !         has niter?
    !           /   \
    !      yes /     \ no
    !         /       \
    !       LEARN   has nsweeps?
    !                /   \
    !           yes /     \ no
    !              /       \
    !         SIM     EVAL
    !

    if (niter_gd > 0 .or. niter_agd > 0) then

       mode = 'LEARN'

       ! LEARN mode needs niter, nsweeps and data. prm/rst is optional

       if ( mc_nsweeps <= 0 ) then
          error_code = 10
          write(error_string,*) mc_nsweeps
          return
       end if

       if ( ww_file /= "" ) then
          inquire( file = ww_file, exist = file_exists )
          if ( .not. file_exists ) then
             error_code = 8
             write(error_string,*) trim(ww_file)
             return
          end if
          call units_open(ww_file,'old',uwgt,err)
          if( err /= 0 ) then
             error_code = 19
             write(error_string,*) trim(ww_file)
             return
          end if
       end if

       if ( uwgt > 0 .and. wid >= 0.0_kflt ) then
          error_code = 14
          return
       end if

    else

       if (mc_nsweeps > 0) then

          mode = 'SIM'

          ! simulation mode needs nsweeps and prm/rst. data is optional (wont read data). niter is incompatible.

       else

          mode = 'EVAL'

          ! EVAL mode needs data and prm/rst. nsweeps and niter are incompatible.

       end if

    end if

    if (uprm == 0 .and. urst == 0) then
       
       if (trim(mode) == 'EVAL' .or. trim(mode) == 'SIM') then 
          ! SIM / EVAL need at least one between prm and rst
          error_code = 36
          return
       end if
       
    end if

    if( data_file == '' ) then

       if (trim(mode) == 'EVAL' .or. trim(mode) == 'LEARN') then 
          error_code = 3
          return
       end if

    else

       inquire( file = data_file, exist = file_exists )
       if ( .not. file_exists ) then
          error_code = 7
          write(error_string,*) trim(data_file)
          return
       end if

       ! open data file
       call units_open(data_file,'old',udata,err)
       if( err /= 0 ) then
          error_code = 19
          write(error_string,*) trim(data_file)
          return
       end if

    end if

  end subroutine command_line_read

end module command_line
