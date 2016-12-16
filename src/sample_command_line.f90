! Copyright (c) 2016 Simone Marsili
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in
! the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
! the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
! A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
! ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

module sample_command_line
  use kinds
  use constants
  implicit none
  private
  public :: command_line_read
  character(len=1), parameter                  :: nl=achar(10)
  character(len=long_string_size)              :: syntax = nl//& 
       '                                       mcdca-sample (v0.3.1)                                   '//nl//&
       '                                    ==========================                                 '//nl//&
       nl//&
       nl//&
       'Option                         Description                                     (Default Value)  '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       ' (-h|--help)                   print this help message                         (None)           '//nl//&
       nl//&
       ' (-p|--prm) <path_to_file>     parameters file                                 (None)           '//nl//&   
       '        OR'//nl//&   
       ' (-r|--rst) <path_to_file>     restart file                                    (None)           '//nl//&   
       nl//&
       ' (-n|--nsweeps) <int>          num. of MC sweeps                               (0)              '//nl//&
       nl//&
       ' (-u|--nupdate) <int>          stride (as num. of sweeps) for averages updates (10)             '//nl//&
       nl//&
       ' --fasta <path_to_file>        data file (MSA format)                          (None)           '//nl//&
       '        OR'//nl//&   
       ' --raw <path_to_file>          data file ("raw" format)                        (None)           '//nl//&
       !       '        OR'//nl//&   
       !       ' --table <path_to_file>        data file ("table" format)                      (None)           '//nl//&
       nl//&
       ' (-w|--weights) <path_to_file> weights file                                    (None)           '//nl//&
       nl//&
       ' (-s|--seq) <path_to_file>     starting sequence file (SIM)                    (None)           '//nl//&
       nl//&
       ' --wid <float>                 %id threshold for weights calculation           (-1)             '//nl//&
       nl//&
       ' --learn-gd <int>              num. of gradient descent steps                  (0)              '//nl//&
       nl//&
       ' (--learn|--learn-agd) <int>   num. of accelerated gradient descent steps      (0)              '//nl//&
       nl//&
       ' (-l|--lambda) <float>         (scaled) regularization parameter               (0.01)           '//nl//&
       nl//&
       ' --random_seed <int>           initialize the random seed                      (0)              '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       nl//&
       nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       ' For more information and usage examples, please check the project github repository:           '//nl//&
       ' https://github.com/simomarsili/mcDCA                                                           '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       '                                                                                                    '
contains

  subroutine command_line_read(uprm,urst,useq,rseed,beta,mc_nsweeps,nupdate,&
                               error_code)
    use units, only: units_open,units_open_unf
    use arguments, only: read_opt,read_arg
    integer,                    intent(inout) :: uprm
    integer,                    intent(inout) :: urst
    integer,                    intent(inout) :: useq
    integer,                    intent(inout) :: rseed
    real(kflt),                 intent(inout) :: beta
    integer,                    intent(inout) :: mc_nsweeps
    integer,                    intent(inout) :: nupdate
    integer,                    intent(out)   :: error_code
    integer                         :: err
    character(len=long_string_size) :: prm_file
    character(len=long_string_size) :: rst_file
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
       write(0,*) trim(syntax)
       error_code = 1
       return
    end if

    ! default 
    prm_file = ''
    rst_file = ''
    seq_file = ''

    iarg = 0
    args_loop: do while(iarg < nargs)
       call read_opt(iarg,nargs,arg,err)
       select case(trim(arg))
       case('-h','--help')
          ! print help and exit
          write(0,*) trim(syntax)
          error_code = 1
          return
       case('-p','--prm')
          ! prm file
          call read_arg(iarg,nargs,prm_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('-r','--rst')
          ! rst file
          call read_arg(iarg,nargs,rst_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('-s','--seq')
          ! seq file
          call read_arg(iarg,nargs,seq_file,err)
          if (err == 1) then
             write(0,*) 'ERROR ! missing argument: '//trim(arg)//' <filename>'
             error_code = 1
             return
          end if
       case('-t','--temp')
          call read_arg(iarg,nargs,beta,err)
          if ( beta < 0.0) then
             write(0,*) 'ERROR ! check beta'
             error_code = 1
             return
          end if
          beta = 1.0_kflt / beta 
       case('-n','--nsweeps')
          call read_arg(iarg,nargs,mc_nsweeps,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nsweeps'
             error_code = 1
             return
          end if
       case('--random_seed')
          call read_arg(iarg,nargs,rseed,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check random_seed'
             error_code = 1
             return
          end if
       case('-u','--nupdate')
          call read_arg(iarg,nargs,nupdate,err)
          if ( err/= 0 ) then
             write(0,*) 'ERROR ! check nupdate'
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

    if (prm_file == "" .and. rst_file == "") then
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

  end subroutine command_line_read

end module sample_command_line
