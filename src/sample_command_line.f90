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

contains

  subroutine command_line_read(uprm,urst,useq,rseed,beta,mc_nsweeps,&
                               error_code,error_string)
    use units, only: units_open,units_open_unf
    integer,                    intent(inout) :: uprm
    integer,                    intent(inout) :: urst
    integer,                    intent(inout) :: useq
    integer,                    intent(inout) :: rseed
    real(kflt),                 intent(inout) :: beta
    integer,                    intent(inout) :: mc_nsweeps
    integer,                    intent(out)   :: error_code
    character(len=*),           intent(out)   :: error_string
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
    error_string = ''

    call get_command(cmd)
    nargs = command_argument_count()

    ! call with no args
    if (nargs == 0) then
       error_code = -1
       return
    end if

    ! default 
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
       case default
          error_code = 2
          write(error_string,*) trim(arg)
          return
       end select
       iarg = iarg + 1
    end do args_loop

    if (prm_file /= "" .and. rst_file /= "") then
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

    if (uprm == 0 .and. urst == 0) then
       
       ! need at least one between prm and rst
       error_code = 36
       return
       
    end if

  end subroutine command_line_read

end module sample_command_line