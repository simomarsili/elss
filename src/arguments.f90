! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module arguments
  use kinds
  use constants
  implicit none
  private
  public :: read_opt, read_opt_arg
  
  interface read_opt_arg
     module procedure read_arg_int
     module procedure read_arg_float_single
     module procedure read_arg_float_double
     module procedure read_arg_string
  end interface read_opt_arg

contains
  
  subroutine read_opt(iarg, nargs, arg, err)
    ! read an option from position iarg+1
    integer, intent(inout)        :: iarg
    integer, intent(in)           :: nargs
    character(len=*), intent(out) :: arg
    integer, intent(out)          :: err

    err = 0
    iarg = iarg + 1
    if (iarg > nargs) then
       err = 1 
       return
    end if
    call get_command_argument(iarg, arg)
    if (arg(1:1) /= '-') err = 2

  end subroutine read_opt

  subroutine read_arg_int(iarg, nargs, arg, err)
    ! read an integer from position iarg+1
    integer, intent(inout)                       :: iarg
    integer, intent(in)                          :: nargs
    integer, intent(out)                         :: arg
    integer, intent(out)                         :: err
    character(len=long_string_size) :: string

    err = 0
    iarg = iarg + 1
    if (iarg > nargs) then
       err = 1
       return
    end if
    call get_command_argument(iarg,string)
    if (string(1:1) == '-') then
       err = 2
       return
    end if
    read(string,*,iostat=err) arg
    if (err /= 0) then
       err = 3
       return
    end if

  end subroutine read_arg_int

  subroutine read_arg_float_single(iarg,nargs,arg,err)
    ! read a single precision float from position iarg+1
    integer, intent(inout)                       :: iarg
    integer, intent(in)                          :: nargs
    real(4), intent(out)                         :: arg
    integer, intent(out)                         :: err
    character(len=long_string_size) :: string

    err = 0
    iarg = iarg + 1
    if (iarg > nargs) then
       err = 1
       return
    end if
    call get_command_argument(iarg,string)
    if (string(1:1) == '-') then
       err = 2
       return
    end if
    read(string,*,iostat=err) arg
    if (err /= 0) then
       err = 3
       return
    end if
    
  end subroutine read_arg_float_single

  subroutine read_arg_float_double(iarg, nargs, arg, err)
    ! read a double precision float from position iarg+1
    integer, intent(inout)                       :: iarg
    integer, intent(in)                          :: nargs
    real(8), intent(out)                         :: arg
    integer, intent(out)                         :: err
    character(len=long_string_size)              :: string

    err = 0
    iarg = iarg + 1
    if (iarg > nargs) then
       err = 1
       return
    end if
    call get_command_argument(iarg, string)
    if (string(1:1) == '-') then
       err = 2
       return
    end if
    read(string,*,iostat=err) arg
    if (err /= 0) then
       err = 3
       return
    end if
    
  end subroutine read_arg_float_double

  subroutine read_arg_string(iarg, nargs, arg, err)
    ! read a string from position iarg+1
    integer, intent(inout)                       :: iarg
    integer, intent(in)                          :: nargs
    character(len=*), intent(out)                :: arg
    integer, intent(out)                         :: err
    character(len=long_string_size) :: string

    err = 0
    iarg = iarg + 1
    if (iarg > nargs) then
       err = 1
       return
    end if
    call get_command_argument(iarg, string)
    if (string(1:1) == '-') then
       err = 2
       return
    end if
    read(string,*,iostat=err) arg
    if (err /= 0) then
       err = 3
       return
    end if
    
  end subroutine read_arg_string
  
end module arguments
