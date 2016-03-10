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

module units
  use kinds
  implicit none
  private 
  public :: units_initialize
  public :: units_open, units_open_unf
  integer(kint) :: nunits

contains

  subroutine units_initialize
    implicit none
    
   nunits = 10 
    
  end subroutine units_initialize
  
  subroutine units_open(filename,status,fileunit,err)
    implicit none
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: status
    integer,          intent(out) :: fileunit
    integer,          intent(out) :: err
    
    nunits = nunits + 1
    fileunit = nunits

    open(unit=fileunit,file=filename,status=status,iostat=err)

  end subroutine units_open

  subroutine units_open_unf(filename,status,fileunit,err)
    implicit none
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: status
    integer,          intent(out) :: fileunit
    integer,          intent(out) :: err
    
    nunits = nunits + 1
    fileunit = nunits

    open(unit=fileunit,file=filename,status=status,access='stream',form='unformatted',iostat=err)

  end subroutine units_open_unf
  
end module units
