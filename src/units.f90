! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

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
