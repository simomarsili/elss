! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module random
  use kinds
  implicit none
  private
  public :: random_initialize
  public :: random_data

  interface random_data
     module procedure random_data_1d
     module procedure random_data_2d
  end interface random_data
  
contains

  subroutine random_initialize(rseed, pid)
    ! Random Numbers In Scientific Computing: An Introduction by Katzgrabber
    implicit none
    integer, intent(in) :: rseed
    integer, intent(in) :: pid
    integer              :: i, n, seedgen, clock
    integer, allocatable :: seed(:)
    real(kflt)           :: rnd
    
    call random_seed(size = n)
    allocate(seed(n))
    if (rseed == 0) then 
       call system_clock(clock)
    else
       clock = rseed
    end if
    seedgen = abs(mod((clock * 181) * ((pid - 83) * 359), 104729))
    seed = seedgen + 37 * [(i - 1, i = 1, n)]
    call random_seed(PUT = seed)
    deallocate(seed)
    ! well discard the first one
    call random_number(rnd)
    
  end subroutine random_initialize

  subroutine random_data_1d(nc, data)
    integer, intent(in)    :: nc
    integer, intent(inout) :: data(:)
    integer    :: iv, ic
    real(kflt) :: rnd

    do iv = 1, size(data)
       call random_number(rnd)
       ic = int(rnd * nc) + 1
       data(iv) = ic
    end do
    
  end subroutine random_data_1d
  
  subroutine random_data_2d(nc, data)
    integer, intent(in)    :: nc
    integer, intent(inout) :: data(:, :)
    integer :: data_shape(2)
    integer    :: nr
    integer    :: iv, ic, ir
    real(kflt) :: rnd

    data_shape = shape(data)
    nr = data_shape(2)
    do ir = 1, nr
       call random_data_1d(nc, data(:, ir))
    end do
    
  end subroutine random_data_2d
  
end module random
