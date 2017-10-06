! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module random
  use kinds
  implicit none
  private
  public :: random_initialize
  public :: random_seq
  
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

  subroutine random_seq(nv, nc, seq)
    integer, intent(in)    :: nv, nc
    integer, intent(inout) :: seq(:)
    integer    :: err
    integer    :: iv, ic
    real(kflt) :: rnd

    do iv = 1, nv
       call random_number(rnd)
       ic = int(rnd * nc) + 1
       seq(iv) = ic
    end do

  end subroutine random_seq


  
end module random
