! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module random
  use kinds
  implicit none
  private
  public :: random_initialize
  public :: random_seq
  
contains 

  subroutine random_initialize(rseed,pid)
    ! Random Numbers In Scientific Computing: An Introduction by Katzgrabber
    implicit none
    integer, intent(IN) :: rseed
    integer, intent(IN) :: pid
    integer              :: n
    integer,allocatable  :: seed(:)
    integer(kint_single) :: seedgen,clock
    integer              :: i
    real(kflt)           :: rnd
    
    call random_seed(size = n)
    allocate(seed(n))
    
    if (rseed == 0) then 
       call system_clock(clock)
    else
       clock = rseed
    end if
    seedgen = abs( mod((clock*181)*((pid-83)*359), 104729) )
    seed = seedgen + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)
    deallocate(seed)
    
    ! discard the first one
    call random_number(rnd)
    
  end subroutine random_initialize

  subroutine random_seq(nv,nc,seq)
    integer, intent(in)    :: nv,nc
    integer, intent(inout) :: seq(:)
    integer    :: err
    integer    :: iv,ic
    real(kflt) :: rnd

    do iv = 1,nv
       call random_number(rnd)
       ic = int(rnd*nc) + 1
       seq(iv) = ic
    end do

  end subroutine random_seq


  
end module random
