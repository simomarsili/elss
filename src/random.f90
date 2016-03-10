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
