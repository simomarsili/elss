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

module averages
  use kinds
  implicit none
  private
  public :: averages_initialize
  public :: averages_update

contains

  subroutine averages_initialize(freq_single,freq_pair)
    real(kflt), intent(out) :: freq_single(:,:)
    real(kflt), intent(out) :: freq_pair(:,:,:)

    freq_single = 0.0_kflt
    freq_pair = 0.0_kflt

  end subroutine averages_initialize

  subroutine averages_update(seq,freq_single,freq_pair)
    integer,    intent(in)    :: seq(:)
    real(kflt), intent(inout) :: freq_single(:,:)
    real(kflt), intent(inout) :: freq_pair(:,:,:)
    integer :: iv,jv
    integer :: is,js
    integer :: k,nv

    nv = size(seq)

    do iv = 1,nv
       is = seq(iv)
       freq_single(is,iv) = freq_single(is,iv) + 1.0_kflt
    end do

    k = 0
    do jv = 1,nv-1
       js = seq(jv)
       do iv = jv+1,nv
          is = seq(iv)
          k = k + 1
          freq_pair(is,js,k) = freq_pair(is,js,k) + 1.0_kflt
       end do
    end do

  end subroutine averages_update

end module averages
