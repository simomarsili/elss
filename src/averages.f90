! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module averages
  use kinds
  implicit none
  private
  public :: averages_update

contains

  subroutine averages_update(seq, freq_single, freq_pair)
    integer,    intent(in)    :: seq(:)
    real(kflt), intent(inout) :: freq_single(:, :)
    real(kflt), intent(inout) :: freq_pair(:, :, :)
    integer :: iv, jv
    integer :: is, js
    integer :: k, nv

    nv = size(seq)

    do iv = 1, nv
       is = seq(iv)
       freq_single(is, iv) = freq_single(is, iv) + 1.0_kflt
    end do

    k = 0
    do jv = 1, nv-1
       js = seq(jv)
       do iv = jv+1, nv
          is = seq(iv)
          k = k + 1
          freq_pair(is, js, k) = freq_pair(is, js, k) + 1.0_kflt
       end do
    end do

  end subroutine averages_update

end module averages
