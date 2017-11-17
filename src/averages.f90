! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module averages
  use kinds
  implicit none
  private
  public :: update_counts

contains

  subroutine update_counts(w, seq, freq_single, freq_pair)
    real(kflt), intent(in)    :: w
    integer,    intent(in)    :: seq(:)
    real(kflt), intent(inout) :: freq_single(:, :)
    real(kflt), intent(inout) :: freq_pair(:, :, :)
    integer :: iv, jv
    integer :: is, js
    integer :: k, nv

    nv = size(seq)

    do iv = 1, nv
       is = seq(iv)
       freq_single(is, iv) = freq_single(is, iv) + w
    end do

    k = 0
    do jv = 1, nv-1
       js = seq(jv)
       do iv = jv+1, nv
          is = seq(iv)
          k = k + 1
          freq_pair(is, js, k) = freq_pair(is, js, k) + w
       end do
    end do

  end subroutine update_counts

end module averages
