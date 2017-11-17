! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module averages
  use kinds
  implicit none
  private
  public :: update_single, update_pair, update_counts

contains

  subroutine update_single(w, seq, single)
    real(kflt), intent(in)    :: w
    integer,    intent(in)    :: seq(:)
    real(kflt), intent(inout) :: single(:, :)
    integer :: nv, iv, is

    nv = size(seq)

    do iv = 1, nv
       is = seq(iv)
       single(is, iv) = single(is, iv) + w
    end do

  end subroutine update_single

  subroutine update_pair(w, seq, pair)
    real(kflt), intent(in)    :: w
    integer,    intent(in)    :: seq(:)
    real(kflt), intent(inout) :: pair(:, :, :)
    integer :: iv, jv
    integer :: is, js
    integer :: k, nv

    nv = size(seq)

    k = 0
    do jv = 1, nv-1
       js = seq(jv)
       do iv = jv+1, nv
          is = seq(iv)
          k = k + 1
          pair(is, js, k) = pair(is, js, k) + w
       end do
    end do

  end subroutine update_pair

  subroutine update_counts(w, seq, single, pair)
    real(kflt), intent(in)    :: w
    integer,    intent(in)    :: seq(:)
    real(kflt), intent(inout) :: single(:, :)
    real(kflt), intent(inout) :: pair(:, :, :)
    integer :: iv, jv
    integer :: is, js
    integer :: k, nv

    nv = size(seq)

    do iv = 1, nv
       is = seq(iv)
       single(is, iv) = single(is, iv) + w
    end do

    k = 0
    do jv = 1, nv-1
       js = seq(jv)
       do iv = jv+1, nv
          is = seq(iv)
          k = k + 1
          pair(is, js, k) = pair(is, js, k) + w
       end do
    end do

  end subroutine update_counts

end module averages
