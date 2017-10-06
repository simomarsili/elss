! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module scoring
  use kinds
  use constants
  implicit none
  private 

  public :: gauge
  public :: compute_scores
  public :: print_scores

contains 

  subroutine gauge(nvars,nclasses,fields,couplings)
    integer,    intent(in)    :: nvars,nclasses       ! nvars, nclasses
    real(kflt), intent(inout) :: fields(nclasses,nvars)
    real(kflt), intent(inout) :: couplings(nclasses,nclasses,nvars*(nvars-1)/2)
    real(kflt) :: mi(nclasses),mj(nclasses),mm
    integer    :: iv,jv
    integer    :: is,js
    integer    :: k

    k = 0
    do jv = 1,nvars-1
       do iv = jv+1,nvars
          k = k + 1
          
          mm = sum(couplings(:,:,k)) / real(nclasses**2)
          
          do is = 1,nclasses
             mi(is) = sum(couplings(is,:,k)) / real(nclasses)
             mj(is) = sum(couplings(:,is,k)) / real(nclasses)
          end do
          
          do js = 1,nclasses
             do is = 1,nclasses
                couplings(is,js,k) = & 
                     couplings(is,js,k) - mi(is) - mj(js) + mm
             end do
          end do
          
       end do
    end do
    
  end subroutine gauge

  subroutine compute_scores(nvars,nclasses,couplings,scores)
    integer,    intent(in)  :: nvars,nclasses       ! nvars, nclasses
    real(kflt), intent(in)  :: couplings(nclasses,nclasses,nvars*(nvars-1)/2)
    real(kflt), intent(out) :: scores(nvars,nvars)

    integer :: iv,jv
    integer :: is,js
    integer :: k

    scores = 0.0_kflt
    k = 0
    do jv = 1,nvars-1
       do iv = jv+1,nvars
          k = k + 1
          scores(iv,jv) = sqrt(sum(couplings(:,:,k)**2))
          scores(jv,iv) = scores(iv,jv)
       end do
    end do

    call apc_correction(nvars,scores)
    
  end subroutine compute_scores

  subroutine apc_correction(nvars,scores)
    integer,    intent(in)    :: nvars
    real(kflt), intent(inout) :: scores(nvars,nvars)
    real(kflt) :: mi(nvars),mj(nvars),mm
    integer    :: iv,jv
    integer    :: k

    mm = sum(scores) / real(nvars**2)
    do iv = 1,nvars
       mi(iv) = sum(scores(iv,:)) / real(nvars)
       mj(iv) = sum(scores(:,iv)) / real(nvars)
    end do

    do jv = 1,nvars
       do iv = jv,nvars
          scores(iv,jv) = scores(iv,jv) - mi(iv)*mj(jv)/mm
          scores(jv,iv) = scores(iv,jv)
       end do
    end do

  end subroutine apc_correction

  subroutine print_scores(nvars,scores)
    integer,    intent(in)    :: nvars
    real(kflt), intent(inout) :: scores(nvars,nvars)
    integer :: iv,jv
    
    do iv = 1,nvars-1
       do jv = iv+1,nvars
          write(444,*) iv,jv,scores(iv,jv)
       end do
    end do
    
  end subroutine print_scores


end module scoring
