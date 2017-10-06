! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module cost
  use kinds
  implicit none
  private 
  public :: compute_gradient
  public :: compute_cost
  
contains

  subroutine compute_gradient(fdata,fmodel,lambda,prm,grd)
    real(kflt), intent(in)      :: fdata(:)
    real(kflt), intent(in)      :: fmodel(:)
    real(kflt), intent(in)      :: lambda
    real(kflt), intent(in)      :: prm(:)
    real(kflt), intent(out)     :: grd(:)
    
    ! the gradient is: F - f + l * p
    grd = fmodel - fdata + lambda * prm
    
  end subroutine compute_gradient

  subroutine compute_cost(fdata,prm,lambda,edata,ereg)
    real(kflt), intent(in)  :: fdata(:)
    real(kflt), intent(in)  :: prm(:)
    real(kflt), intent(in)  :: lambda
    real(kflt), intent(out) :: edata,ereg

    edata = -sum(fdata * prm)

    ereg = 0.5_kflt * lambda * sum(prm**2)

  end subroutine compute_cost

end module cost
