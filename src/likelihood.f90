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

module likelihood
  use kinds
  implicit none
  private 
  public :: likelihood_compute_gradient
  public :: likelihood_compute_energies
  
contains

  subroutine likelihood_compute_gradient(fdata,fmodel,lambda,prm,grd)
    real(kflt), intent(in)      :: fdata(:)
    real(kflt), intent(in)      :: fmodel(:)
    real(kflt), intent(in)      :: lambda
    real(kflt), intent(in)      :: prm(:)
    real(kflt), intent(out)     :: grd(:)
    
    ! the gradient is: F - f + l * p
    grd = fdata - fmodel + lambda * prm 
    
  end subroutine likelihood_compute_gradient

  subroutine likelihood_compute_energies(fdata,prm,lambda,edata,ereg)
    real(kflt), intent(in)  :: fdata(:)
    real(kflt), intent(in)  :: prm(:)
    real(kflt), intent(in)  :: lambda
    real(kflt), intent(out) :: edata,ereg

    edata = sum(fdata * prm)

    ereg = 0.5_kflt * lambda * sum(prm**2)

  end subroutine likelihood_compute_energies

end module likelihood
