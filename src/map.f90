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

module map
  use kinds
  use constants
  use mpi_wrapper
  implicit none
  private

  public :: map_learn

  real(kflt) :: eps_map = 0.01_kflt
  real(kflt) :: data_energy,regularization_energy,&
       data_energy1,regularization_energy1

contains

  subroutine map_learn(nvars,nclasses,niter_agd,niter_gd,lambda,mc_nsweeps,beta,&
       nupdate,data_format,ulog,fdata,seq,seqs_table,prm,fmodel)
    use dump, only: dump_array_file,dump_prm_file
    integer,          intent(inout)  :: nvars,nclasses
    integer,          intent(in)     :: niter_agd, niter_gd
    real(kflt),       intent(in)     :: lambda
    integer,          intent(in)     :: mc_nsweeps
    real(kflt),       intent(in)     :: beta
    integer,          intent(in)     :: nupdate
    character(len=*), intent(in)     :: data_format  ! data format ('raw', 'table', 'protein')
    integer,          intent(in)     :: ulog
    real(kflt),       intent(in)     :: fdata(:)
    integer,          intent(inout)  :: seq(:)
    integer,          intent(inout)  :: seqs_table(:,:)
    real(kflt),       intent(inout)  :: prm(:)
    real(kflt),       intent(out)    :: fmodel(:)
    character(len=long_string_size) :: filename
    integer                         :: tot_iter
    integer                         :: dim1,dim2,dimm
    integer                         :: err,unt

    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2

    tot_iter = 0

    write(ulog,'(a)') &
         '# iter, '//&
         '    max_err,   mean_err,'//&
         '     ll_ene,   data_ene,'//&
         '  prior_ene,   cpu_time,       pacc'

    if (niter_agd > 0) call map_agd(nvars,nclasses,seq,seqs_table,prm,fmodel,&
         fdata,data_format,ulog,beta,lambda,&
         niter_agd,mc_nsweeps,tot_iter,nupdate)
    
    if (niter_gd > 0) call map_gd(nvars,nclasses,seq,seqs_table,prm,fmodel,&
         fdata,data_format,ulog,beta,lambda,&
         niter_gd,mc_nsweeps,tot_iter,nupdate)

    tot_iter = tot_iter + 1
    
    if(iproc == 0) then
       write(filename,*) tot_iter
       filename = 'prm'
!       call dump_array_file(filename,nvars,nclasses,&
!            prm(1:dim1),prm(dim1+1:dimm),err)
       call dump_prm_file(filename,data_format,nvars,nclasses,&
            prm(1:dim1),prm(dim1+1:dimm),err)
       if( err /= 0 ) then 
          call mpi_wrapper_finalize(err)
          stop
       end if
    end if

  end subroutine map_learn

  subroutine map_gd(nvars,nclasses,seq,seqs_table,prm,fmodel,fdata,data_format,ulog,beta,lambda,&
                    niter,mc_nsweeps,tot_iter,nupdate)
    use mcmc, only:       mcmc_update_energy
    use dump, only:       dump_rst
    use likelihood, only: likelihood_compute_energies,likelihood_compute_gradient
    integer,    intent(inout)    :: nvars,nclasses
    integer,    intent(inout)    :: seq(:)
    integer,    intent(inout)    :: seqs_table(:,:)
    real(kflt), intent(inout)    :: prm(:)
    real(kflt), intent(out)      :: fmodel(:)
    real(kflt), intent(in)       :: fdata(:)
    character(len=*), intent(in) :: data_format
    integer,    intent(in)       :: ulog
    real(kflt), intent(in)       :: beta
    real(kflt), intent(in)       :: lambda
    integer,    intent(in)       :: niter
    integer,    intent(in)       :: mc_nsweeps
    integer,    intent(inout)    :: tot_iter
    integer,    intent(in)       :: nupdate
    integer                         :: iter
    character(len=long_string_size) :: filename
    integer                         :: nind, err
    integer                         :: dim1, dim2, dimm
    real(kflt), allocatable         :: grd(:)
    real                            :: elapsed
    real(kflt)                      :: facc
    real(kflt)                      :: grd_nrm ! grd_nrm to be minimized; magnitude of the gradient
    real(kflt)                      :: max_err


    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2
    if (iproc == 0) allocate(grd(dimm),stat=err)

    do iter = 1,niter

       tot_iter = tot_iter + 1
       
       ! compute model frequencies
       call map_compute_fmodel(nvars,nclasses,seq,prm,beta,iter,mc_nsweeps,nupdate,fmodel,elapsed,facc)

       ! each chain knows the coordinates of other chains
       seqs_table(:,iproc+1) = seq
       CALL mpi_allgather(seq, nvars, MPI_INTEGER, seqs_table, nvars, MPI_INTEGER, MPI_COMM_WORLD, err)

       if(iproc == 0) then
          
          ! collect seqs from other procs
          ! dump prms
          if(mod(iter,1) == 0 .or. iter == 1) then
             call dump_rst('rst','replace',data_format,nvars,nclasses,nproc,seqs_table,prm,err)
             if( err /= 0 ) then 
                if ( iproc == 0 ) write(0,*) "error opening file rst", err
                call mpi_wrapper_finalize(err)
                stop
             end if
          end if

          ! print reconstr. err. 
          max_err = sqrt(maxval((fmodel - fdata - lambda * prm)**2))
          grd_nrm = sqrt(sum((fmodel - fdata - lambda * prm)**2)/size(prm))

          ! compute likelihood gradient 
          call likelihood_compute_gradient(fdata,fmodel,lambda,prm,grd)

          ! update parameters
          prm = prm - eps_map * grd

          call likelihood_compute_energies(fdata,prm,lambda,data_energy,regularization_energy)

          write(ulog,'(i6,1x,7f14.6)') &
               iter, max_err, grd_nrm, &
               (data_energy + regularization_energy)/float(nvars), & 
               data_energy/float(nvars), regularization_energy/float(nvars), elapsed, facc
          flush(ulog)

       end if

       call float_bcast(size(prm),prm)
       
       call mcmc_update_energy(nvars,nclasses,seq,prm(1:dim1),prm(dim1+1:dimm))

    end do

    if (iproc == 0) deallocate(grd)

  end subroutine map_gd

  subroutine map_agd(nvars,nclasses,seq,seqs_table,prm,fmodel,fdata,data_format,ulog,beta,lambda,niter,&
                     mc_nsweeps,tot_iter,nupdate)
    use mcmc, only:       mcmc_update_energy
    use dump, only:       dump_rst
    use likelihood, only: likelihood_compute_energies,likelihood_compute_gradient
    
    integer,    intent(inout)    :: nvars,nclasses
    integer,    intent(inout)    :: seq(:)
    integer,    intent(inout)    :: seqs_table(:,:)
    real(kflt), intent(inout)    :: prm(:)
    real(kflt), intent(out)      :: fmodel(:)
    real(kflt), intent(in)       :: fdata(:)
    character(len=*), intent(in) :: data_format
    integer,    intent(in)       :: ulog
    real(kflt), intent(in)       :: beta
    real(kflt), intent(in)       :: lambda
    integer,    intent(in)       :: niter
    integer,    intent(in)       :: mc_nsweeps
    integer,    intent(inout)    :: tot_iter
    integer,    intent(in)       :: nupdate
    integer                         :: iter
    character(len=long_string_size) :: filename
    integer                         :: nind, err
    integer                         :: dim1, dim2, dimm
    real(kflt), allocatable         :: grd(:)
    real(kflt), allocatable         :: prm1(:)
    real(kflt)                      :: tfista,tpfista,mnest
    real                            :: elapsed
    real(kflt)                      :: facc
    real(kflt)                      :: grd_nrm ! grd_nrm to be minimized; magnitude of the gradient
    real(kflt)                      :: max_err

    
    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2

    if (iproc == 0) then 
       allocate(prm1(size(prm)),stat=err)
       allocate(grd(size(prm)),stat=err)
       prm1 = prm
    end if

    tfista = 1.0_kflt

    do iter = 1,niter

       tot_iter = tot_iter + 1
       
       ! compute model frequencies
       call map_compute_fmodel(nvars,nclasses,seq,prm,beta,iter,mc_nsweeps,nupdate,fmodel,elapsed,facc)

       ! each chain knows the coordinates of other chains
       seqs_table(:,iproc+1) = seq
       CALL mpi_allgather(seq, nvars, MPI_INTEGER, seqs_table, nvars, MPI_INTEGER, MPI_COMM_WORLD, err)

       if(iproc == 0) then

          ! dump prms
          if(mod(iter,1) == 0 .or. iter == 1) then
             ! in AGD prm estimates are stored in prm1
             ! while averages are taken at prm 
             call dump_rst('rst','replace',data_format,nvars,nclasses,nproc,seqs_table,prm1,err)
             if( err /= 0 ) then 
                if ( iproc == 0 ) write(0,*) "error opening file rst", err
                call mpi_wrapper_finalize(err)
                stop
             end if
          end if
          
          ! print reconstr. err. 
          max_err = sqrt(maxval((fmodel - fdata - lambda * prm)**2))
          grd_nrm = sqrt(sum((fmodel - fdata - lambda * prm)**2)/size(prm))
          
          ! compute likelihood gradient 
          call likelihood_compute_gradient(fdata,fmodel,lambda,prm,grd)

          ! update parameters
          tpfista = 0.5_kflt * (1.0_kflt + sqrt(1.0_kflt + 4.0_kflt * tfista**2))
          mnest = (tfista - 1.0_kflt) / tpfista
          
          grd = - eps_map * grd + prm - prm1
          prm1 = prm1 + grd
          prm = prm1 + mnest * grd
          
          tfista = tpfista
          
          call likelihood_compute_energies(fdata,prm,lambda,data_energy,regularization_energy)

          call likelihood_compute_energies(fdata,prm1,lambda,data_energy1,regularization_energy1)

          write(ulog,'(i6,1x,7f14.6)') &
               iter, max_err, grd_nrm, &
               (data_energy1 + regularization_energy1)/float(nvars), & 
               data_energy1/float(nvars), regularization_energy1/float(nvars), elapsed, facc
          flush(ulog)

       end if

       call float_bcast(size(prm),prm)
       
       call mcmc_update_energy(nvars,nclasses,seq,prm(1:dim1),prm(dim1+1:dimm))

    end do

    if (iproc == 0) then 
       deallocate(prm1)
       deallocate(grd)
    end if

  end subroutine map_agd

  subroutine map_nag(nvars,nclasses,seq,seqs_table,prm,fmodel,fdata,data_format,ulog,beta,lambda,niter,&
       mc_nsweeps,tot_iter,nupdate)
    ! modified NAG algorithm (gradient is evaluated at the current parameters)
    use mcmc, only:       mcmc_update_energy
    use dump, only:       dump_rst
    use likelihood, only: likelihood_compute_energies,likelihood_compute_gradient
    
    integer,    intent(inout)    :: nvars,nclasses
    integer,    intent(inout)    :: seq(:)
    integer,    intent(inout)    :: seqs_table(:,:)
    real(kflt), intent(inout)    :: prm(:)
    real(kflt), intent(out)      :: fmodel(:)
    real(kflt), intent(in)       :: fdata(:)
    character(len=*), intent(in) :: data_format
    integer,    intent(in)       :: ulog
    real(kflt), intent(in)       :: beta
    real(kflt), intent(in)       :: lambda
    integer,    intent(in)       :: niter
    integer,    intent(in)       :: mc_nsweeps
    integer,    intent(inout)    :: tot_iter
    integer,    intent(in)       :: nupdate
    integer                         :: iter
    character(len=long_string_size) :: filename
    integer                         :: nind, err
    integer                         :: dim1, dim2, dimm
    real(kflt), allocatable         :: grd(:)
    real(kflt), allocatable         :: prm1(:)
    real(kflt)                      :: tfista,tpfista,mnest
    real                            :: elapsed
    real(kflt)                      :: facc
    real(kflt)                      :: grd_nrm ! grd_nrm to be minimized; magnitude of the gradient
    real(kflt)                      :: max_err

    
    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2

    if (iproc == 0) then 
       allocate(prm1(size(prm)),stat=err)
       allocate(grd(size(prm)),stat=err)
       prm1 = prm
    end if

    tfista = 1.0_kflt

    do iter = 1,niter

       tot_iter = tot_iter + 1
       
       ! compute model frequencies
       call map_compute_fmodel(nvars,nclasses,seq,prm,beta,iter,mc_nsweeps,nupdate,fmodel,elapsed,facc)

       ! each chain knows the coordinates of other chains
       seqs_table(:,iproc+1) = seq
       CALL mpi_allgather(seq, nvars, MPI_INTEGER, seqs_table, nvars, MPI_INTEGER, MPI_COMM_WORLD, err)

       if(iproc == 0) then

          ! dump prms
          if(mod(iter,1) == 0 .or. iter == 1) then
             ! in AGD prm estimates are stored in prm1
             ! while averages are taken at prm 
             call dump_rst('rst','replace',data_format,nvars,nclasses,nproc,seqs_table,prm,err)
             if( err /= 0 ) then 
                if ( iproc == 0 ) write(0,*) "error opening file rst", err
                call mpi_wrapper_finalize(err)
                stop
             end if
          end if
          
          ! print reconstr. err. 
          max_err = sqrt(maxval((fmodel - fdata - lambda * prm)**2))
          grd_nrm = sqrt(sum((fmodel - fdata - lambda * prm)**2)/size(prm))
          
          ! compute likelihood gradient 
          call likelihood_compute_gradient(fdata,fmodel,lambda,prm,grd)

          ! update parameters
          tpfista = 0.5_kflt * (1.0_kflt + sqrt(1.0_kflt + 4.0_kflt * tfista**2))
          mnest = (tfista - 1.0_kflt) / tpfista

          prm = prm + mnest**2 * prm1 - (1.0_kflt+mnest) * eps_map * grd
          prm1 = mnest * prm1 - eps_map * grd
          
          tfista = tpfista
          
          call likelihood_compute_energies(fdata,prm,lambda,data_energy,regularization_energy)

          call likelihood_compute_energies(fdata,prm1,lambda,data_energy1,regularization_energy1)

          write(ulog,'(i6,1x,7f14.6)') &
               iter, max_err, grd_nrm, &
               (data_energy + regularization_energy)/float(nvars), & 
               data_energy/float(nvars), regularization_energy/float(nvars), elapsed, facc
          flush(ulog)

       end if

       call float_bcast(size(prm),prm)
       
       call mcmc_update_energy(nvars,nclasses,seq,prm(1:dim1),prm(dim1+1:dimm))

    end do

    if (iproc == 0) then 
       deallocate(prm1)
       deallocate(grd)
    end if

  end subroutine map_nag

  subroutine map_adam(nvars,nclasses,seq,seqs_table,prm,fmodel,fdata,data_format,ulog,beta,lambda,niter,&
                     mc_nsweeps,tot_iter,nupdate)
    use mcmc, only:       mcmc_update_energy
    use dump, only:       dump_rst
    use likelihood, only: likelihood_compute_energies,likelihood_compute_gradient
    
    integer,    intent(inout)    :: nvars,nclasses
    integer,    intent(inout)    :: seq(:)
    integer,    intent(inout)    :: seqs_table(:,:)
    real(kflt), intent(inout)    :: prm(:)
    real(kflt), intent(out)      :: fmodel(:)
    real(kflt), intent(in)       :: fdata(:)
    character(len=*), intent(in) :: data_format
    integer,    intent(in)       :: ulog
    real(kflt), intent(in)       :: beta
    real(kflt), intent(in)       :: lambda
    integer,    intent(in)       :: niter
    integer,    intent(in)       :: mc_nsweeps
    integer,    intent(inout)    :: tot_iter
    integer,    intent(in)       :: nupdate
    integer                         :: iter
    character(len=long_string_size) :: filename
    integer                         :: nind, err
    integer                         :: dim1, dim2, dimm
    real(kflt), allocatable         :: grd(:)
    real(kflt), allocatable         :: prm1(:)
    real(kflt), allocatable         :: prm2(:)
    real(kflt)                      :: tfista,tpfista,mnest
    real                            :: elapsed
    real(kflt)                      :: facc
    real(kflt)                      :: grd_nrm ! grd_nrm to be minimized; magnitude of the gradient
    real(kflt)                      :: max_err
!    real(kflt)                      :: mu=0.9_kflt,nu=0.999_kflt
    real(kflt)                      :: mu=0.9_kflt,nu=0.999_kflt
!    real(kflt)                      :: mu=0.9_kflt,nu=0.9_kflt
!    real(kflt)                      :: mu=0.8_kflt,nu=0.9_kflt

    
    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2

    eps_map = 0.001_kflt

    if (iproc == 0) then 
       allocate(prm1(size(prm)),stat=err)
       allocate(prm2(size(prm)),stat=err)
       allocate(grd(size(prm)),stat=err)
       prm1 = 0.0_kflt
       prm2 = 0.0_kflt
    end if

    tfista = 1.0_kflt

    do iter = 1,niter

       tot_iter = tot_iter + 1
       
       ! compute model frequencies
       call map_compute_fmodel(nvars,nclasses,seq,prm,beta,iter,mc_nsweeps,nupdate,fmodel,elapsed,facc)

       ! each chain knows the coordinates of other chains
       seqs_table(:,iproc+1) = seq
       CALL mpi_allgather(seq, nvars, MPI_INTEGER, seqs_table, nvars, MPI_INTEGER, MPI_COMM_WORLD, err)

       if(iproc == 0) then

          ! dump prms
          if(mod(iter,1) == 0 .or. iter == 1) then
             ! in AGD prm estimates are stored in prm1
             ! while averages are taken at prm 
             call dump_rst('rst','replace',data_format,nvars,nclasses,nproc,seqs_table,prm,err)
             if( err /= 0 ) then 
                if ( iproc == 0 ) write(0,*) "error opening file rst", err
                call mpi_wrapper_finalize(err)
                stop
             end if
          end if
          
          ! print reconstr. err. 
          max_err = sqrt(maxval((fmodel - fdata - lambda * prm)**2))
          grd_nrm = sqrt(sum((fmodel - fdata - lambda * prm)**2)/size(prm))
          
          ! compute likelihood gradient 
          call likelihood_compute_gradient(fdata,fmodel,lambda,prm,grd)

          ! update parameters
          tpfista = 0.5_kflt * (1.0_kflt + sqrt(1.0_kflt + 4.0_kflt * tfista**2))
          mnest = (tfista - 1.0_kflt) / tpfista

          prm1 = mu*prm1 + (1.0_kflt - mu) * grd
          prm2 = nu*prm2 + (1.0_kflt - nu) * grd**2

          !          prm1 = prm1 / (1.0_kflt - mu**iter)
          !          prm2 = prm2 / (1.0_kflt - nu**iter)
          
!          prm = prm - (eps_map/sqrt(iter*1.0_kflt)) * (prm1 / (1.0_kflt - mu**iter)) / &
!               (sqrt(prm2 / (1.0_kflt - nu**iter)) + 1.e-8_kflt)

          prm = prm - eps_map * (prm1 / (1.0_kflt - mu**iter)) / &
               (sqrt(prm2 / (1.0_kflt - nu**iter)) + 1.e-8_kflt)

          
          tfista = tpfista
          
          call likelihood_compute_energies(fdata,prm,lambda,data_energy,regularization_energy)

          write(ulog,'(i6,1x,7f14.6)') &
               iter, max_err, grd_nrm, &
               (data_energy + regularization_energy)/float(nvars), & 
               data_energy/float(nvars), regularization_energy/float(nvars), elapsed, facc
          flush(ulog)

       end if

       call float_bcast(size(prm),prm)
       
       call mcmc_update_energy(nvars,nclasses,seq,prm(1:dim1),prm(dim1+1:dimm))

    end do

    if (iproc == 0) then 
       deallocate(prm1)
       deallocate(prm2)
       deallocate(grd)
    end if

  end subroutine map_adam

  subroutine map_compute_fmodel(nvars,nclasses,seq,prm,beta,iter,&
                                mc_nsweeps,nupdate,fmodel,elapsed,facc)
    ! given a set of prms, compute the model frequencies via MC simulation
    use mcmc, only: mcmc_simulate
    integer,    intent(in)     :: nvars,nclasses
    integer,    intent(inout)  :: seq(:)
    real(kflt), intent(in)     :: prm(:)
    real(kflt), intent(in)     :: beta
    integer,    intent(in)     :: iter
    integer,    intent(in)     :: mc_nsweeps
    integer,    intent(in)     :: nupdate
    real(kflt), intent(out)    :: fmodel(:)
    real,       intent(out)    :: elapsed
    real(kflt), intent(out)    :: facc    
    integer     :: dim1,dim2,err
    real        :: start,finish
    logical     :: hot_start


    dim1 = nclasses*nvars
    dim2 = nclasses**2*nvars*(nvars-1)/2

    call mpi_wrapper_barrier(err)
    call cpu_time(start)

    hot_start = .false. 
    call mcmc_simulate(nvars,nclasses,seq,&
         prm(1:dim1),prm(dim1+1:dim1+dim2),'raw',&
         fmodel(1:dim1),fmodel(dim1+1:dim1+dim2),&
         beta,mc_nsweeps,hot_start,nupdate,-1,facc)

    call mpi_wrapper_barrier(err)
    call cpu_time(finish)

    elapsed = finish - start

    call float_reduce(size(fmodel),fmodel)
    fmodel = fmodel / sum(fmodel(:nclasses))

  end subroutine map_compute_fmodel

end module map
