! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module map
  use kinds
  use constants
  use mpi_wrapper
  implicit none
  private

  public :: map_learn

  real(kflt) :: data_energy, regularization_energy

contains

  subroutine map_learn(algorithm,rate,nvars,nclasses,niter,lambda,mc_nsweeps,beta,&
       nupdate,data_type,ulog,fdata,prefix,seq,seqs_table,prm,fmodel)
    use dump, only: dump_chk
    character(len=*), intent(in)  :: algorithm
    real(kflt), intent(in) :: rate
    integer,          intent(inout)  :: nvars,nclasses
    integer,          intent(in)     :: niter
    real(kflt),       intent(in)     :: lambda
    integer,          intent(in)     :: mc_nsweeps
    real(kflt),       intent(in)     :: beta
    integer,          intent(in)     :: nupdate
    character(len=*), intent(in)     :: data_type  ! data format ('unk', 'bio', 'protein', 'nuc_acid')
    integer,          intent(in)     :: ulog
    real(kflt),       intent(in)     :: fdata(:)
    character(len=*), intent(in)     :: prefix
    integer,          intent(inout)  :: seq(:)
    integer,          intent(inout)  :: seqs_table(:,:)
    real(kflt),       intent(inout)  :: prm(:)
    real(kflt),       intent(out)    :: fmodel(:)
    character(len=long_string_size) :: chk_file
    integer                         :: tot_iter
    integer                         :: dim1,dim2,dimm
    integer                         :: err

    if (len_trim(prefix) == 0) then
       chk_file = 'chk'
    else
       chk_file = trim(prefix)//'.chk'
    end if

    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2

    tot_iter = 0

    if (iproc == 0) then
       write(ulog,'(a)') &
            '# iter, '//&
            ' max_abs_grad, sum(grad**2),'//&
            '     data_ene, prior_energy,'//&
            '     cpu_time'
    end if

    call map_all(algorithm,rate,nvars,nclasses,seq,seqs_table,prm,fmodel,&
         fdata,data_type,ulog,beta,lambda,&
         niter,mc_nsweeps,tot_iter,nupdate)
    
    tot_iter = tot_iter + 1
    
    if(iproc == 0) then
       call dump_chk(chk_file,'replace',data_type,nclasses,&
            seqs_table,prm,err)
    end if

  end subroutine map_learn

  subroutine map_all(algorithm,eps_map,nvars,nclasses,seq,seqs_table,prm,fmodel,fdata,&
       data_type,ulog,beta,lambda,niter,&
       mc_nsweeps,tot_iter,nupdate)
    use mcmc, only:       mcmc_update_energy
    use dump, only:       dump_chk
    use cost, only: compute_cost,compute_gradient

    character(len=*), intent(in)  :: algorithm
    real(kflt), intent(in)       :: eps_map
    integer,    intent(inout)    :: nvars,nclasses
    integer,    intent(inout)    :: seq(:)
    integer,    intent(inout)    :: seqs_table(:,:)
    real(kflt), intent(inout)    :: prm(:)
    real(kflt), intent(out)      :: fmodel(:)
    real(kflt), intent(in)       :: fdata(:)
    character(len=*), intent(in) :: data_type
    integer,    intent(in)       :: ulog
    real(kflt), intent(in)       :: beta
    real(kflt), intent(in)       :: lambda
    integer,    intent(in)       :: niter
    integer,    intent(in)       :: mc_nsweeps
    integer,    intent(inout)    :: tot_iter
    integer,    intent(in)       :: nupdate
    integer                         :: iter
    integer                         :: err
    integer                         :: dim1, dim2, dimm
    real(kflt), allocatable         :: grd(:)
    real(kflt), allocatable         :: prm1(:)
    real(kflt), allocatable         :: prm2(:)
    real                            :: elapsed
    real(kflt)                      :: facc
    real(kflt),parameter            :: gamma1=0.9_kflt,gamma2=0.999_kflt
    real(kflt)                      :: g1,g2
    integer, parameter              :: nt0=10


    dim1 = nvars*nclasses
    dim2 = nvars*(nvars-1)*nclasses**2/2
    dimm = dim1 + dim2
    if (iproc == 0) allocate(grd(dimm),stat=err)
    
    ! initialize and set alg. dependent defaults
    dimm = size(prm)
    select case(trim(algorithm))
    case('momentum')
       if (iproc == 0) then
          allocate(prm1(dimm),stat=err)
          prm1 = prm
       end if
    case('nag')
       if (iproc == 0) then 
          allocate(prm1(dimm),stat=err)
          prm1 = prm
       end if
    case('adam')
       ! the "standard" value for adam should be eps_map = 1.e-3
       if (iproc == 0) then 
          allocate(prm1(dimm),stat=err)
          allocate(prm2(dimm),stat=err)
          prm1 = 0.0_kflt
          prm2 = 0.0_kflt
       end if
    case default
       if (iproc == 0) then 
          allocate(prm1(dimm),stat=err)
          allocate(prm2(dimm),stat=err)
          prm1 = 0.0_kflt
          prm2 = 0.0_kflt
       end if
    end select

    do iter = 1,niter

       tot_iter = tot_iter + 1
       
       ! compute model frequencies
       call map_compute_fmodel(nvars,nclasses,seq,prm,beta,iter,mc_nsweeps,nupdate,fmodel,elapsed,facc)

       ! each chain knows the coordinates of other chains
       seqs_table(:,iproc+1) = seq
       CALL mpi_allgather(seq, nvars, MPI_INTEGER, seqs_table, nvars, MPI_INTEGER, MPI_COMM_WORLD, err)

       if(iproc == 0) then
          if(mod(iter,1) == 0 .or. iter == 1) then
             call dump_chk('chk','replace',data_type,nclasses,seqs_table,prm,err)
             if( err /= 0 ) then 
                if ( iproc == 0 ) write(0,*) "error opening file chk", err
                call mpi_wrapper_finalize(err)
                stop
             end if
          end if
          
          ! compute gradient  of the cost function
          call compute_gradient(fdata,fmodel,lambda,prm,grd)
          
          write(ulog,'(i6,1x,5f14.6)') &
               iter, sqrt(maxval(grd**2)), sum(grd**2), &
               data_energy/float(nvars), regularization_energy/float(nvars), elapsed
          flush(ulog)

          ! update parameters and momentum terms
          ! see: https://github.com/lisa-lab/pylearn2/pull/136#issuecomment-10381617
          ! http://colinraffel.com/wiki/stochastic_optimization_techniques
          ! http://caffe.berkeleyvision.org/tutorial/solver.html
          ! https://arxiv.org/abs/1609.04747
          select case(trim(algorithm))
          case('gd')
             prm = prm - eps_map * grd
          case('momentum')
             ! this is not SGD with momentum. correct this. 
             prm1 = gamma1 * prm1 + (1.0_kflt -gamma1) * grd
             prm = prm - eps_map * prm1 / (1.0_kflt - gamma1**(iter+1))
          case('nag')
             ! same implementation as in lasagne v0.2
             ! velocity := momentum * velocity - learning_rate * gradient
             ! param := param + momentum * velocity - learning_rate * gradient
             ! """
             ! The classic formulation of Nesterov momentum (or Nesterov accelerated gradient)
             ! requires the gradient to be evaluated at the predicted next position in parameter space.
             ! Here, we use the formulation described at
             ! https://github.com/lisa-lab/pylearn2/pull/136#issuecomment-10381617,
             ! which allows the gradient to be evaluated at the current parameters (eqs 6,7).
             ! """
             prm1 = gamma1 * prm1 - eps_map * grd
             prm = prm + gamma1 * prm1 - eps_map * grd
          case('adam')
             ! update the (exp. decaying) average of the gradient
             prm1 = gamma1*prm1 + (1.0_kflt - gamma1) * grd
             ! update the (exp. decaying) average of the squared gradient
             prm2 = gamma2*prm2 + (1.0_kflt - gamma2) * grd**2
             ! update parameters using the 'unbiased' averages
             prm = prm - eps_map * (prm1 / (1.0_kflt - gamma1**(iter+1))) / &
                  (sqrt(prm2 / (1.0_kflt - gamma2**(iter+1))) + 1.e-8_kflt)
          case('adam-modified')
             if(real(iter) <= gamma1*nt0/(1.0_kflt-gamma1)) then
                g1 = real(iter) / real(iter+nt0)                
             else
                g1 = gamma1
             end if
             if(real(iter) <= gamma2*nt0/(1.0_kflt-gamma2)) then
                g2 = real(iter) / real(iter+nt0)                
             else
                g2 = gamma2
             end if
             prm1 = g1*prm1 + (1.0_kflt - g1) * grd
             prm2 = g2*prm2 + (1.0_kflt - g2) * grd**2
             prm = prm - eps_map * prm1 / (sqrt(prm2) + 1.e-8_kflt)
          case default
             ! default is adam 
             prm1 = gamma1*prm1 + (1.0_kflt - gamma1) * grd
             prm2 = gamma2*prm2 + (1.0_kflt - gamma2) * grd**2
             prm = prm - eps_map * (prm1 / (1.0_kflt - gamma1**(iter+1))) / &
                  (sqrt(prm2 / (1.0_kflt - gamma2**(iter+1))) + 1.e-8_kflt)
          end select
          
          call compute_cost(fdata,prm,lambda,data_energy,regularization_energy)

       end if

       call float_bcast(size(prm),prm)
       
       call mcmc_update_energy(nvars,nclasses,seq,prm(1:dim1),prm(dim1+1:dimm))

    end do

    if (iproc == 0) then
       deallocate(grd)
       select case(trim(algorithm))
       case('momentum')
          deallocate(prm1)
       case('nag')
          deallocate(prm1)
       case('adam')
          deallocate(prm1)
          deallocate(prm2)
       case default
          deallocate(prm1)
          deallocate(prm2)
       end select
    end if

  end subroutine map_all

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
