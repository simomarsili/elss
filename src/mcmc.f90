! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module mcmc
  use kinds
  use constants
  implicit none
  private

  public :: mcmc_update_energy
  public :: mcmc_compute_energy
  public :: mcmc_simulate

  public :: etot

  real(kflt) :: etot
  real(kflt) :: efields
  real(kflt) :: ecouplings

contains

  subroutine mcmc_simulate(nvars,nclasses,seq,fields,couplings,&
                           data_format,freq_single,freq_pair,beta,&
                           mc_nsweeps,hot_start,nupdate,utrj,facc)
    use dump, only:     dump_seq, dump_energies, dump_fasta
    use averages, only: averages_initialize,averages_update
    integer,          intent(in)    :: nvars,nclasses       ! nvars, nclasses
    integer,          intent(inout) :: seq(:)
    real(kflt),       intent(in)    :: fields(nclasses,nvars)
    real(kflt),       intent(in)    :: couplings(nclasses,nclasses,nvars*(nvars-1)/2)
    character(len=*), intent(in)    :: data_format  ! data format ('raw', 'table', 'FASTA')
    real(kflt),       intent(out)   :: freq_single(nclasses,nvars)
    real(kflt),       intent(out)   :: freq_pair(nclasses,nclasses,nvars*(nvars-1)/2)
    real(kflt),       intent(in)    :: beta                ! inverse temperature
    integer,          intent(in)    :: mc_nsweeps          ! total number of MC sweeps
    logical,          intent(in)    :: hot_start
    integer,          intent(in)    :: nupdate          ! stride for average update
    integer,          intent(in)    :: utrj
    real(kflt),       intent(out)   :: facc                ! p(acc)
    integer :: mc_step             ! current MCMC step
    integer :: mc_nsteps           ! total MCMC steps
    integer :: nmoves              ! stride for average update and
    integer :: na,nacc
    integer :: err

    !    if ( dump_trj ) then
    !       mc_nsteps = nvars * 1000000 !!!! PAPER VALUE
    !    end if

    ! initialize averages
    call averages_initialize(freq_single,freq_pair)

    ! initialize energy values (efields, ecouplings)
    call mcmc_update_energy(nvars,nclasses,seq,fields,couplings)

    ! set total mc steps
    mc_nsteps = mc_nsweeps*nvars

    ! set strides for averages update
    nmoves = nvars*nupdate

    ! mc_nsteps is the total number of mc steps for the run
    ! set mc_nsteps as a multiple of nmoves
    mc_nsteps = mc_nsteps - mod(mc_nsteps,nmoves)

    ! hot start: first mc_nsteps are discarded - no print at all
    !call mcmc_move(nvars,nclasses,seq,fields,couplings,beta,mc_nsteps,na)
    if (hot_start) call mcmc_move(nvars,nclasses,seq,fields,couplings,beta,mc_nsteps,na)

    mc_step = 0
    nacc = 0

    ! starting configuration
    if (utrj > 0) then 
       if (trim(data_format) == 'raw' .or. trim(data_format) == 'table') then 
          call dump_seq(utrj,seq,mc_step/nvars,etot,efields,ecouplings)
       else if (trim(data_format) == 'FASTA') then 
          call dump_fasta(utrj,seq,mc_step/nvars,etot,efields,ecouplings)
       end if
    end if

    mc_loop: do
       ! nmoves steps
       call mcmc_move(nvars,nclasses,seq,fields,couplings,beta,nmoves,na)
       mc_step = mc_step + nmoves
       nacc = nacc + na

       ! update averages every nmoves
       call averages_update(seq,freq_single,freq_pair)

       if (utrj > 0) then 
          if (trim(data_format) == 'raw' .or. trim(data_format) == 'table') then 
             call dump_seq(utrj,seq,mc_step/nvars,etot,efields,ecouplings)
          else if (trim(data_format) == 'FASTA') then 
             call dump_fasta(utrj,seq,mc_step/nvars,etot,efields,ecouplings)
          end if
       end if

       ! every mc_nsteps:
       if (mc_step == mc_nsteps) exit mc_loop

    end do mc_loop

    facc = real(nacc)/real(mc_step)

  end subroutine mcmc_simulate

  subroutine mcmc_compute_energy(nvars,nclasses,seq,fields,couplings,&
                                eh,ej,et)
    integer,    intent(in)  :: nvars,nclasses
    integer,    intent(in)  :: seq(:)
    real(kflt), intent(in)  :: fields(nclasses,nvars)
    real(kflt), intent(in)  :: couplings(nclasses,nclasses,nvars*(nvars-1)*nclasses**2)
    real(kflt), intent(out) :: eh,ej,et
    integer    :: iv,jv
    integer    :: is,js
    integer    :: k
    real(kflt) :: en
    
    eh = 0.0_kflt
    do iv = 1,nvars
       is = seq(iv)
       eh = eh + fields(is,iv)
    end do

    ej = 0.0_kflt
    k = 0
    do jv = 1,nvars-1
       js = seq(jv)
       do iv = jv+1,nvars
          is = seq(iv)
          k = k + 1
          ej = ej + couplings(is,js,k)
       end do
    end do

    eh = -eh
    ej = -ej
    et = eh + ej

  end subroutine mcmc_compute_energy

  subroutine mcmc_update_energy(nvars,nclasses,seq,fields,couplings)
    integer,    intent(in) :: nvars,nclasses
    integer,    intent(in) :: seq(:)
    real(kflt), intent(in) :: fields(nclasses,nvars)
    real(kflt), intent(in) :: couplings(nclasses,nclasses,nvars*(nvars-1)*nclasses**2)
    integer    :: iv,jv
    integer    :: is,js
    integer    :: k
    real(kflt) :: en
    
    efields = 0.0_kflt
    do iv = 1,nvars
       is = seq(iv)
       efields = efields + fields(is,iv)
    end do

    ecouplings = 0.0_kflt
    k = 0
    do jv = 1,nvars-1
       js = seq(jv)
       do iv = jv+1,nvars
          is = seq(iv)
          k = k + 1
          ecouplings = ecouplings + couplings(is,js,k)
       end do
    end do

    efields = -efields
    ecouplings = -ecouplings
    etot = efields + ecouplings

  end subroutine mcmc_update_energy

  subroutine mcmc_move(nvars,nclasses,seq,fields,couplings,beta,&
                       nsteps,nacc)
    integer,    intent(in)    :: nvars,nclasses
    integer,    intent(inout) :: seq(:)
    real(kflt), intent(in)    :: fields(:,:)
    real(kflt), intent(in)    :: couplings(:,:,:)
    real(kflt), intent(in)    :: beta
    integer,    intent(in)    :: nsteps
    integer,    intent(out)   :: nacc
    integer                 :: step
    real(kflt), allocatable :: rnd_array(:)
    integer                 :: i,j,k
    integer                 :: v,v2,s0,s1,s2
    integer                 :: iv,jv
    integer                 :: kr
    integer                 :: kv,ks0,ks1
    integer                 :: is,js
    real(kflt)              :: rnd
    real(kflt)              :: en1,en0,de,deh,dek
    real(kflt)              :: delta,boltz
    integer                 :: kstart,dk
    integer                 :: err

    allocate(rnd_array(nsteps*3),stat=err)
    call random_number(rnd_array)

    delta = 0.0_kflt
    kr = 0
    nacc = 0
    do step = 1,nsteps

       kr = kr + 1
       rnd = rnd_array(kr)
       v = int(rnd*nvars) + 1
       s0 = seq(v)

       kr = kr + 1
       rnd = rnd_array(kr)
       s1 = int(rnd*nclasses) + 1

       if ( s1 == s0 ) then
          nacc = nacc + 1
          cycle
       end if

       deh = 0.0_kflt
       dek = 0.0_kflt

       ! compute deh
       deh = fields(s1,v) - fields(s0,v)

       ! compute dek
       dek = 0.0_kflt
       ! j is the column index
       k = v - 1
       do j = 1,v-1
          js = seq(j) 
          dek = dek + couplings(s1,js,k) - couplings(s0,js,k)
          k = k + nvars - j - 1
       end do
       ! i is the row index
       do i = v+1,nvars
          is = seq(i) 
          k = k + 1
          dek = dek + couplings(is,s1,k) - couplings(is,s0,k)
       end do

       ! sum for de
       deh = -deh
       dek = -dek
       de = deh + dek

       ! metropolis
       if(de <= 0.0_kflt) then
          ! accept
          nacc = nacc + 1
          etot = etot + de
          efields = efields + deh
          ecouplings = ecouplings + dek
          seq(v) = s1
       else
          boltz = beta*de
          kr = kr + 1
          rnd = rnd_array(kr)
          if(rnd <= exp(-boltz)) then
             ! accept
             nacc = nacc + 1
             etot = etot + de
             efields = efields + deh
             ecouplings = ecouplings + dek
             seq(v) = s1
          end if
       end if

    end do

    deallocate(rnd_array)

  end subroutine mcmc_move

end module mcmc

!--------- a note on lower-triangular-packed storage mode (with no diagonal)

! in the code a lot of matrices are stored like this: 
! given a symmetric matrix A:
!
! a11 a12 a13 a14
! a21 a22 a23 a24
! a31 a32 a33 a34
! a41 a42 a43 a44
! we store the out-of-diagonal elements in a one-dimensional array: 
!
! a21 a31 a41 a32 a42 a43
! the element (i,j) correspond to the element k where: 
! k = (2n-j)*(j-1)/2+i-j (j < i)
!
! the elements can be accessed sequentially like this:
!    k = 0
!    do j = 1,n-1       ! loop on columns
!       js = seq(j)
!       do i = j+1,n  ! loop on rows
!          is = seq(i)
!          k = k + 1  ! k = (2n-j)*(j-1)/2+i-j
!       end do
!    end do


