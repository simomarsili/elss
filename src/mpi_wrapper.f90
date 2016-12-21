! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module mpi_wrapper
  use kinds
!  use mpi
  implicit none

  include 'mpif.h'

  integer :: iproc
  integer :: nproc

contains

  subroutine mpi_wrapper_initialize(err)
    !   initialize the MPI environment
    implicit none
    ! local variables
    integer, intent(out) :: err

    CALL MPI_Init(err)

    CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, err)

    CALL MPI_Comm_rank(MPI_COMM_WORLD, iproc, err)

  end subroutine mpi_wrapper_initialize

  subroutine mpi_wrapper_finalize(err)
    implicit none
    integer, intent(out) :: err

    call mpi_finalize(err)

  end subroutine mpi_wrapper_finalize

  subroutine mpi_wrapper_barrier(err)
    implicit none
    integer, intent(out) :: err
    
    call mpi_barrier(MPI_COMM_WORLD,err)
    
  end subroutine mpi_wrapper_barrier

  subroutine float_allreduce(dim,a)
    integer,    intent(in)    :: dim
    real(kflt), intent(inout) :: a(dim)
    real(kflt) :: b(dim)
    integer    :: err

    select case(kflt)
    case(4)
       call MPI_ALLREDUCE(a,b,dim,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,err)
    case(8)
       call MPI_ALLREDUCE(a,b,dim,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,err)
    end select

    a = b 
    
  end subroutine float_allreduce

  subroutine int_allreduce(dim,a)
    integer, intent(in)    :: dim
    integer, intent(inout) :: a(dim)
    integer :: b(dim)
    integer :: err

    call MPI_ALLREDUCE(a,b,dim,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,err)

    a = b 
    
  end subroutine int_allreduce

  subroutine float_reduce(dim,a)
    integer, intent(in)       :: dim
    real(kflt), intent(inout) :: a(dim)
    real(kflt) :: b(dim)
    integer    :: err

    b = a

    select case(kflt)
    case(4)
       call MPI_REDUCE(b,a,dim,MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,err)
    case(8)
       call MPI_REDUCE(b,a,dim,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,err)
    end select
    
  end subroutine float_reduce

  subroutine int_reduce(dim,a)
    integer, intent(in)    :: dim
    integer, intent(inout) :: a(dim)
    integer :: b(dim)
    integer :: err

    b = a

    call MPI_REDUCE(b,a,dim,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,err)
    
  end subroutine int_reduce

  subroutine float_bcast(dim,a)
    integer, intent(in)       :: dim
    real(kflt), intent(inout) :: a(dim)
    integer :: err

    select case(kflt)
    case(4)
       call MPI_BCAST(a,dim,MPI_REAL4,0,MPI_COMM_WORLD,err)
    case(8)
       call MPI_BCAST(a,dim,MPI_REAL8,0,MPI_COMM_WORLD,err)
    end select
    
  end subroutine float_bcast

  subroutine int_bcast(dim,a)
    integer, intent(in)    :: dim
    integer, intent(inout) :: a(dim)
    integer :: err

    call MPI_BCAST(a,dim,MPI_INTEGER,0,MPI_COMM_WORLD,err)
    
  end subroutine int_bcast

end module mpi_wrapper
