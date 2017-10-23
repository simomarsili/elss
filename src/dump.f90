! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module dump
  use kinds
  use constants
  use units, only: units_open,units_open_unf
  implicit none
  private

  public :: read_chk
  public :: dump_chk
  public :: dump_seq
  public :: dump_energies

  interface read_chk
     module procedure read_chk_file
     module procedure read_chk_unit
  end interface read_chk

  interface dump_chk
     module procedure dump_chk_file
     module procedure dump_chk_unit
  end interface dump_chk

contains

  subroutine dump_energies(unt,etot,efields,ecouplings)
    integer, intent(in) :: unt
    real(kflt), intent(in) :: etot,efields,ecouplings

    write(unt,'(3(f12.3,1x))') etot,efields,ecouplings
    flush(unt)
    
  end subroutine dump_energies

  subroutine dump_seq(data_type, unt, seq, time, etot, eh, ej)
    character(len=*), intent(in) :: data_type
    integer, intent(in) :: unt
    integer, intent(in) :: seq(:)
    integer, intent(in) :: time
    real(kflt), intent(in) :: etot, eh, ej

    select case(trim(data_type))
    case ('int')
       call dump_int(unt, seq, time, etot, eh, ej)
    case ('bio', 'protein', 'nuc_acid')
       call dump_fasta(unt, seq, time, etot, eh, ej)
    end select

  end subroutine dump_seq

  subroutine dump_int(unt,seq,time,etot,eh,ej)
    integer, intent(in) :: unt
    integer, intent(in) :: seq(:)
    integer, intent(in) :: time
    real(kflt), intent(in) :: etot,eh,ej

    write(unt,'(a,1x,i8,3(1x,f12.3))') '>', time, etot, eh, ej
    write(unt,'(10000i3)') seq
    flush(unt)
    
  end subroutine dump_int

  subroutine dump_fasta(unt,seq,time,etot,eh,ej)
    use fasta, only: fasta_alphabet
    use constants
    integer, intent(in) :: unt
    integer, intent(in) :: seq(:)
    integer, intent(in) :: time 
    real(kflt), intent(in) :: etot,eh,ej
    integer :: i,n,si
    character(len=long_string_size) :: string

    n = size(seq)
    string=''
    do i = 1,n
       si = seq(i)
       string = trim(string)//fasta_alphabet(si:si)
    end do
    write(unt,'(a,1x,i8,3(1x,f12.3))') '>', time, etot, eh, ej
    write(unt,'(a)') trim(string)
    flush(unt)
    
  end subroutine dump_fasta

  subroutine read_chk_unit(unt,data_type,nvars,nclasses,iproc,nproc,seq,prm,error_code)
    ! read a checkpoint file 
    use random, only: random_seq
    use fasta, only: set_fasta_alphabet
    integer,          intent(in)                 :: unt
    character(len=*), intent(inout)              :: data_type
    integer,          intent(inout)              :: nvars,nclasses
    integer,          intent(in)                 :: iproc
    integer,          intent(in)                 :: nproc
    integer,          intent(inout), allocatable :: seq(:)
    real(kflt),       intent(inout), allocatable :: prm(:)
    integer,          intent(out) :: error_code
    integer, allocatable    :: dummy(:)
    integer                 :: p,np,err
    integer                 :: nv,nc

    error_code = 0

    read(unt) data_type
    read(unt) nv
    read(unt) nc
    read(unt) np

    if (nvars == 0 .and. nclasses == 0) then 
       nvars = nv
       nclasses = nc
       ! no data have been read
       ! => allocate and initialize
       if (allocated(seq) .or. allocated(prm)) then 
          ! arrays are already allocated
          ! => exit
          write(0,*) 'ERROR ! arrays are already allocated'
          error_code = 1
          return
       end if
       allocate(seq(nvars),stat=err)
       allocate(prm(nvars*nclasses + nvars*(nvars - 1)*nclasses**2/2),stat=err)
       seq = 0 
       prm = 0.0_kflt
    else 
       ! after reading data
       if (nv /= nvars .or. nc /= nclasses ) then 
          ! n. of variables/classes are inconsistent with data
          ! => exit
          write(0,*) 'ERROR ! n. of variables/classes is inconsistent with data'
          error_code = 1
          return
       end if
    end if

    allocate(dummy(nvars),stat=err)
    if (np == 0) then 
       ! no seqs found in checkpoint 
       ! => randomize 
       call random_seq(nvars,nclasses,seq)
    else
       if (np >= nproc) then
          ! if there are more samples than chains, read the starting points from the chk file
          do p = 1,np
             read(unt) dummy
             if (p == iproc+1) seq = dummy
          end do
       else
          ! else, randomize
          call random_seq(nvars,nclasses,seq)
       end if
    end if
          
    read(unt) prm

    deallocate(dummy)

    select case (trim(data_type))
    case('bio', 'protein', 'nuc_acid')
       call set_fasta_alphabet(data_type)
    case default
       ! do nothing
    end select
    
  end subroutine read_chk_unit

  subroutine read_chk_file(filename,data_type,nvars,nclasses,iproc,nproc,seq,&
                           prm,error_code)
    ! should read both a filename or a unit
    use random, only: random_seq
    character(len=*), intent(in)                 :: filename
    character(len=*), intent(inout)              :: data_type
    integer,          intent(inout)              :: nvars,nclasses
    integer,          intent(in)                 :: iproc
    integer,          intent(in)                 :: nproc
    integer,          intent(inout), allocatable :: seq(:)
    real(kflt),       intent(inout), allocatable :: prm(:)
    integer,          intent(out)                :: error_code
    integer :: unt,err

    error_code = 0

    call units_open_unf(filename,'old',unt,err)
    if( err /= 0 ) then
       write(0,*) 'ERROR ! error opening file '//trim(filename)
       error_code = 1
       return
    end if

    call read_chk_unit(unt,data_type,nvars,nclasses,iproc,nproc,seq,prm,error_code)
    
    close(unt)
    
  end subroutine read_chk_file

  subroutine dump_chk_unit(unt,data_type,nvars,nclasses,nproc,seqs_table,prm,&
                           error_code)
    ! dump a checkpoint file 
    integer,          intent(in)  :: unt
    character(len=*), intent(in)  :: data_type
    integer,          intent(in)  :: nvars,nclasses,nproc
    integer,          intent(in)  :: seqs_table(:,:)
    real(kflt),       intent(in)  :: prm(:)
    integer,          intent(out) :: error_code
    integer :: p

    error_code = 0

    write(unt) data_type
    write(unt) nvars
    write(unt) nclasses
    write(unt) nproc
    do p = 1,nproc
       write(unt) seqs_table(:,p)
    end do
    write(unt) prm

  end subroutine dump_chk_unit

  subroutine dump_chk_file(filename,status,data_type,nvars,nclasses,nproc,&
                           seqs_table,prm,error_code)
    ! dump a checkpoint file 
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: status
    character(len=*), intent(in)  :: data_type
    integer,          intent(in)  :: nvars,nclasses,nproc
    integer,          intent(in)  :: seqs_table(:,:)
    real(kflt),       intent(in)  :: prm(:)
    integer,          intent(out) :: error_code
    integer :: unt,err

    error_code = 0

    call units_open_unf(filename,status,unt,err)
    if( err /= 0 ) then
       write(0,*) 'ERROR ! error opening file '//trim(filename)
       error_code = 1
       return
    end if

    call dump_chk_unit(unt,data_type,nvars,nclasses,nproc,seqs_table,prm,error_code)

    close(unt)

  end subroutine dump_chk_file

end module dump
