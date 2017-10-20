! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module dump
  use kinds
  use constants
  use units, only: units_open,units_open_unf
  implicit none
  private

  public :: read_prm_unit
  public :: read_array_file,read_array_unit
  public :: dump_array_file,dump_array_unit
  public :: dump_prm_file,dump_prm_unit
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

  subroutine read_prm_unit(unt,nvars,nclasses,&
                                    prm,data_type,error_code)
    use parser, only: parser_nfields,remove_comments
    use random, only: random_seq
    use fasta, only: set_fasta_alphabet
    integer,    intent(in)                 :: unt
    integer,    intent(inout)              :: nvars, nclasses
    real(kflt), intent(inout), allocatable :: prm(:)
    character(len=*), intent(inout)        :: data_type
    integer,    intent(out)                :: error_code
    integer                         :: iv,jv
    integer                         :: is,js
    integer                         :: k
    integer                         :: index
    real(kflt)                      :: x
    character(len=long_string_size) :: line, newline
    integer                         :: nfields
    integer                         :: err
    integer                         :: nv,nc

    error_code = 0

    do 
       read(unt,'(a)') line
       call remove_comments(line)
       if (len_trim(line) == 0) then 
          cycle 
       else
          read(line,*) data_type, nv, nc
          exit
       end if
    end do

    select case(trim(data_type))
    case('protein')
       if (nc /= 21) then
          write(0,*) 'ERROR ! in prm file, num. of classes must be 21'
          error_code = 1
          return
       end if
    case('nuc_acid')
       if (nc /= 6) then
          write(0,*) 'ERROR ! in prm file, num. of classes must be 6'
          error_code = 1
          return
       end if
    end select
       
    if (nvars == 0 .and. nclasses == 0) then 
       nvars = nv
       nclasses = nc
       ! no data have been read
       ! => allocate and initialize
       if (allocated(prm)) then 
          ! arrays are already allocated
          ! => exit
          write(0,*) 'ERROR ! prm is already allocated'
          error_code = 1
          return
       end if
       allocate(prm(nvars*nclasses + nvars*(nvars-1)*nclasses**2/2),stat=err)
       prm = 0.0_kflt
    else 
       ! after reading data
       if (nv /= nvars .or. nc /= nclasses ) then 
          ! n. of variables/classes are inconsistent with data
          ! => exit
          write(0,*) 'ERROR ! inconsistent dimensionality nvars'
          error_code = 1
          return
       end if
    end if

    do 
       read(unt,'(a)',iostat=err) line
       if (err < 0) exit
       call remove_comments(line)
       if (len_trim(line) == 0) cycle
       call parser_nfields(line,newline,nfields)
       select case(nfields)
       case(3)
          read(line,*,iostat=err) iv,is,x
          if (err < 0) exit
          if (err > 0) then
             write(0,*) 'ERROR ! reading fields in prm file'
             error_code = 1
          end if
          index = (iv-1)*nclasses + is
          prm(index) = x
       case(5)
          ! k = (2n-j)*(j-1)+i-j
          read(line,*,iostat=err) iv,jv,is,js,x
          if (err < 0) exit
          if (err > 0) then
             write(0,*) 'ERROR ! reading couplings in prm file'
             error_code = 1
             return
          end if
          if (jv > iv) then 
             k = jv
             jv = iv
             iv = k 
             k = js
             js = is
             is = k 
          end if
          k = (2*nvars - jv) * (jv - 1)/2 + iv - jv
          index = nvars*nclasses + (k-1)*nclasses**2 + (js-1)*nclasses + is
          prm(index) = x
       case default
          write(0,*) 'ERROR ! inconsistent num. of lines in prm file'
          error_code = 1
          return
       end select
    end do

    select case (trim(data_type))
    case('bio', 'protein', 'nuc_acid')
       call set_fasta_alphabet(data_type)
    case default
       ! do nothing
    end select

  end subroutine read_prm_unit

  subroutine dump_prm_file(filename,data_type,nvars,nclasses,nproc,seqs,&
                       array1,array2,error_code)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: data_type
    integer,          intent(in)  :: nvars,nclasses,nproc
    integer,          intent(in)  :: seqs(:,:)
    real(kflt),       intent(in)  :: array1(nclasses,nvars)
    real(kflt),       intent(in)  :: array2(nclasses,nclasses,nvars*(nvars-1)/2)
    integer,          intent(out) :: error_code
    integer :: unt,err

    error_code = 0

    call units_open(filename,'unknown',unt,err)

    if( err /= 0 ) then
       write(0,*) 'ERROR ! error opening file '//trim(filename)
       error_code = 1
       return
    end if

    call dump_prm_unit(unt,data_type,nvars,nclasses,nproc,seqs,&
                       array1,array2,error_code)

    close(unt)

  end subroutine dump_prm_file

  subroutine dump_prm_unit(unt,data_type,nvars,nclasses,nproc,seqs,&
       array1,array2,error_code)
    ! dump a checkpoint file 
    integer,          intent(in)  :: unt
    character(len=*), intent(in)  :: data_type
    integer,          intent(in)  :: nvars,nclasses,nproc
    integer,          intent(in)  :: seqs(:,:)
    real(kflt),       intent(in)  :: array1(nclasses,nvars)
    real(kflt),       intent(in)  :: array2(nclasses,nclasses,nvars*(nvars-1)/2)
    integer,          intent(out) :: error_code
    integer :: iv,jv,is,js,k,p
    
    error_code = 0
    write(unt,'(a)') '# <data_type> <nvars> <nclasses> <nseq> <nproc>'
    write(unt,'(a,3(1x,i4))') trim(data_type), nvars, nclasses, nproc
    if (nproc > 0) then
       do p = 1,nproc
          write(unt,'(a,1000(1x,i2))') '# ',seqs(:,p)
       end do
    end if
    do iv = 1,nvars
       write(unt,*) iv,array1(:,iv)
    end do
    k = 0
    do jv = 1,nvars-1
       do iv = jv+1,nvars
          k = k + 1
          write(unt,*) jv, iv, array2(:,:,k)
       end do
    end do

  end subroutine dump_prm_unit

  subroutine read_array_unit(unt,nvars,nclasses,array1,array2,error_code)
    integer,    intent(in)    :: unt
    integer,    intent(inout) :: nvars, nclasses
    real(kflt), intent(inout) :: array1(nclasses,nvars)
    real(kflt), intent(inout) :: array2(nclasses,nclasses,nvars*(nvars-1)/2)
    integer,    intent(out)   :: error_code
    integer :: iv,jv
    integer :: is,js
    integer :: k,a,b,c,d

    error_code = 0

    read(unt,*) nvars,nclasses

    do iv = 1,nvars
       do is = 1,nclasses
          read(unt,*) a,b,array1(is,iv)
       end do
    end do

    k = 0
    do jv = 1,nvars-1
       do iv = jv+1,nvars
          k = k + 1
          do is = 1,nclasses
             do js = 1,nclasses
                read(unt,*) a,b,c,d,array2(is,js,k)
             end do
          end do
       end do
    end do

  end subroutine read_array_unit

  subroutine read_array_file(filename,nvars,nclasses,array1,array2,&
                             error_code)
    character(len=*), intent(in)    :: filename
    integer,          intent(inout) :: nvars, nclasses
    real(kflt),       intent(inout) :: array1(nclasses,nvars)
    real(kflt),       intent(inout) :: array2(nclasses,nclasses,nvars*(nvars-1)/2)
    integer,          intent(out)   :: error_code
    integer :: unt,err

    error_code = 0

    call units_open(filename,'unknown',unt,err)
    if( err /= 0 ) then
       write(0,*) 'ERROR ! error opening file '//trim(filename)
       error_code = 1
       return
    end if

    call read_array_unit(unt,nvars,nclasses,array1,array2,err)

    close(unt)

  end subroutine read_array_file

  subroutine dump_array_unit(unt,nvars,nclasses,array1,array2,error_code)
    integer,    intent(in)  :: unt
    integer,    intent(in)  :: nvars, nclasses
    real(kflt), intent(in)  :: array1(nclasses,nvars)
    real(kflt), intent(in)  :: array2(nclasses,nclasses,nvars*(nvars-1)/2)
    integer,    intent(out) :: error_code
    integer :: iv,jv
    integer :: is,js
    integer :: k

    error_code = 0

    write(unt,*) nvars,nclasses

    do iv = 1,nvars
       do is = 1,nclasses
          write(unt,'(2i4,f10.5)') iv,is,array1(is,iv)
       end do
    end do

    k = 0
    do jv = 1,nvars-1
       do iv = jv+1,nvars
          k = k + 1
          do is = 1,nclasses
             do js = 1,nclasses
                write(unt,'(4i4,f10.5)') iv,jv,is,js,array2(is,js,k)
             end do
          end do
       end do
    end do

  end subroutine dump_array_unit

  subroutine dump_array_file(filename,nvars,nclasses,array1,array2,&
       error_code)
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: nvars, nclasses
    real(kflt),       intent(in)  :: array1(nclasses,nvars)
    real(kflt),       intent(in)  :: array2(nclasses,nclasses,nvars*(nvars-1)/2)
    integer,          intent(out) :: error_code
    integer :: unt,err

    error_code = 0

    call units_open(filename,'unknown',unt,err)

    if( err /= 0 ) then
       write(0,*) 'ERROR ! error opening file '//trim(filename)
       error_code = 1
       return
    end if

    call dump_array_unit(unt,nvars,nclasses,array1,array2,err)

    close(unt)

  end subroutine dump_array_file

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
       if (nproc == 1 .or. np == nproc) then 
          ! single process or n. of seqs in checkpoint matches nproc
          ! => each proc reads a different seq
          do p = 1,np
             read(unt) dummy
             if (p == iproc+1) seq = dummy
          end do
       else
          ! if nproc is larger than 1 and different from the num. of seqs in the checkpoint 
          ! => randomize
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
