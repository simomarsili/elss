! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module fasta
  use constants
  implicit none
  private 
  public :: fasta_read
  public :: set_fasta_alphabet
  public :: fasta_alphabet
  integer, parameter :: kprot=21, knuc=6
  character(len=kprot) :: protein_alphabet='ACDEFGHIKLMNPQRSTVWY-'
  character(len=knuc)  :: nuc_acid_alphabet='ACGTU-'
  character(len=string_size) :: fasta_alphabet

contains

  subroutine set_fasta_alphabet(data_type)
    character(len=*), intent(in) :: data_type
    select case(data_type)
    case('protein')
       fasta_alphabet = protein_alphabet
    case('nuc_acid')
       fasta_alphabet = nuc_acid_alphabet
    case default
    end select
  end subroutine set_fasta_alphabet

  function is_nuc_acid(string) result(res)
    ! check that a string contain only natural nucleotides symbols
    character(len=*), intent(in)  :: string
    logical :: res
    
    res = .false.
    if (verify(trim(string),trim(nuc_acid_alphabet)) == 0) res = .true.
    
  end function is_nuc_acid
  
  function is_protein(string) result(res)
    ! check that a string contain only natural amino acids symbols
    character(len=*), intent(in)  :: string
    logical :: res
    
    res = .false.
    if (verify(trim(string),trim(protein_alphabet)) == 0) res = .true.
    
  end function is_protein

  subroutine fasta_string2seq(string,seq,err)
    ! read a sequence, return an array of integer
    character(len=*), intent(in)  :: string
    integer,          intent(out) :: seq(:)
    integer,          intent(out) :: err
    integer :: i

    seq = [(index(trim(fasta_alphabet),string(i:i)), i=1,len_trim(string))]

    if (any(seq == 0)) err = 1
    
  end subroutine fasta_string2seq

  subroutine fasta_read(unt,seqs,data_type,error_code,error_string)
    ! read n sequences from unit udata
    integer,              intent(in)  :: unt
    integer, allocatable, intent(out) :: seqs(:,:)
    character(len=*),     intent(out) :: data_type
    integer,              intent(out) :: error_code
    character(len=*),     intent(out) :: error_string
    character(len=long_string_size) :: header
    character(len=long_string_size) :: string
    character(len=long_string_size) :: line
    integer, allocatable :: seq(:)
    integer                         :: k,ns,nv,nnuc,nprot
    integer                         :: err

    ! count seqs and Ls
    ns = 0
    nnuc = 0
    nprot = 0
    do
       read(unt,'(a)',iostat=err) line
       ! exit at end of file or empty line
       if ( err < 0 .or. len_trim(line) == 0 ) exit
       if ( line(1:1) == ">" ) then
          ns = ns + 1
       else
          if (is_nuc_acid(line)) nnuc = nnuc + 1
          if (is_protein(line)) nprot = nprot + 1
       end if
    end do
    rewind(unt)

    write(*,*) ns,nprot,nnuc

    if (nnuc > nprot) then
       data_type = 'nuc_acid'
    else
       data_type = 'protein'
    end if

    call set_fasta_alphabet(data_type)

    string = ""
    k = 1
    read(unt,'(a)',iostat=err) header
    if (header(1:1) /= '>') then 
       error_code = 43
    end if
    do
       read(unt,'(a)',iostat=err) line
       if (err < 0) then 
          if (k == 1) then 
             nv = len_trim(string)
             allocate(seq(nv),stat=err)
             allocate(seqs(nv,ns),stat=err)
          end if
          call fasta_string2seq(string,seq,error_code)
          if (error_code > 0) return
          seqs(:,k) = seq
          exit 
       end if
       if (line(1:1) == ">") then 
          header = line
          if (k == 1) then 
             nv = len_trim(string)
             allocate(seq(nv),stat=err)
             allocate(seqs(nv,ns),stat=err)
          end if
          call fasta_string2seq(string,seq,error_code)
          if (error_code > 0) return
          seqs(:,k) = seq
          string = ""
          k = k + 1
          cycle
       else
          string = trim(string)//trim(line)
       end if
    end do

!    do ns = 1,nseqs
!       write(*,'(1000i2)') seqs(:,k)
!    end do

  end subroutine fasta_read

end module fasta
