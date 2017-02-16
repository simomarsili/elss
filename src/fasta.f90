! Copyright (C) 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module fasta
  use constants
  implicit none
  private 
  public :: protein_alphabet
  public :: fasta_read
  character(len=1) :: protein_alphabet(21) = &
       ['A','C','D','E','F','G','H','I','K',&
       'L','M','N','P','Q','R','S','T','V','W','Y','-']
  character(len=1) :: nuc_acid_alphabet(6) = &
       ['A','C','G','T','U','-']
  character(len=21) :: protein_set='ACDEFGHIKLMNPQRSTVWY-'
  character(len=5)  :: nuc_acid_set='ACGTU-'
  integer          :: aamap(239) ! 239 is n. of ascii codes (0 is null)

contains

  subroutine fasta_init(data_type)
    character(len=*), intent(in) :: data_type
    ! set the aa -> class mapping
    integer :: k,nn
    
    aamap = 0
    select case(data_type)
    case('protein')
       do k = 1,21
          nn = iachar(protein_alphabet(k))
          aamap(nn) = k
       end do
    case('nuc_acid')
       do k = 1,6
          nn = iachar(nuc_acid_alphabet(k))
          aamap(nn) = k
       end do
    end select
    
  end subroutine fasta_init

  function is_nuc_acid(string) result(res)
    character(len=*), intent(in)  :: string
    logical :: res

    res = .false.
    if (verify(string,nuc_acid_set) == 0) res = .true.
    
  end function is_nuc_acid
  
  function is_protein(string) result(res)
    character(len=*), intent(in)  :: string
    logical :: res

    res = .false.
    if (verify(string,protein_set) == 0) res = .true.
    
  end function is_protein

  function aa_to_class(aa) result(c)
    ! takes an a.a. as input, return its index
    character(len=1), intent(in)  :: aa
    integer :: c
    
    c = aamap(iachar(aa))

  end function aa_to_class

  subroutine fasta_string2seq(string,seq,error_code,error_string)
    ! read a sequence, return an array of integer
    character(len=*), intent(in)  :: string
    integer,          intent(out) :: seq(:)
    integer,          intent(out) :: error_code
    character(len=*), intent(out) :: error_string
    integer :: i,n,index
    integer :: err

    error_string = ''
    n = len_trim(string)
    do i = 1,n
       index = aa_to_class(string(i:i))
       if ( index == 0 ) then ! symbol is not in protein_alphabet
          error_code = 24
          error_string = string(i:i)
          return
       else
          seq(i) = index
       end if
    end do
    
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

    data_type = 'protein'
    if (nnuc > nprot) data_type = 'nuc_acid'

    call fasta_init(data_type)

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
          call fasta_string2seq(string,seq,error_code,error_string)
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
          call fasta_string2seq(string,seq,error_code,error_string)
          if (error_code > 0) return
          seqs(:,k) = seq
          string = ""
          k = k + 1
          cycle
       else
          string = trim(string)//trim(line)
       end if
    end do

  end subroutine fasta_read

end module fasta
