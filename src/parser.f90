! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module parser
  use constants
  implicit none
  private
  public :: parser_nfields, cleanline, remove_comments, parse_next_line
  ! possible delimiters are: space, tab, comma, colon, semicolon
  !  character(len=1) :: delimiters(5)=[" ", achar(9), ",", ":", ";"]   
  ! accepted delimiters are space and tab
  character(len=1) :: delimiters(2) = [" ", achar(9)]
  
contains

  subroutine parse_next_line(unt, line, parsed_line, nfields, err)
    integer, intent(in) :: unt
    character(len=*), intent(out) :: line
    character(len=*), intent(out) :: parsed_line
    integer, intent(out) :: nfields
    integer, intent(out) :: err

    do
       read(unt,'(a)', iostat=err) line
       if (err < 0) exit
       call remove_comments(line)
       if (len_trim(line) == 0) cycle
       call parser_nfields(line, parsed_line, nfields)
       exit
    end do

  end subroutine parse_next_line
  
  subroutine parser_nfields(line, parsed, nfields)
    character(len=*), intent(in)  :: line
    integer,          intent(out) :: nfields
    character(len=*), intent(out) :: parsed
    integer :: i, n, toks
    
    toks = 0
    nfields = 0
    parsed = ""

    n = len_trim(line)
    i = 1
    outer: do while(i <= n)
       do while(any(delimiters == line(i:i)))
          i = i + 1
          if (n < i) exit outer
       enddo
       toks = toks + 1
       nfields = toks
       if(nfields == 1) then 
          parsed=trim(parsed)//line(i:i)
       else
          parsed=trim(parsed)//" "//line(i:i)
       end if
       do
          i = i + 1
          if (n < i) exit outer
          if (any(delimiters == line(i:i))) exit
          parsed=trim(parsed)//line(i:i)
       enddo
    enddo outer
    
  end subroutine parser_nfields

  subroutine cleanline(line, parsed)
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: parsed
    character(len=1) :: a
    integer :: i, n
    
    parsed = line
    n = len_trim(line)
    i = 1
    do while(i <= n)
       a = line(i:i)
       if (any(delimiters == a)) then
          parsed(i:i) = " "
       else
          parsed(i:i) = a
       end if
       i = i + 1
    end do

  end subroutine cleanline

  subroutine remove_comments(line)
    character(len=*), intent(inout)  :: line
    character(len=1) :: a
    integer :: i, n
    
    n = len_trim(line)
    i = 1
    do while(i <= n)
       a = line(i:i)
       if (a == "#") then 
          line(i:n)=''
          exit
       end if
       i = i + 1
    end do

  end subroutine remove_comments


end module parser
