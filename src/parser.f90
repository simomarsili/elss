! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module parser
  implicit none
  private
  public :: parser_nfields, cleanline, remove_comments
  ! possible delimiters are: space, tab, comma, colon, semicolon
  !  character(len=1) :: delimiters(5)=[" ", achar(9), ",", ":", ";"]
  character(len=1) :: delimiters(2)=[" ", achar(9)] ! accepted delimiters are space and tab
  
contains

  subroutine parser_nfields(line,parsed,nfields)
    character(len=*),intent(in)  :: line
    integer,         intent(out) :: nfields
    character(len=*),intent(out) :: parsed
    integer i, n, toks
    
    i = 1
    n = len_trim(line)
    toks = 0
    nfields = 0
    parsed = ""

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

  subroutine cleanline(line,parsed)
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: parsed
    character(len=1) :: a
    
    integer :: i,n
    
    n = len_trim(line)
    parsed = line
    do while(i <= n)
       a = line(i:i)
       if ( any(delimiters == a) ) then
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
    integer :: i,n
    
    n = len_trim(line)
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
