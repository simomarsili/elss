! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module data
  use kinds
  use constants
  use parser, only: parser_nfields
  implicit none
  private

  public :: data_read
  public :: data_average

  real(kflt), allocatable                      :: ws(:)

contains

  subroutine data_read(iproc, udata, uwgt, wid, nvars, nclasses, data_type, &
                       ndata, neff, data, error_code)
    use fasta, only: fasta_read, fasta_alphabet
    integer,              intent(in)    :: iproc
    integer,              intent(in)    :: udata
    integer,              intent(in)    :: uwgt
    real(kflt),           intent(in)    :: wid
    integer,              intent(out)   :: nvars
    integer,              intent(out)   :: nclasses
    character(len=*),     intent(out)   :: data_type
    integer,              intent(out)   :: ndata
    real(kflt),           intent(out)   :: neff
    integer, allocatable, intent(out)   :: data(:,:)
    integer,              intent(out)   :: error_code
    integer                                      :: err
    character(len=long_string_size)              :: line, newline
    integer                                      :: nfields
    integer                                      :: i,nn,cmin

    error_code = 0
    select case(trim(data_type))
    case('int')
       
       ! get number of fields
       read(udata, '(a)', iostat=err) line
       call parser_nfields(line, newline, nfields)
       rewind(udata)
       
       ! set nvars
       nvars = nfields
       
       ! count data lines
       ndata = 0
       do
          read(udata, '(a)', iostat=err) line
          ! exit at end of file or empty line
          if( err < 0 .or. len_trim(line) == 0) exit
          ndata = ndata + 1
       end do
       rewind(udata)

       ! allocate memory for data
       allocate(data(nvars, ndata), stat=err)

       ! read data
       do i = 1, ndata
          read(udata, *, iostat=err) data(:, i)
          if(err > 0) then 
             error_code = 1
             return
          end if
       end do

       ! smallest index should be ZERO
       cmin = minval(data)
       if (cmin /= 0) then
          if (iproc == 0) write(0,*) 'cmin: ', cmin
          if (cmin < 0)then
             if (iproc == 0) write(0, *) "ERROR: class indices should start from zero."
             error_code = 3
             return
          else
             if (iproc == 0) write(0, *) "WARNING: class indices should start from zero."
          end if
       end if
       data = data + 1
       ! set n. of classes per variable as the max value in data
       nclasses = maxval(data)

    case('bio', 'protein', 'nuc_acid')

       ! read sequences from MSA
       call fasta_read(udata, data, data_type, error_code)
       if (error_code /= 0) return
       ndata = size(data, 2)
       nvars = size(data, 1)
       nclasses = len_trim(fasta_alphabet)

    end select

    allocate(ws(ndata), stat=err)
    ! initialize weigths to one
    ws = 1.0_kflt

    if (uwgt > 0) then 
       ! count lines
       nn = 0
       do
          read(uwgt, '(a)', iostat=err) line
          ! exit at end of file or empty line
          if( err < 0 .or. len_trim(line) == 0) exit
          nn = nn + 1
       end do

       rewind(uwgt)
       if (nn /= ndata) then
          error_code = 2
          return
       end if

       do i = 1,ndata
          read(uwgt, *, iostat=err) ws(i)
          if(err > 0) then 
             error_code = 2
             return
          end if
       end do

    end if

    if (wid > 0.0_kflt) then 

       call data_reweight(data, wid, iproc)

    end if

    neff = sum(ws)
  end subroutine data_read

  subroutine data_average(nvars, nclasses, ndata, neff, data, &
       data_freq_single, data_freq_pair)
    use averages, only: update_counts
    integer,    intent(in)    :: nvars
    integer,    intent(in)    :: nclasses
    integer,    intent(in)    :: ndata
    real(kflt), intent(inout) :: neff
    integer,    intent(in)    :: data(nvars, ndata)
    real(kflt), intent(out)   :: data_freq_single(nclasses, nvars)
    real(kflt), intent(out)   :: data_freq_pair(nclasses, nclasses, nvars*(nvars-1)/2)
    integer             :: err
    integer,allocatable :: seq(:)
    integer             :: k, iv, jv
    logical             :: pseudocount = .false.

    allocate(seq(nvars), stat=err)

    ! take averages
    data_freq_single = 0.0_kflt
    data_freq_pair = 0.0_kflt
    do k = 1, ndata
       call update_counts(ws(k), data(:,k), data_freq_single, &
                                 data_freq_pair)
    end do

    if (pseudocount) then
       ! regularize (+1 sequence)
       neff = neff + 1.0_kflt
       do iv = 1,nvars
          data_freq_single(:, iv) = data_freq_single(:, iv) + 1.0_kflt / real(nclasses)
       end do
       k = 0
       do jv = 1, nvars - 1
          do iv = jv + 1, nvars
             k = k + 1
             data_freq_pair(:, :, k) = data_freq_pair(:, :, k) + 1.0_kflt / real(nclasses**2)
          end do
       end do
    end if

    data_freq_single = data_freq_single / neff
    data_freq_pair = data_freq_pair / neff

    deallocate(seq)

  end subroutine data_average

  subroutine data_reweight(data, wid, iproc)
    use units, only: units_open
    integer,    intent(in) :: data(:,:)
    real(kflt), intent(in) :: wid
    integer,    intent(in) :: iproc
    integer              :: err
    integer              :: nv, ns
    integer              :: id, jd
    integer, allocatable :: x(:), y(:)
    integer              :: thr, uwgt
    logical              :: dump_weights = .false.

    nv = size(data, 1)
    ns = size(data, 2)
    
    allocate(x(nv), y(nv), stat=err)
    
    thr = nint(nv * wid * 0.01_kflt)
    ws = 1.0_kflt
    do id = 1, ns - 1
       if (iproc == 0 .and. mod(ns, 10000)==0) then
          write(0,'(a,f8.1)') &
               'computing weigths: ', 100.0 * real(id) / real(ns)
       end if
       x = data(:, id)
       do jd = id + 1, ns
          y = data(:, jd)
          if(count(x == y) >= thr) then
             ws(id) = ws(id) + 1.0_kflt
             ws(jd) = ws(jd) + 1.0_kflt
          end if
       end do
    end do
    
    ws = 1.0_kflt / ws

    if (dump_weights) then
       if (iproc == 0) then 
          call units_open('wgt', 'unknown', uwgt, err)
          do id = 1, ns
             write(uwgt, *) ws(id)
          end do
          close(uwgt)
       end if
    end if
    
    deallocate(x, y)

  end subroutine data_reweight

end module data
