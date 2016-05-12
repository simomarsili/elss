! Copyright (c) 2016 Simone Marsili
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in
! the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
! the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
! A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
! ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

module errors
  use kinds
  use constants
  implicit none
  private 
  public :: dump_error
  character(len=1), parameter                  :: nl=achar(10)
  character(len=long_string_size)              :: syntax = nl//& 
       '                                       elss (0.2.1)                                            '//nl//&
       '                                    ==================                                         '//nl//&
       '                                                                                               '//nl//&
       ' "elss" is a Monte Carlo (MC) code for the analysis and the inference of energy landscapes in  '//nl//&
       ' protein sequence spaces. elss can either be used to simulate a trajectory with a user-defined '//nl//&
       ' energy function, or to infer a data-driven statistical model from a multiple sequence         '//nl//&
       ' alignment (MSA) via maximum a posteriori (MAP) estimation. The energy function and its        '//nl//&
       ' parameters control both the frequencies of amino acids at the different positions along the   '//nl//& 
       ' chain and their correlations.                                                                 '//nl//& 
       nl//&
       nl//&
       'Option                         Description                                     (Default Value)  '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       ' (-h|--help)                   print this help message                         (None)         '//nl//&
       nl//&
       ' (-p|--prm) <path_to_file>     parameters file                                 (None)           '//nl//&   
       '        OR'//nl//&   
       ' (-r|--rst) <path_to_file>     restart file                                    (None)           '//nl//&   
       nl//&
       ' (-n|--nsweeps) <int>          num. of MC sweeps                               (0)              '//nl//&
       nl//&
       ' (-u|--nupdate) <int>          stride (as num. of sweeps) for averages updates (10)             '//nl//&
       nl//&
       ' --fasta <path_to_file>        data file (MSA format)                          (None)           '//nl//&
       '        OR'//nl//&   
       ' --raw <path_to_file>          data file ("raw" format)                        (None)           '//nl//&
!       '        OR'//nl//&   
!       ' --table <path_to_file>        data file ("table" format)                      (None)           '//nl//&
       nl//&
       ' (-w|--weights) <path_to_file> weights file                                    (None)           '//nl//&
       nl//&
       ' (-s|--seq) <path_to_file>     starting sequence file (SIM)                    (None)           '//nl//&
       nl//&
       ' --wid <float>                 %id threshold for weights calculation           (-1)             '//nl//&
       nl//&
       ' --learn-gd <int>              num. of gradient descent steps                  (0)              '//nl//&
       nl//&
       ' (--learn|--learn-agd) <int>   num. of accelerated gradient descent steps      (0)              '//nl//&
       nl//&
       ' (-l|--lambda) <float>         (scaled) regularization parameter               (0.01)           '//nl//&
       nl//&
       ' --random_seed <int>           initialize the random seed                      (0)              '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       nl//&
       nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       ' For more information and usage examples, please check the project github repository:           '//nl//&
       ' https://github.com/simomarsili/elss                                                            '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       '                                                                                                    '
  character(len=string_size), dimension(-1:45) :: err_msg = & 
       [ &
       ! -1 print syntax only
       '                                                                                         ', & 
       !  0 no errors
       '                                                                                         ', & 
       !  1 
       'elss: error                                                                              ', & 
       !  2
       'elss: invalid option                                                                     ', & 
       !  3
       'elss: EVAL/MAP mode: please specify a data_file ([-i|--raw|--table|--fasta] <data_file>) ', & 
       !  4
       'elss: the [-i|--raw|--table|--fasta] option needs an argument <data_file>)               ', & 
       !  5
       'elss: the [-w|--weights] option needs an argument <weights_file>)                        ', & 
       !  6
       'elss: the [-r|--rst] option needs an argument <rst_file>)                                ', & 
       !  7
       'elss: cannot access <data_file>: No such file or directory                               ', & 
       !  8
       'elss: cannot access <weights_file>: No such file or directory                            ', & 
       !  9 
       'elss: cannot access <rst_file>: No such file or directory                                ', & 
       !  10
       'elss: please check nsweeps value ([-n|--nsweeps] <nsweeps>)                              ', & 
       !  11
       'elss: please check niter value ([--agd] <niter>)                                         ', & 
       !  12
       'elss: please check niter value ([--gd] <niter>)                                          ', & 
       !  13
       'elss: please check wid value (--wid <wid>)                                               ', & 
       !  14 
       'elss: you must use either the [-w|--weights] or the --wid option                         ', & 
       !  15 
       'elss: please check <nupdate> value ([-u|--nupdate] <nupdate>)                            ', & 
       !  16 
       'elss: you must specify a single data_file ([-i|--raw|--table|--fasta] <data_file>)       ', & 
       !  17 
       'elss: you must specify a single weights_file ([-w|--weights] <weights_file>)             ', & 
       !  18 
       'elss: you must specify a single rst_file ([-r|--rst] <rst_file>)                         ', & 
       !  19 
       'elss: error opening file                                                                 ', & 
       !  20 
       'elss: error reading data_file                                                            ', & 
       !  21 
       'elss: error reading weights_file                                                         ', & 
       !  22 
       'elss: error reading rst_file                                                             ', & 
       !  23 
       'elss: data_file and weights file must have the same number of lines                      ', & 
       !  24 
       'elss: unkown amino acid                                                                  ',  & 
       !  25 
       'elss: in rst file, nvars is not consistent with data                                     ',  & 
       !  26 
       'elss: you must specify a single prm_file ([-p|--prm] <prm_file>)                         ', & 
       !  27
       'elss: the [-p|--prm] option needs an argument <prm_file>)                                ', & 
       !  28 
       'elss: cannot access <prm_file>: No such file or directory                                ', & 
       !  29 
       'elss: error reading fields in <prm_file>                                                 ', & 
       !  30 
       'elss: error reading couplings in <prm_file>                                              ', & 
       !  31 
       'elss: unexpected number of fields in <prm_file>                                          ', & 
       !  32 
       'elss: error opening array file                                                           ', & 
       !  33 
       'elss: unexpcted memory allocation                                                        ', & 
       !  34 
       'elss: error opening rst file                                                             ', & 
       !  35 
       'elss: rst and prm reads cannot go together                                               ', & 
       !  36 
       'elss: SIM/EVAL modes require a rst file or a prm file                                    ', & 
       !  37 
       'elss: MAP mode: please specify a nsweeps value ([-n|--nsweeps] <nsweeps>)                ', & 
       !  38 
       'elss: MAP mode: lambda hyperparameter cannot be negative ([-l|--lambda] <lamnda>)        ', & 
       !  39 
       'elss: in prm file: in protein, num. of classes must be 21 ([-p|--prm] <prm_file>)        ', & 
       !  40 
       'elss: you must specify a single seq_file ([-s|--seq] <seq_file>)                         ', & 
       !  41
       'elss: the [-s|--seq] option needs an argument <seq_file>)                                ', & 
       !  42 
       'elss: cannot access <seq_file>: No such file or directory                                ', & 
       !  43 
       'elss: not a valid fasta file                                                             ', & 
       !  44
       'elss: please check random seed value (--random_seed <rseed>)                             ', & 
       !  45
       'elss: check temperature factor [-t|--temp] <temp>)                                       ' & 
       ]

contains
  
  subroutine dump_error(error_code,string)
    integer, intent(in)          :: error_code
    character(len=*), intent(in) :: string

    if (error_code == 0) then 
       ! do nothing
       return
    else if (error_code == -1) then 
       write(0,'(a)') trim(syntax)
       return
    else if (error_code > 0) then 
       if (len_trim(string) > 0) then 
          write(0,'(a)') trim(err_msg(error_code))//": "//trim(string)
       else
          write(0,'(a)') trim(err_msg(error_code))
       end if       
    end if

  end subroutine dump_error

end module errors

