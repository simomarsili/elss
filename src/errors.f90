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
       '                                       mcDCA (0.2.1)                                           '//nl//&
       '                                    ==================                                         '//nl//&
       '                                                                                               '//nl//&
       ' mcDCA  is a Monte Carlo (MC) code for the analysis and the inference of energy landscapes in  '//nl//&
       ' protein sequence spaces. mcDAC can either be used to simulate a trajectory with a user-defined'//nl//&
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
       ' https://github.com/simomarsili/mcDCA                                                           '//nl//&
       '------------------------------------------------------------------------------------------------'//nl//&
       '                                                                                                    '
  character(len=string_size), dimension(-1:46) :: err_msg = & 
       [ &
       ! -1 print syntax only
       '                                                                                          ', & 
       !  0 no errors
       '                                                                                          ', & 
       !  1 
       'mcDCA: error                                                                              ', & 
       !  2
       'mcDCA: invalid option                                                                     ', & 
       !  3
       'mcDCA: EVAL/MAP mode: please specify a data_file ([-i|--raw|--table|--fasta] <data_file>) ', & 
       !  4
       'mcDCA: the [-i|--raw|--table|--fasta] option needs an argument <data_file>)               ', & 
       !  5
       'mcDCA: the [-w|--weights] option needs an argument <weights_file>)                        ', & 
       !  6
       'mcDCA: the [-r|--rst] option needs an argument <rst_file>)                                ', & 
       !  7
       'mcDCA: cannot access <data_file>: No such file or directory                               ', & 
       !  8
       'mcDCA: cannot access <weights_file>: No such file or directory                            ', & 
       !  9 
       'mcDCA: cannot access <rst_file>: No such file or directory                                ', & 
       !  10
       'mcDCA: please check nsweeps value ([-n|--nsweeps] <nsweeps>)                              ', & 
       !  11
       'mcDCA: please check niter value ([--agd] <niter>)                                         ', & 
       !  12
       'mcDCA: please check niter value ([--gd] <niter>)                                          ', & 
       !  13
       'mcDCA: please check wid value (--wid <wid>)                                               ', & 
       !  14 
       'mcDCA: you must use either the [-w|--weights] or the --wid option                         ', & 
       !  15 
       'mcDCA: please check <nupdate> value ([-u|--nupdate] <nupdate>)                            ', & 
       !  16 
       'mcDCA: you must specify a single data_file ([-i|--raw|--table|--fasta] <data_file>)       ', & 
       !  17 
       'mcDCA: you must specify a single weights_file ([-w|--weights] <weights_file>)             ', & 
       !  18 
       'mcDCA: you must specify a single rst_file ([-r|--rst] <rst_file>)                         ', & 
       !  19 
       'mcDCA: error opening file                                                                 ', & 
       !  20 
       'mcDCA: error reading data_file                                                            ', & 
       !  21 
       'mcDCA: error reading weights_file                                                         ', & 
       !  22 
       'mcDCA: error reading rst_file                                                             ', & 
       !  23 
       'mcDCA: data_file and weights file must have the same number of lines                      ', & 
       !  24 
       'mcDCA: unkown amino acid                                                                  ',  & 
       !  25 
       'mcDCA: in rst file, nvars is not consistent with data                                     ',  & 
       !  26 
       'mcDCA: you must specify a single prm_file ([-p|--prm] <prm_file>)                         ', & 
       !  27
       'mcDCA: the [-p|--prm] option needs an argument <prm_file>)                                ', & 
       !  28 
       'mcDCA: cannot access <prm_file>: No such file or directory                                ', & 
       !  29 
       'mcDCA: error reading fields in <prm_file>                                                 ', & 
       !  30 
       'mcDCA: error reading couplings in <prm_file>                                              ', & 
       !  31 
       'mcDCA: unexpected number of fields in <prm_file>                                          ', & 
       !  32 
       'mcDCA: error opening array file                                                           ', & 
       !  33 
       'mcDCA: unexpcted memory allocation                                                        ', & 
       !  34 
       'mcDCA: error opening rst file                                                             ', & 
       !  35 
       'mcDCA: rst and prm reads cannot go together                                               ', & 
       !  36 
       'mcDCA: SIM/EVAL modes require a rst file or a prm file                                    ', & 
       !  37 
       'mcDCA: MAP mode: please specify a nsweeps value ([-n|--nsweeps] <nsweeps>)                ', & 
       !  38 
       'mcDCA: MAP mode: lambda hyperparameter cannot be negative ([-l|--lambda] <lamnda>)        ', & 
       !  39 
       'mcDCA: in prm file: in protein, num. of classes must be 21 ([-p|--prm] <prm_file>)        ', & 
       !  40 
       'mcDCA: you must specify a single seq_file ([-s|--seq] <seq_file>)                         ', & 
       !  41
       'mcDCA: the [-s|--seq] option needs an argument <seq_file>)                                ', & 
       !  42 
       'mcDCA: cannot access <seq_file>: No such file or directory                                ', & 
       !  43 
       'mcDCA: not a valid fasta file                                                             ', & 
       !  44
       'mcDCA: please check random seed value (--random_seed <rseed>)                             ', & 
       !  45
       'mcDCA: check temperature factor [-t|--temp] <temp>                                        ', &
       !  46
       'mcdca-sample: needs one between rst file (-r <rst>) or prm file (-p <prm>) as argument    ' & 
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

