What do quality scores mean?
============================

Phred Quality Scores
--------------------

We can calculate the *phred quality score* from the probability of
sequencing error (i.e. the base call is wrong) using:

Alternatively, we can rearrange to calculate the probability of error
from the *phred quality score* using:

Where *Q* is the *phred quality score* and *p* is the probability of
error (i.e. the base call is wrong)

``` python
def error_prob(quality):
    qval = quality
    return 10**(qval/-10.0)

print ("{0:^5}  {1:^10}".format("Phred", "Prob of"))
print ("{0:^5}  {1:^10}".format("score", "Error"))
for phred in range(0,42):
    print ("{0:^5}  {1:03.5f}".format(phred, error_prob(phred)))
```

    ## Phred   Prob of  
    ## score    Error   
    ##   0    1.00000
    ##   1    0.79433
    ##   2    0.63096
    ##   3    0.50119
    ##   4    0.39811
    ##   5    0.31623
    ##   6    0.25119
    ##   7    0.19953
    ##   8    0.15849
    ##   9    0.12589
    ##  10    0.10000
    ##  11    0.07943
    ##  12    0.06310
    ##  13    0.05012
    ##  14    0.03981
    ##  15    0.03162
    ##  16    0.02512
    ##  17    0.01995
    ##  18    0.01585
    ##  19    0.01259
    ##  20    0.01000
    ##  21    0.00794
    ##  22    0.00631
    ##  23    0.00501
    ##  24    0.00398
    ##  25    0.00316
    ##  26    0.00251
    ##  27    0.00200
    ##  28    0.00158
    ##  29    0.00126
    ##  30    0.00100
    ##  31    0.00079
    ##  32    0.00063
    ##  33    0.00050
    ##  34    0.00040
    ##  35    0.00032
    ##  36    0.00025
    ##  37    0.00020
    ##  38    0.00016
    ##  39    0.00013
    ##  40    0.00010
    ##  41    0.00008

Ascii Codes
-----------

In FASTQ files, phred scores are represented using characters. Each
character on the keyboard can be represented by a number, called an
ascii code.

``` python
print ("{0:^8}  {1:^8}".format("Character", "ASCII #"))
for i in range(33,90):
    print("{0:^8}  {1:^8}".format(chr(i),i))
```

    ## Character  ASCII # 
    ##    !         33   
    ##    "         34   
    ##    #         35   
    ##    $         36   
    ##    %         37   
    ##    &         38   
    ##    '         39   
    ##    (         40   
    ##    )         41   
    ##    *         42   
    ##    +         43   
    ##    ,         44   
    ##    -         45   
    ##    .         46   
    ##    /         47   
    ##    0         48   
    ##    1         49   
    ##    2         50   
    ##    3         51   
    ##    4         52   
    ##    5         53   
    ##    6         54   
    ##    7         55   
    ##    8         56   
    ##    9         57   
    ##    :         58   
    ##    ;         59   
    ##    <         60   
    ##    =         61   
    ##    >         62   
    ##    ?         63   
    ##    @         64   
    ##    A         65   
    ##    B         66   
    ##    C         67   
    ##    D         68   
    ##    E         69   
    ##    F         70   
    ##    G         71   
    ##    H         72   
    ##    I         73   
    ##    J         74   
    ##    K         75   
    ##    L         76   
    ##    M         77   
    ##    N         78   
    ##    O         79   
    ##    P         80   
    ##    Q         81   
    ##    R         82   
    ##    S         83   
    ##    T         84   
    ##    U         85   
    ##    V         86   
    ##    W         87   
    ##    X         88   
    ##    Y         89

Phred Encodings
---------------

There are several different ways to encode phred scores with ascii
characters. The two most common are called phred+33 and phred+64. The
names are strange until you understand how then encoding works.

### Phred+33

To use the phred+33 encoding, take the phred quality score, add 33 to
it, then use the ascii character corresponding to the sum. For example,
using the phred+33 encoding, a quality score of 30 would be represented
with the ascii character with the ascii code of 63 (30 + 33), which is
‘?’.

### Phred+64

The phred+64 encoding works the same as the phred+33 encoding, except
you add 64 to the phred score to determine the ascii code of the quality
character. You will only find phred+64 encoding on older data, which was
sequenced several years ago. The tricky part is that there is no
indication in the FASTQ file as to which encoding was used, you have to
make an educated guess.

``` python
def error_prob(quality):
    qval = quality
    return 10**(qval/-10.0)
    
print ("{0:^5}  {1:^8}  {2:^8}  {3:^8}".format("Phred",  "Prob of", "Phred+33", "Phred+64"))
print ("{0:^5}  {1:^8}  {2:^8}  {3:^8}".format("score",  "Error", "Ascii", "Ascii"))
for phred in range(0,42):
    print ("{0:^5}  {1:03.5f}  {2:^8}  {3:^8}".format(phred, error_prob(phred), chr(phred+33), chr(phred+64)))
```

    ## Phred  Prob of   Phred+33  Phred+64
    ## score   Error     Ascii     Ascii  
    ##   0    1.00000     !         @    
    ##   1    0.79433     "         A    
    ##   2    0.63096     #         B    
    ##   3    0.50119     $         C    
    ##   4    0.39811     %         D    
    ##   5    0.31623     &         E    
    ##   6    0.25119     '         F    
    ##   7    0.19953     (         G    
    ##   8    0.15849     )         H    
    ##   9    0.12589     *         I    
    ##  10    0.10000     +         J    
    ##  11    0.07943     ,         K    
    ##  12    0.06310     -         L    
    ##  13    0.05012     .         M    
    ##  14    0.03981     /         N    
    ##  15    0.03162     0         O    
    ##  16    0.02512     1         P    
    ##  17    0.01995     2         Q    
    ##  18    0.01585     3         R    
    ##  19    0.01259     4         S    
    ##  20    0.01000     5         T    
    ##  21    0.00794     6         U    
    ##  22    0.00631     7         V    
    ##  23    0.00501     8         W    
    ##  24    0.00398     9         X    
    ##  25    0.00316     :         Y    
    ##  26    0.00251     ;         Z    
    ##  27    0.00200     <         [    
    ##  28    0.00158     =         \    
    ##  29    0.00126     >         ]    
    ##  30    0.00100     ?         ^    
    ##  31    0.00079     @         _    
    ##  32    0.00063     A         `    
    ##  33    0.00050     B         a    
    ##  34    0.00040     C         b    
    ##  35    0.00032     D         c    
    ##  36    0.00025     E         d    
    ##  37    0.00020     F         e    
    ##  38    0.00016     G         f    
    ##  39    0.00013     H         g    
    ##  40    0.00010     I         h    
    ##  41    0.00008     J         i

Why +33?
--------

ASCII 33 is the first “normal” ASCII character that. 1 through 32
include whitespace and non-printing characters, which cannot be
identified by eye)

``` python
print ("{0:^8}  {1:^8}".format("Character", "ASCII #"))
for i in range(0,40):
    print("{0:^8}  {1:^8}".format(chr(i),i))
```

    ## Character  ASCII # 
    ##    
    ##             1    
    ##             2    
    ##             3    
    ##             4    
    ##             5    
    ##             6    
    ##             7    
    ##             8    
    ##               9    
    ##    
    ##          10   
    ##             11   
    ##             12   
    ##             13   
    ##             14   
    ##             15   
    ##             16   
    ##             17   
    ##             18   
    ##             19   
    ##             20   
    ##             21   
    ##             22   
    ##             23   
    ##             24   
    ##             25   
    ##             26   
    ##             27   
    ##             28   
    ##             29   
    ##             30   
    ##             31   
    ##              32   
    ##    !         33   
    ##    "         34   
    ##    #         35   
    ##    $         36   
    ##    %         37   
    ##    &         38   
    ##    '         39
