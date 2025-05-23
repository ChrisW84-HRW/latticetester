
========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 1e-06.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  943      3324     11209     19955     33615 
LLL99999              820      4475     33616     91597    160583 
BKZ99999-10           710      4350     38384    135419    337773 
L5+9+BKZ-10           901      3826     24731    100375    280424 
LLL99999+BB           908      5062     39989    230422  46109593 
BKZ99999-6+BB         789      5057     43433    204356  14545500 
BKZ99999-8+BB         709      4887     44246    211109  12502489 
BKZ99999-10+BB        679      4752     45000    212070  10807650 
BKZ99999-12+BB        667      4622     45504    222949  10009146 
BKZ999-6+BB           672      4589     43586    208297  13612553 
BKZ999-8+BB           666      4558     44708    212735  14055653 
BKZ999-10+BB          657      4517     45131    218260  11661822 
BKZ999-12+BB          655      4496     45695    218987  11069404 
L8+BKZ-10+BB          950      5014     35341    181441  13332749 
L5+9+BKZ-10+BB        943      4676     31575    179011  12970786 
L5+9+BKZ-12+BB        781      4446     31831    192277  10618784 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5           1.040622168e+19 4.396980005e+21 2.306416649e+23 9.0331178e+23 1.208432454e+24 
LLL99999       1.040622168e+19 4.205880845e+21 1.146463186e+23 4.33286451e+23 9.83282428e+23 
BKZ99999-10    1.040622168e+19 4.205880845e+21 1.144617381e+23 4.167450442e+23 8.663775181e+23 
L5+9+BKZ-10    1.040622168e+19 4.205880845e+21 1.143017526e+23 4.129312932e+23 8.747641965e+23 
LLL99999+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-6+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-8+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-6+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-8+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-12+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L8+BKZ-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        654      94287   35676811 
BKZ99999-6+BB           5         16        566      51105   11069270 
BKZ99999-8+BB           5         16        553      52183    9426568 
BKZ99999-10+BB          5         16        573      45066    8071678 
BKZ99999-12+BB          5         16        575      44060    7434917 
BKZ999-6+BB             5         16        553      54488   10298185 
BKZ999-8+BB             5         16        559      53281   10650118 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        576      48460   10116347 
L5+9+BKZ-10+BB          5         16        566      46857    9832036 
L5+9+BKZ-12+BB          5         16        559      44873    7941386 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          4 
BKZ99999-6+BB           1          1          1          1          2 
BKZ99999-8+BB           1          1          1          1          2 
BKZ99999-10+BB          1          1          1          1          2 
BKZ99999-12+BB          1          1          1          1          2 
BKZ999-6+BB             1          1          1          1          2 
BKZ999-8+BB             1          1          1          1          2 
BKZ999-10+BB            1          1          1          1          2 
BKZ999-12+BB            1          1          1          1          2 
L8+BKZ-10+BB            1          1          1          1          2 
L5+9+BKZ-10+BB          1          1          1          1          2 
L5+9+BKZ-12+BB          1          1          1          1          2 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          5          9         24 
BKZ99999-6+BB           1          1          4         10         27 
BKZ99999-8+BB           1          1          3          9         29 
BKZ99999-10+BB          1          1          4          8         24 
BKZ99999-12+BB          1          1          4         11         34 
BKZ999-6+BB             1          1          4         11         30 
BKZ999-8+BB             1          1          4          9         36 
BKZ999-10+BB            1          1          4         15         21 
BKZ999-12+BB            1          1          4          9         25 
L8+BKZ-10+BB            1          1          4          9         30 
L5+9+BKZ-10+BB          1          1          4         10         22 
L5+9+BKZ-12+BB          1          1          4          9         27 

Total time for everything: 185.654606 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 1e-06.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  652      2083      5425      7386     11048 
LLL99999              575      2504     13442     32544     46866 
BKZ99999-10           562      2493     17486     67449    167862 
L5+9+BKZ-10           660      2752     17644     65824    161940 
LLL99999+BB           645      3175     19411    211171  74706138 
BKZ99999-6+BB         600      3152     22303    144974  20187485 
BKZ99999-8+BB         539      3020     23071    150297  15539536 
BKZ99999-10+BB        517      2914     23765    155551  16512458 
BKZ99999-12+BB        500      2840     24346    152889  14400252 
BKZ999-6+BB           501      2779     22469    148945  22141230 
BKZ999-8+BB           493      2772     23128    153519  16677538 
BKZ999-10+BB          496      2735     23692    157140  17743893 
BKZ999-12+BB          493      2702     24443    150735  13385732 
L8+BKZ-10+BB          637      3308     22829    148009  14200000 
L5+9+BKZ-10+BB        721      3495     23625    157733  14030575 
L5+9+BKZ-12+BB        666      3359     24130    160991  12271449 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.16      14.28      11.94 
L5+9+BKZ-10      36395.12      231.4      25.14      14.26      11.86 
LLL99999+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ99999-6+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-8+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-10+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ99999-12+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ999-6+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-8+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-10+BB     36395.12      231.4      25.14       14.1       11.5 
BKZ999-12+BB     36395.12      231.4      25.14       14.1       11.5 
L8+BKZ-10+BB     36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-10+BB   36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-12+BB   36395.12      231.4      25.14       14.1       11.5 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        793     129459   58161623 
BKZ99999-6+BB           5         16        726      61994   15542899 
BKZ99999-8+BB           5         16        684      61143   11909034 
BKZ99999-10+BB          5         16        701      57957   12646862 
BKZ99999-12+BB          5         16        693      53608   11006050 
BKZ999-6+BB             5         16        716      66232   17066547 
BKZ999-8+BB             5         16        690      63917   12828353 
BKZ999-10+BB            5         16        697      58726   13659907 
BKZ999-12+BB            5         16        680      52185   10183027 
L8+BKZ-10+BB            5         16        687      56562   10830150 
L5+9+BKZ-10+BB          5         16        680      61158   10721944 
L5+9+BKZ-12+BB          5         16        668      56555    9309786 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          7 
BKZ99999-6+BB           1          1          1          2          6 
BKZ99999-8+BB           1          1          1          2          6 
BKZ99999-10+BB          1          1          1          1          5 
BKZ99999-12+BB          1          1          1          2          5 
BKZ999-6+BB             1          1          1          2          6 
BKZ999-8+BB             1          1          1          2          6 
BKZ999-10+BB            1          1          1          1          5 
BKZ999-12+BB            1          1          1          2          5 
L8+BKZ-10+BB            1          1          1          2          5 
L5+9+BKZ-10+BB          1          1          1          2          5 
L5+9+BKZ-12+BB          1          1          1          2          5 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          5         14         25 
BKZ99999-6+BB           1          1          4         11         40 
BKZ99999-8+BB           1          1          5         11         21 
BKZ99999-10+BB          1          1          4         11         45 
BKZ99999-12+BB          1          1          4          8         24 
BKZ999-6+BB             1          1          3         14         24 
BKZ999-8+BB             1          1          4         11         21 
BKZ999-10+BB            1          1          4         12         31 
BKZ999-12+BB            1          1          3          9         33 
L8+BKZ-10+BB            1          1          4          8         25 
L5+9+BKZ-10+BB          1          1          4         14         23 
L5+9+BKZ-12+BB          1          1          4         12         20 

Total time for everything: 254.748349 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 2992     17085     85279    186800    340264 
LLL99999             3314     28847    284335    885948   1682416 
BKZ99999-10          3325     30293    342985   1382131   3813978 
L5+9+BKZ-10          3319     22123    213088   1001398   3065160 
LLL99999+BB          3516     30077    297582   1913732 422451622 
BKZ99999-6+BB        3470     31245    339095   1770468 134462642 
BKZ99999-8+BB        3412     31339    346798   1845630 114818916 
BKZ99999-10+BB       3384     31351    355284   1886831  99610505 
BKZ99999-12+BB       3382     31244    360932   2016390  92580690 
BKZ999-6+BB          3371     31014    337847   1803019 125382161 
BKZ999-8+BB          3348     31168    350272   1864974 129476075 
BKZ999-10+BB         3339     31201    356205   1950036 107497580 
BKZ999-12+BB         3335     31189    362009   1959490 102080223 
L8+BKZ-10+BB         3600     27295    268610   1613923 110236039 
L5+9+BKZ-10+BB       3514     23317    225507   1557721 135101433 
L5+9+BKZ-12+BB       3405     23208    230238   1630088  99854746 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5           1.040622168e+19 4.396980005e+21 2.306416649e+23 9.0331178e+23 1.208432454e+24 
LLL99999       1.040622168e+19 4.205880845e+21 1.146463186e+23 4.33286451e+23 9.83282428e+23 
BKZ99999-10    1.040622168e+19 4.205880845e+21 1.144617381e+23 4.167450442e+23 8.663775181e+23 
L5+9+BKZ-10    1.040622168e+19 4.205880845e+21 1.143917209e+23 4.166116336e+23 8.805377564e+23 
LLL99999+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-6+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-8+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-6+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-8+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-12+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L8+BKZ-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        654      94287   35676811 
BKZ99999-6+BB           5         16        566      51105   11069270 
BKZ99999-8+BB           5         16        553      52183    9426568 
BKZ99999-10+BB          5         16        573      45066    8071678 
BKZ99999-12+BB          5         16        575      44060    7434917 
BKZ999-6+BB             5         16        553      54488   10293946 
BKZ999-8+BB             5         16        559      53194   10650118 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        573      45842    8996138 
L5+9+BKZ-10+BB          5         16        562      49758   11185848 
L5+9+BKZ-12+BB          5         16        564      45664    8079564 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          4 
BKZ99999-6+BB           1          1          1          1          2 
BKZ99999-8+BB           1          1          1          1          2 
BKZ99999-10+BB          1          1          1          1          2 
BKZ99999-12+BB          1          1          1          1          2 
BKZ999-6+BB             1          1          1          1          2 
BKZ999-8+BB             1          1          1          1          2 
BKZ999-10+BB            1          1          1          1          2 
BKZ999-12+BB            1          1          1          1          2 
L8+BKZ-10+BB            1          1          1          1          2 
L5+9+BKZ-10+BB          1          1          1          1          2 
L5+9+BKZ-12+BB          1          1          1          1          2 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          5          9         24 
BKZ99999-6+BB           1          1          4         10         27 
BKZ99999-8+BB           1          1          3          9         29 
BKZ99999-10+BB          1          1          4          8         24 
BKZ99999-12+BB          1          1          4         11         34 
BKZ999-6+BB             1          1          4         11         30 
BKZ999-8+BB             1          1          4          9         36 
BKZ999-10+BB            1          1          4         15         21 
BKZ999-12+BB            1          1          4          9         25 
L8+BKZ-10+BB            1          1          4         10         17 
L5+9+BKZ-10+BB          1          1          4          9         23 
L5+9+BKZ-12+BB          1          1          4         11         21 

Total time for everything: 1713.00353 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 1613      6312     19454     29031     42729 
LLL99999             1660      8968     61677    165377    241469 
BKZ99999-10          1730      9959     93390    413467   1147713 
L5+9+BKZ-10          1959      9723     85421    381013   1061006 
LLL99999+BB          1843      9968     75143   1569132 687076890 
BKZ99999-6+BB        1860     10770     93799    975037 185340381 
BKZ99999-8+BB        1812     10797     99778   1018511 147462500 
BKZ99999-10+BB       1791     10793    105691   1057103 150663394 
BKZ99999-12+BB       1770     10758    112533   1052217 130050402 
BKZ999-6+BB          1762     10526     93556   1007696 202583060 
BKZ999-8+BB          1753     10623    100009   1054403 155664224 
BKZ999-10+BB         1755     10659    105500   1073148 158816667 
BKZ999-12+BB         1748     10629    112657   1035480 125855747 
L8+BKZ-10+BB         1977     10750     97724   1034137 137540115 
L5+9+BKZ-10+BB       2119     10751     97703   1045133 153613152 
L5+9+BKZ-12+BB       2073     10645    102190   1085147 113774566 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.16      14.28      11.96 
L5+9+BKZ-10      36395.12      231.4      25.14      14.34      11.92 
LLL99999+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ99999-6+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-8+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-10+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ99999-12+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ999-6+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-8+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-10+BB     36395.12      231.4      25.14       14.1       11.5 
BKZ999-12+BB     36395.12      231.4      25.14       14.1       11.5 
L8+BKZ-10+BB     36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-10+BB   36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-12+BB   36395.12      231.4      25.14       14.1       11.5 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        793     129460   58178894 
BKZ99999-6+BB           5         16        727      61892   15470688 
BKZ99999-8+BB           5         16        683      61332   12334879 
BKZ99999-10+BB          5         16        681      58005   12572744 
BKZ99999-12+BB          5         16        669      54628   10821560 
BKZ999-6+BB             5         16        716      65612   16916204 
BKZ999-8+BB             5         16        689      64748   12989651 
BKZ999-10+BB            5         16        686      58810   13281881 
BKZ999-12+BB            5         16        669      52984   10469421 
L8+BKZ-10+BB            5         16        671      60083   11459510 
L5+9+BKZ-10+BB          5         16        689      60103   12842552 
L5+9+BKZ-12+BB          5         16        674      58309    9463081 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          3          8 
BKZ99999-6+BB           1          1          1          2          6 
BKZ99999-8+BB           1          1          1          2          6 
BKZ99999-10+BB          1          1          1          2          5 
BKZ99999-12+BB          1          1          1          2          5 
BKZ999-6+BB             1          1          1          2          6 
BKZ999-8+BB             1          1          1          2          6 
BKZ999-10+BB            1          1          1          2          6 
BKZ999-12+BB            1          1          1          2          5 
L8+BKZ-10+BB            1          1          1          2          6 
L5+9+BKZ-10+BB          1          1          1          2          5 
L5+9+BKZ-12+BB          1          1          1          2          5 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          5         14         25 
BKZ99999-6+BB           1          1          4         11         40 
BKZ99999-8+BB           1          1          5         11         21 
BKZ99999-10+BB          1          1          4         11         45 
BKZ99999-12+BB          1          1          3          8         24 
BKZ999-6+BB             1          1          3         14         24 
BKZ999-8+BB             1          1          4         11         21 
BKZ999-10+BB            1          1          4         12         31 
BKZ999-12+BB            1          1          3         10         33 
L8+BKZ-10+BB            1          1          3         13         26 
L5+9+BKZ-10+BB          1          1          4         14         25 
L5+9+BKZ-12+BB          1          1          3         10         25 

Total time for everything: 2366.691793 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                13720     79835    370115    676384   1052621 
LLL99999            15315    136642   1340935   4086383   7264513 
BKZ99999-10         15357    142402   1600088   6385520  17244085 
L5+9+BKZ-10         14344     99230    960039   4641511  14217589 
LLL99999+BB         15889    139949   1393955  10953059 3098692970 
BKZ99999-6+BB       15829    145209   1583743   9253582 965852457 
BKZ99999-8+BB       15641    145650   1614814   9583479 818782991 
BKZ99999-10+BB      15517    145829   1647322   9623425 704652469 
BKZ99999-12+BB      15538    145910   1664453  10180777 653303600 
BKZ999-6+BB         15449    144935   1578939   9496454 894326231 
BKZ999-8+BB         15457    145397   1632458   9764053 920392824 
BKZ999-10+BB        15404    145686   1652667   9981140 760992925 
BKZ999-12+BB        15380    145579   1671099   9961188 724199532 
L8+BKZ-10+BB        15820    124143   1219088   8259982 879219312 
L5+9+BKZ-10+BB      15054    102665   1008410   8023320 856794930 
L5+9+BKZ-12+BB      14725    102579   1021777   8676823 694499592 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5           1.040622168e+19 4.396980005e+21 2.306416649e+23 9.0331178e+23 1.208432454e+24 
LLL99999       1.040622168e+19 4.205880845e+21 1.146463186e+23 4.33286451e+23 9.83282428e+23 
BKZ99999-10    1.040622168e+19 4.205880845e+21 1.144617381e+23 4.221671674e+23 8.791854587e+23 
L5+9+BKZ-10    1.040622168e+19 4.205880845e+21 1.143017526e+23 4.140416639e+23 8.879795872e+23 
LLL99999+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-6+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-8+BB  1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ99999-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-6+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-8+BB    1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
BKZ999-12+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L8+BKZ-10+BB   1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-10+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 
L5+9+BKZ-12+BB 1.040622168e+19 4.205880845e+21 1.141502843e+23 4.075344092e+23 8.290890941e+23 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        654      94287   35676673 
BKZ99999-6+BB           5         16        566      51105   11069270 
BKZ99999-8+BB           5         16        553      51967    9426568 
BKZ99999-10+BB          5         16        573      45135    8071678 
BKZ99999-12+BB          5         16        575      44060    7434917 
BKZ999-6+BB             5         16        553      54488   10298185 
BKZ999-8+BB             5         16        559      53351   10562930 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        576      48460   10116347 
L5+9+BKZ-10+BB          5         16        566      46857    9832892 
L5+9+BKZ-12+BB          5         16        559      44873    7941386 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          4 
BKZ99999-6+BB           1          1          1          1          2 
BKZ99999-8+BB           1          1          1          1          2 
BKZ99999-10+BB          1          1          1          1          2 
BKZ99999-12+BB          1          1          1          1          2 
BKZ999-6+BB             1          1          1          1          2 
BKZ999-8+BB             1          1          1          1          2 
BKZ999-10+BB            1          1          1          1          2 
BKZ999-12+BB            1          1          1          1          2 
L8+BKZ-10+BB            1          1          1          1          2 
L5+9+BKZ-10+BB          1          1          1          1          2 
L5+9+BKZ-12+BB          1          1          1          1          2 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          5          9         24 
BKZ99999-6+BB           1          1          4         10         27 
BKZ99999-8+BB           1          1          3          9         29 
BKZ99999-10+BB          1          1          4          8         24 
BKZ99999-12+BB          1          1          4         11         34 
BKZ999-6+BB             1          1          4         11         30 
BKZ999-8+BB             1          1          4          9         36 
BKZ999-10+BB            1          1          4         15         21 
BKZ999-12+BB            1          1          4          9         25 
L8+BKZ-10+BB            1          1          4          9         30 
L5+9+BKZ-10+BB          1          1          4         10         22 
L5+9+BKZ-12+BB          1          1          4          9         27 

Total time for everything: 12165.35814 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.
Safety margin on the BB bounds: epsBounds = 0.1e-5.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 8845     43126    152912    228015    335502 
LLL99999             9240     67498    563412   1576207   2330187 
BKZ99999-10          9186     72708    785528   3670170  10433339 
L5+9+BKZ-10          9427     63173    678849   3439052   9688762 
LLL99999+BB          9701     70740    625219  11027593 5100843952 
BKZ99999-6+BB        9636     75436    781220   7226666 1357873337 
BKZ99999-8+BB        9389     75706    815832   7796159 1019176165 
BKZ99999-10+BB       9304     75860    841585   7931776 1072554320 
BKZ99999-12+BB       9305     75993    866260   7720035 930906621 
BKZ999-6+BB          9293     75141    778348   7441295 1613411220 
BKZ999-8+BB          9286     75566    815661   7867235 1138625542 
BKZ999-10+BB         9244     75647    839688   7987847 1164112942 
BKZ999-12+BB         9184     75476    865278   7548353 862024962 
L8+BKZ-10+BB         9858     70647    734468   7175372 891125598 
L5+9+BKZ-10+BB      10135     66532    732722   7921935 975071285 
L5+9+BKZ-12+BB       9849     66385    753475   8048495 841009395 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.22      14.42      12.22 
L5+9+BKZ-10      36395.12      231.4      25.16       14.3      11.96 
LLL99999+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ99999-6+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-8+BB    36395.12      231.4      25.14       14.1       11.5 
BKZ99999-10+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ99999-12+BB   36395.12      231.4      25.14       14.1       11.5 
BKZ999-6+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-8+BB      36395.12      231.4      25.14       14.1       11.5 
BKZ999-10+BB     36395.12      231.4      25.14       14.1       11.5 
BKZ999-12+BB     36395.12      231.4      25.14       14.1       11.5 
L8+BKZ-10+BB     36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-10+BB   36395.12      231.4      25.14       14.1       11.5 
L5+9+BKZ-12+BB   36395.12      231.4      25.14       14.1       11.5 

Average number of calls to the recursive BB procedure `tryZ`:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             5         16        793     129459   58165667 
BKZ99999-6+BB           5         16        730      62300   15576914 
BKZ99999-8+BB           5         16        692      64362   11737003 
BKZ99999-10+BB          5         16        690      58951   12300726 
BKZ99999-12+BB          5         16        693      54522   10757205 
BKZ999-6+BB             5         16        724      66459   18443096 
BKZ999-8+BB             5         16        698      65244   13073535 
BKZ999-10+BB            5         16        686      58382   13358260 
BKZ999-12+BB            5         16        680      51461    9964228 
L8+BKZ-10+BB            5         16        687      55530   10274708 
L5+9+BKZ-10+BB          5         16        680      62255   11224780 
L5+9+BKZ-12+BB          5         16        668      57075    9674228 

Average number of visited leaves in the BB procedure:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          1          2          9 
BKZ99999-6+BB           1          1          1          2          6 
BKZ99999-8+BB           1          1          1          2          6 
BKZ99999-10+BB          1          1          1          2          5 
BKZ99999-12+BB          1          1          1          1          4 
BKZ999-6+BB             1          1          1          2          7 
BKZ999-8+BB             1          1          1          1          5 
BKZ999-10+BB            1          1          1          1          5 
BKZ999-12+BB            1          1          1          1          4 
L8+BKZ-10+BB            1          1          1          1          4 
L5+9+BKZ-10+BB          1          1          1          2          5 
L5+9+BKZ-12+BB          1          1          1          2          4 

Largest absolute value of z_j in the BB procedures:
Num. dimensions:        5         10         20         30         40   

LLL99999+BB             1          1          5         14         25 
BKZ99999-6+BB           1          1          4         11         40 
BKZ99999-8+BB           1          1          5         14         21 
BKZ99999-10+BB          1          1          4         12         45 
BKZ99999-12+BB          1          1          4          8         22 
BKZ999-6+BB             1          1          3         14         24 
BKZ999-8+BB             1          1          4         11         21 
BKZ999-10+BB            1          1          4         12         31 
BKZ999-12+BB            1          1          3          9         33 
L8+BKZ-10+BB            1          1          4          8         25 
L5+9+BKZ-10+BB          1          1          4         14         19 
L5+9+BKZ-12+BB          1          1          4         12         20 

Total time for everything: 17107.16821 seconds

