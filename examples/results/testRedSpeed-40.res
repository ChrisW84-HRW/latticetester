
========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  953      3358     10640     20013     33378 
LLL99999              824      4499     32762     91721    160783 
BKZ99999-10           713      4385     39077    135529    339142 
L5+9+BKZ-10           853      3879     25837    100563    281666 
LLL99999+BB           882      5170     41001    229244  45531188 
BKZ99999-6+BB         774      5140     43546    204185  14367102 
BKZ99999-8+BB         702      4982     44491    210964  12351859 
BKZ99999-10+BB        679      4849     45260    212129  10683752 
BKZ99999-12+BB        676      4711     45751    222903   9896012 
BKZ999-6+BB           675      4592     43902    207671  13450940 
BKZ999-8+BB           668      4556     44876    212578  13887090 
BKZ999-10+BB          666      4500     45288    218465  11521399 
BKZ999-12+BB          661      4416     45754    218982  10938427 
L8+BKZ-10+BB          951      4916     35659    181245  13168286 
L5+9+BKZ-10+BB        943      4545     31799    178916  12816808 
L5+9+BKZ-12+BB        796      4302     31924    192589  10496292 

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

Total time for everything: 183.472758 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = double
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                  658      2122      5351      7448     11200 
LLL99999              577      2522     13243     32444     46542 
BKZ99999-10           568      2498     17400     67854    168313 
L5+9+BKZ-10           662      2776     17491     66511    162156 
LLL99999+BB           641      3177     19337    208866  73829682 
BKZ99999-6+BB         605      3130     22415    145056  19956686 
BKZ99999-8+BB         554      2979     23338    150444  15361535 
BKZ99999-10+BB        545      2873     23757    155579  16326048 
BKZ99999-12+BB        526      2829     24363    153148  14245856 
BKZ999-6+BB           514      2781     22419    149050  21910016 
BKZ999-8+BB           499      2778     23084    153446  16499192 
BKZ999-10+BB          490      2767     23617    157149  17551249 
BKZ999-12+BB          483      2724     24427    151007  13238024 
L8+BKZ-10+BB          640      3328     22880    147975  14035185 
L5+9+BKZ-10+BB        715      3498     23810    157815  13870062 
L5+9+BKZ-12+BB        660      3350     24106    161311  12135620 

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
BKZ99999-10+BB          5         16        701      57957   12646861 
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

Total time for everything: 251.911254 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = xdouble
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 4045     17173     62829    109917    170264 
LLL99999             4014     26994    218726    625799   1047157 
BKZ99999-10          3958     27751    258606    953582   2387281 
L5+9+BKZ-10          3829     21424    158917    680682   1891515 
LLL99999+BB          3966     28179    232087   1730965 465249578 
BKZ99999-6+BB        3855     29045    260504   1445876 145103805 
BKZ99999-8+BB        3710     28909    265854   1498872 124773771 
BKZ99999-10+BB       3616     28883    271118   1492993 107472118 
BKZ99999-12+BB       3519     28850    274530   1575038  99809490 
BKZ999-6+BB          3466     28711    259946   1480887 136914974 
BKZ999-8+BB          3423     28797    268580   1514611 140310704 
BKZ999-10+BB         3390     28799    272066   1544040 116355818 
BKZ999-12+BB         3365     28696    275535   1548746 109987704 
L8+BKZ-10+BB         3911     25924    205226   1296692 119612225 
L5+9+BKZ-10+BB       3905     22823    171483   1271940 147061372 
L5+9+BKZ-12+BB       3722     22450    174381   1289876 108459400 

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
BKZ999-6+BB             5         16        553      54488   10298185 
BKZ999-8+BB             5         16        559      53281   10650118 
BKZ999-10+BB            5         16        577      46995    8732897 
BKZ999-12+BB            5         16        577      46194    8237243 
L8+BKZ-10+BB            5         16        573      45842    8987982 
L5+9+BKZ-10+BB          5         16        562      49378   11235264 
L5+9+BKZ-12+BB          5         16        564      45733    8170195 

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

Total time for everything: 1850.803672 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = xdouble
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 2420      8558     25394     35299     46760 
LLL99999             2223     11461     74210    199888    296283 
BKZ99999-10          2265     12111    100565    430502   1147841 
L5+9+BKZ-10          2357     11914     91274    392917   1053168 
LLL99999+BB          2306     12437     88087   1700496 748787181 
BKZ99999-6+BB        2217     13004    104538   1052341 202923667 
BKZ99999-8+BB        2127     12911    109343   1092926 154627325 
BKZ99999-10+BB       2062     12769    113251   1113349 164322154 
BKZ99999-12+BB       2035     12721    117700   1084511 142760790 
BKZ999-6+BB          2016     12537    104674   1089136 222491514 
BKZ999-8+BB          2007     12558    109518   1122862 167780585 
BKZ999-10+BB         1988     12570    113129   1127002 176601975 
BKZ999-12+BB         1979     12582    117550   1072598 132539981 
L8+BKZ-10+BB         2269     12893    103195   1082609 149286873 
L5+9+BKZ-10+BB       2433     13013    103672   1098536 165495980 
L5+9+BKZ-12+BB       2351     12773    106291   1122938 121840481 

Average square length of shortest basis vector:
Num. dimensions:        5         10         20         30         40   

LLL5              36783.7        255      46.82      40.64      36.38 
LLL99999         36395.12     231.94      25.38      14.84       13.5 
BKZ99999-10      36395.12      231.4      25.16      14.28      11.94 
L5+9+BKZ-10      36395.12      231.4      25.14      14.34      11.96 
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
BKZ99999-8+BB           5         16        684      61897   11843113 
BKZ99999-10+BB          5         16        690      57957   12646318 
BKZ99999-12+BB          5         16        693      53608   10942933 
BKZ999-6+BB             5         16        716      66232   17066547 
BKZ999-8+BB             5         16        690      64671   12828353 
BKZ999-10+BB            5         16        686      58726   13659907 
BKZ999-12+BB            5         16        680      53002   10134442 
L8+BKZ-10+BB            5         16        667      59441   11457671 
L5+9+BKZ-10+BB          5         16        678      60102   12710804 
L5+9+BKZ-12+BB          5         16        673      58153    9279093 

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
L5+9+BKZ-10+BB          1          1          1          2          6 
L5+9+BKZ-12+BB          1          1          1          2          5 

Total time for everything: 2568.754358 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 3035     17297     85458    187368    338946 
LLL99999             3329     28918    284797    888567   1686514 
BKZ99999-10          3363     30398    343877   1385886   3817393 
L5+9+BKZ-10          3330     22209    213540   1003115   3068463 
LLL99999+BB          3524     30167    297563   1819647 386302887 
BKZ99999-6+BB        3469     31338    339346   1719953 123290614 
BKZ99999-8+BB        3404     31486    347300   1793676 105250741 
BKZ99999-10+BB       3383     31467    355323   1844843  91467179 
BKZ99999-12+BB       3361     31461    361609   1977409  85048003 
BKZ999-6+BB          3356     31060    338210   1750490 115089644 
BKZ999-8+BB          3349     31229    350938   1813853 118804566 
BKZ999-10+BB         3336     31269    356699   1905117  98681275 
BKZ999-12+BB         3340     31253    362360   1915684  93768288 
L8+BKZ-10+BB         3603     27397    268810   1571112 101140206 
L5+9+BKZ-10+BB       3508     23371    225282   1510753 123809668 
L5+9+BKZ-12+BB       3396     23307    230306   1585827  91753185 

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

Total time for everything: 1573.284628 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = quad_float
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 1688      6303     19350     29088     42910 
LLL99999             1651      8948     61732    165202    241830 
BKZ99999-10          1756     10050     93476    413390   1148711 
L5+9+BKZ-10          1954      9737     85535    381094   1062131 
LLL99999+BB          1848      9996     74318   1436015 627796374 
BKZ99999-6+BB        1850     10800     93078    911422 169582542 
BKZ99999-8+BB        1803     10911     99221    955622 134917467 
BKZ99999-10+BB       1787     10827    105101    997571 137852780 
BKZ99999-12+BB       1765     10738    111894    996706 119018978 
BKZ999-6+BB          1756     10554     92950    940372 185376845 
BKZ999-8+BB          1750     10635     99242    988281 142445398 
BKZ999-10+BB         1745     10650    104874   1012592 145312260 
BKZ999-12+BB         1740     10628    112061    980950 115174962 
L8+BKZ-10+BB         1967     10765     97184    972110 125854871 
L5+9+BKZ-10+BB       2106     10748     97094    982891 140498632 
L5+9+BKZ-12+BB       2064     10683    101472   1025504 104138554 

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

Total time for everything: 2165.413662 seconds


========================================================================
Compare different reduction strategies for different real types.
****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 50
PRIMAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                13854     79949    371409    678761   1061260 
LLL99999            15528    136750   1345708   4099084   7286741 
BKZ99999-10         15366    142284   1603932   6402358  17287035 
L5+9+BKZ-10         14437     99219    961444   4655642  14264100 
LLL99999+BB         15967    139832   1397135  10859862 3049167784 
BKZ99999-6+BB       15646    144770   1585693   9208761 951021601 
BKZ99999-8+BB       15553    145436   1617965   9543939 806052802 
BKZ99999-10+BB      15416    145502   1651531   9594035 693951227 
BKZ99999-12+BB      15407    145540   1665896  10147657 643278203 
BKZ999-6+BB         15330    144516   1580123   9446301 880758990 
BKZ999-8+BB         15325    145184   1633165   9719170 905744565 
BKZ999-10+BB        15250    145575   1655547   9936180 748860027 
BKZ999-12+BB        15266    145630   1672171   9921546 712913554 
L8+BKZ-10+BB        15744    124146   1219611   8214792 865692125 
L5+9+BKZ-10+BB      14963    102668   1007987   7985385 843391360 
L5+9+BKZ-12+BB      14629    102622   1021316   8636733 684483991 

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

Total time for everything: 11978.6204 seconds

****************************************************
TestReducersSpeed with m = 1099511627791
Types: Int = NTL::ZZ, Real = NTL::RR
Number of replications (different multipliers a): 50
DUAL lattice,  Norm: L2NORM,  Decomposition: CHOLESKY.

Computing times in microseconds:
Num. dimensions:       5        10        20        30        40  

LLL5                 8921     42978    153277    229367    339620 
LLL99999             9222     67313    564784   1581336   2346651 
BKZ99999-10          9207     72250    786190   3668392  10480849 
L5+9+BKZ-10          9424     63198    679226   3438149   9726791 
LLL99999+BB          9649     70466    625305  10860190 5015246979 
BKZ99999-6+BB        9538     75239    778654   7157284 1335835858 
BKZ99999-8+BB        9334     75520    815026   7732617 1002873181 
BKZ99999-10+BB       9262     75619    839732   7868370 1055914570 
BKZ99999-12+BB       9262     75515    865215   7663929 916316013 
BKZ999-6+BB          9257     74875    778114   7370521 1588438390 
BKZ999-8+BB          9122     75029    814603   7796078 1119907032 
BKZ999-10+BB         9128     75161    838423   7926347 1145449506 
BKZ999-12+BB         9140     75213    864819   7497064 847713047 
L8+BKZ-10+BB         9808     70614    734219   7118721 876711067 
L5+9+BKZ-10+BB      10009     66092    731762   7855054 959404513 
L5+9+BKZ-12+BB       9872     66217    751933   7985934 827608348 

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

Total time for everything: 16831.09101 seconds

