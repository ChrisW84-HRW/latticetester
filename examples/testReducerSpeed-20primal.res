****************************************************
Types: Int = long, Real = double
Tests in dual lattice?  0

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                  406      1823      7884     16909     29011 
LLL99999              382      2773     29426     87017    162156 
BKZ99999-10           415      2866     42071    195892    404089 
L5+8+99+BKZ-10        452      2684     32381    179550    340436 
LLL5+BB               323      1730     15096   1423158   3232433 
LLL8+BB               363      2282     18232    303681    611323 
LLL99999+BB           349      2585     29665    255547    544968 
BKZ99999-6+BB         368      2797     37799    257904    612253 
BKZ99999-8+BB         347      2742     40281    281775    666058 
BKZ99999-10+BB        352      2711     41286    296508    688237 
BKZ99999-12+BB        350      2691     42632    323915    712748 
BKZ999-6+BB           350      2675     37928    262313    617932 
BKZ999-8+BB           348      2680     40243    265174    635983 
BKZ999-10+BB          348      2659     41658    284468    673112 
BKZ999-12+BB          353      2646     43584    315419    742233 
L8+BKZ-10+BB          391      2747     34671    265368    634793 
L5+8+99+BKZ-10+BB        413      2646     32673    274628    642671 
L5+8+99+BKZ-12+BB        403      2554     35611    298064    695012 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5           1.227149089e+11 3.398566431e+12 2.88017716e+13 5.462741053e+13 5.497526682e+13 
LLL99999       1.22211526e+11 3.094930539e+12 2.214367071e+13 5.071462778e+13 5.497526682e+13 
BKZ99999-10    1.22211526e+11 3.091616709e+12 2.189785757e+13 4.932703092e+13 5.497526682e+13 
L5+8+99+BKZ-10 1.22211526e+11 3.091616709e+12 2.191217234e+13 4.904064518e+13 5.497526682e+13 
LLL5+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL8+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL99999+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-6+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-8+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-6+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-8+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-12+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L8+BKZ-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 

Total time for everything: 18.891863 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = double
Tests in dual lattice?  0

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                  550      1898      5969     11750     18965 
LLL99999              491      2731     18457     42935     67177 
BKZ99999-10           466      2866     28474    105867    188163 
L5+8+99+BKZ-10        640      3091     25152    103661    176641 
LLL5+BB               566      2738     19350   1686579   3837485 
LLL8+BB               598      3254     19724    349482    698978 
LLL99999+BB           569      3503     25994    256447    549772 
BKZ99999-6+BB         547      3691     32469    222010    541578 
BKZ99999-8+BB         515      3557     34671    236389    562829 
BKZ99999-10+BB        543      3512     36068    240356    559628 
BKZ99999-12+BB        514      3478     37279    255651    555585 
BKZ999-6+BB           520      3526     33018    228257    540598 
BKZ999-8+BB           501      3515     34866    222499    531374 
BKZ999-10+BB          503      3534     35910    231841    553639 
BKZ999-12+BB          509      3517     37516    251994    567382 
L8+BKZ-10+BB          613      3807     32097    225674    529202 
L5+8+99+BKZ-10+BB        707      3959     32373    232910    567795 
L5+8+99+BKZ-12+BB        654      3760     34711    249905    561857 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5           1.227149089e+11 3.398566431e+12 2.88017716e+13 5.462741053e+13 5.497526682e+13 
LLL99999       1.22211526e+11 3.094930539e+12 2.214367071e+13 5.071462778e+13 5.497526682e+13 
BKZ99999-10    1.22211526e+11 3.091616709e+12 2.189785757e+13 4.932703092e+13 5.497526682e+13 
L5+8+99+BKZ-10 1.22211526e+11 3.091616709e+12 2.191217234e+13 4.904064518e+13 5.497526682e+13 
LLL5+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL8+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL99999+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-6+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-8+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-6+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-8+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-12+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L8+BKZ-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 

Total time for everything: 17.363609 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = xdouble
Tests in dual lattice?  0

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                 1685      6548     18233     31770     47876 
LLL99999             1708      9770     64878    138586    202430 
BKZ99999-10          1646     10487     96594    328907    535140 
L5+8+99+BKZ-10       1994     10394     87457    310607    484236 
LLL5+BB              1613      7618     71020   9627791  21729571 
LLL8+BB              1799      9490     59224   1900037   3846401 
LLL99999+BB          1745     10819     80633   1286790   2831336 
BKZ99999-6+BB        1679     11405    100414    979462   2515232 
BKZ99999-8+BB        1603     11333    107155   1027834   2508769 
BKZ99999-10+BB       1573     11221    111684   1021425   2487699 
BKZ99999-12+BB       1561     11166    118228   1072138   2436870 
BKZ999-6+BB          1552     10968    100877   1014360   2500537 
BKZ999-8+BB          1529     11126    107278    961704   2377085 
BKZ999-10+BB         1531     11007    112156    993250   2474422 
BKZ999-12+BB         1527     11078    119930   1073077   2436270 
L8+BKZ-10+BB         1846     11198    101797    937472   2337489 
L5+8+99+BKZ-10+BB       2065     11560    102413    990681   2524500 
L5+8+99+BKZ-12+BB       1942     11286    108105   1073308   2437624 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5           1.227149089e+11 3.398566431e+12 2.88017716e+13 5.462741053e+13 5.497526682e+13 
LLL99999       1.22211526e+11 3.094930539e+12 2.214367071e+13 5.071462778e+13 5.497526682e+13 
BKZ99999-10    1.22211526e+11 3.091616709e+12 2.189785757e+13 4.932703092e+13 5.497526682e+13 
L5+8+99+BKZ-10 1.22211526e+11 3.091616709e+12 2.185897891e+13 4.935188661e+13 5.497526682e+13 
LLL5+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL8+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL99999+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-6+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-8+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-6+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-8+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-12+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L8+BKZ-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 

Total time for everything: 83.377389 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = quad_float
Tests in dual lattice?  0

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                  918      3764     12754     26626     44247 
LLL99999              980      5826     40792     97798    159972 
BKZ99999-10          1009      6510     67650    244390    445010 
L5+8+99+BKZ-10       1307      6864     62403    232387    411112 
LLL5+BB              1014      4891     54153   7410708  16950570 
LLL8+BB              1110      6073     41695   1415200   2860386 
LLL99999+BB          1097      6775     54515    965724   2136043 
BKZ99999-6+BB        1114      7380     69027    742511   1932848 
BKZ99999-8+BB        1088      7336     74888    773858   1956365 
BKZ99999-10+BB       1076      7253     79438    774823   1929918 
BKZ99999-12+BB       1076      7204     84869    810819   1895828 
BKZ999-6+BB          1075      7019     69667    767351   1921081 
BKZ999-8+BB          1084      7124     75093    724704   1831491 
BKZ999-10+BB         1074      7138     79517    750484   1919827 
BKZ999-12+BB         1070      7137     85792    811268   1895148 
L8+BKZ-10+BB         1253      7372     73466    711255   1809507 
L5+8+99+BKZ-10+BB       1446      7790     74779    754222   1969836 
L5+8+99+BKZ-12+BB       1398      7651     79980    814586   1878662 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5           1.227149089e+11 3.398566431e+12 2.88017716e+13 5.462741053e+13 5.497526682e+13 
LLL99999       1.22211526e+11 3.094930539e+12 2.214367071e+13 5.071462778e+13 5.497526682e+13 
BKZ99999-10    1.22211526e+11 3.091616709e+12 2.189785757e+13 4.932703092e+13 5.497526682e+13 
L5+8+99+BKZ-10 1.22211526e+11 3.091616709e+12 2.185897891e+13 4.933472619e+13 5.497526682e+13 
LLL5+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL8+BB        1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
LLL99999+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-6+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-8+BB  1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ99999-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-6+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-8+BB    1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
BKZ999-12+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L8+BKZ-10+BB   1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-10+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 
L5+8+99+BKZ-12+BB 1.22211526e+11 3.091616709e+12 2.185805651e+13 4.862809212e+13 5.497526682e+13 

Total time for everything: 64.105159 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


