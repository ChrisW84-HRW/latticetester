****************************************************
Types: Int = long, Real = double
Tests in dual lattice?  1

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                  370      1878      9852     25307     49767 
LLL99999              323      2384     24203     74287    135562 
BKZ99999-10           360      2633     36217    178208    409551 
L5+8+99+BKZ-10        416      2740     34029    166709    382358 
LLL5+BB               289      1807     18063    700169  23589291 
LLL8+BB               298      2163     19534    386764  18628396 
LLL99999+BB           288      2318     25966    241766  10151574 
BKZ99999-6+BB         326      2535     34053    251910   6602922 
BKZ99999-8+BB         311      2524     36549    265198   5872619 
BKZ99999-10+BB        302      2456     37678    274214   5443500 
BKZ99999-12+BB        303      2444     39681    308496   4994020 
BKZ999-6+BB           306      2423     34975    246461   6678903 
BKZ999-8+BB           300      2466     36499    252181   5505468 
BKZ999-10+BB          304      2433     37637    274477   5350550 
BKZ999-12+BB          305      2425     39937    315856   5316444 
L8+BKZ-10+BB          333      2594     37061    270520   4832504 
L5+8+99+BKZ-10+BB        377      2667     35789    258438   4328284 
L5+8+99+BKZ-12+BB        379      2488     37194    280478   5307214 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5                 7478        840        456        411        391 
LLL99999             7478        759        341        294        278 
BKZ99999-10          7478        757        337        287        274 
L5+8+99+BKZ-10       7478        757        339        286        268 
LLL5+BB              7478        757        337        284        266 
LLL8+BB              7478        757        337        284        266 
LLL99999+BB          7478        757        337        284        266 
BKZ99999-6+BB        7478        757        337        284        266 
BKZ99999-8+BB        7478        757        337        284        266 
BKZ99999-10+BB       7478        757        337        284        266 
BKZ99999-12+BB       7478        757        337        284        266 
BKZ999-6+BB          7478        757        337        284        266 
BKZ999-8+BB          7478        757        337        284        266 
BKZ999-10+BB         7478        757        337        284        266 
BKZ999-12+BB         7478        757        337        284        266 
L8+BKZ-10+BB         7478        757        337        284        266 
L5+8+99+BKZ-10+BB       7478        757        337        284        266 
L5+8+99+BKZ-12+BB       7478        757        337        284        266 

Total time for everything: 118.987736 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = double
Tests in dual lattice?  1

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                  473      1570      3848      6906     11158 
LLL99999              424      2099     12272     24873     32697 
BKZ99999-10           429      2320     21619     83455    179293 
L5+8+99+BKZ-10        625      2956     24178     89637    189982 
LLL5+BB               533      2390     18894    817527  27486558 
LLL8+BB               511      2668     15979    431168  21799805 
LLL99999+BB           481      2796     18637    234918  11955397 
BKZ99999-6+BB         491      2932     24793    206455   7557285 
BKZ99999-8+BB         458      2814     26470    204459   6543208 
BKZ99999-10+BB        439      2717     27480    205999   6014882 
BKZ99999-12+BB        435      2680     29038    222540   5502030 
BKZ999-6+BB           435      2659     25197    195881   7623976 
BKZ999-8+BB           434      2669     26377    194958   6056208 
BKZ999-10+BB          432      2589     27462    205869   5924767 
BKZ999-12+BB          427      2587     29186    231994   5873044 
L8+BKZ-10+BB          542      3115     28131    205531   5252097 
L5+8+99+BKZ-10+BB        679      3607     29662    208627   5105157 
L5+8+99+BKZ-12+BB        624      3445     31139    221677   6181972 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5                 7478        840        456        411        391 
LLL99999             7478        759        341        294        278 
BKZ99999-10          7478        757        337        287        274 
L5+8+99+BKZ-10       7478        757        338        286        269 
LLL5+BB              7478        757        337        284        266 
LLL8+BB              7478        757        337        284        266 
LLL99999+BB          7478        757        337        284        266 
BKZ99999-6+BB        7478        757        337        284        266 
BKZ99999-8+BB        7478        757        337        284        266 
BKZ99999-10+BB       7478        757        337        284        266 
BKZ99999-12+BB       7478        757        337        284        266 
BKZ999-6+BB          7478        757        337        284        266 
BKZ999-8+BB          7478        757        337        284        266 
BKZ999-10+BB         7478        757        337        284        266 
BKZ999-12+BB         7478        757        337        284        266 
L8+BKZ-10+BB         7478        757        337        284        266 
L5+8+99+BKZ-10+BB       7478        757        337        284        266 
L5+8+99+BKZ-12+BB       7478        757        337        284        266 

Total time for everything: 133.853426 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = xdouble
Tests in dual lattice?  1

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                 1429      5040     11645     19117     29463 
LLL99999             1329      7081     42118     81511    100127 
BKZ99999-10          1350      7827     72876    261810    497326 
L5+8+99+BKZ-10       1732      9025     75595    297123    553044 
LLL5+BB              1404      6246     69952   4544978 153067883 
LLL8+BB              1397      7356     50351   2364217 122901001 
LLL99999+BB          1374      8027     58641   1202230  66633765 
BKZ99999-6+BB        1407      8657     76543    939989  41996294 
BKZ99999-8+BB        1334      8671     83037    904585  35486605 
BKZ99999-10+BB       1299      8415     87717    894118  32307203 
BKZ99999-12+BB       1285      8324     93749    934719  29949098 
BKZ999-6+BB          1267      8085     77661    873870  42020068 
BKZ999-8+BB          1259      8231     83177    865163  33799519 
BKZ999-10+BB         1245      8172     87762    895623  32364440 
BKZ999-12+BB         1232      8136     95675    979332  32199393 
L8+BKZ-10+BB         1469      9012     89740    904426  30485897 
L5+8+99+BKZ-10+BB       1813     10125     90084    951591  34726902 
L5+8+99+BKZ-12+BB       1729      9851     98383    988273  31209800 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5                 7478        840        456        411        391 
LLL99999             7478        759        341        294        278 
BKZ99999-10          7478        757        337        287        274 
L5+8+99+BKZ-10       7478        757        337        287        269 
LLL5+BB              7478        757        337        284        266 
LLL8+BB              7478        757        337        284        266 
LLL99999+BB          7478        757        337        284        266 
BKZ99999-6+BB        7478        757        337        284        266 
BKZ99999-8+BB        7478        757        337        284        266 
BKZ99999-10+BB       7478        757        337        284        266 
BKZ99999-12+BB       7478        757        337        284        266 
BKZ999-6+BB          7478        757        337        284        266 
BKZ999-8+BB          7478        757        337        284        266 
BKZ999-10+BB         7478        757        337        284        266 
BKZ999-12+BB         7478        757        337        284        266 
L8+BKZ-10+BB         7478        757        337        284        266 
L5+8+99+BKZ-10+BB       7478        757        337        284        266 
L5+8+99+BKZ-12+BB       7478        757        337        284        266 

Total time for everything: 740.848011 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


****************************************************
Types: Int = NTL::ZZ, Real = quad_float
Tests in dual lattice?  1

TestReducersSpeed with m = 1048573
Results for `testLoop`
Timings (in microseconds) for different methods for 50 replications 

Num. dimensions:       5        10        20        30        40  

LLL5                  816      2972      7380     13602     22919 
LLL99999              805      4164     25801     50978     65021 
BKZ99999-10           887      4977     51985    193206    400741 
L5+8+99+BKZ-10       1218      6169     53846    223552    451822 
LLL5+BB               937      4010     53450   3512682 119280764 
LLL8+BB               922      4602     35612   1800485  95423253 
LLL99999+BB           926      5057     39448    915086  52076296 
BKZ99999-6+BB         978      5679     53847    721064  32987041 
BKZ99999-8+BB         934      5732     59140    664516  26988956 
BKZ99999-10+BB        922      5617     63208    679147  26333944 
BKZ99999-12+BB        915      5552     68639    707364  23506907 
BKZ999-6+BB           910      5349     54587    684983  33270870 
BKZ999-8+BB           908      5544     59030    659401  24891051 
BKZ999-10+BB          905      5480     64389    683104  28057129 
BKZ999-12+BB          899      5483     70774    737531  24540817 
L8+BKZ-10+BB         1086      6224     64866    688991  23241437 
L5+8+99+BKZ-10+BB       1381      7161     65753    753815  26580812 
L5+8+99+BKZ-12+BB       1315      6999     72738    771169  24131337 

Sums of square lengths of shortest basis vector:
Num. dimensions:          5         10         20         30         40 

LLL5                 7478        840        456        411        391 
LLL99999             7478        759        341        294        278 
BKZ99999-10          7478        757        337        287        274 
L5+8+99+BKZ-10       7478        757        337        287        270 
LLL5+BB              7478        757        337        284        266 
LLL8+BB              7478        757        337        284        266 
LLL99999+BB          7478        757        337        284        266 
BKZ99999-6+BB        7478        757        337        284        266 
BKZ99999-8+BB        7478        757        337        284        266 
BKZ99999-10+BB       7478        757        337        284        266 
BKZ99999-12+BB       7478        757        337        284        266 
BKZ999-6+BB          7478        757        337        284        266 
BKZ999-8+BB          7478        757        337        284        266 
BKZ999-10+BB         7478        757        337        284        266 
BKZ999-12+BB         7478        757        337        284        266 
L8+BKZ-10+BB         7478        757        337        284        266 
L5+8+99+BKZ-10+BB       7478        757        337        284        266 
L5+8+99+BKZ-12+BB       7478        757        337        284        266 

Total time for everything: 577.885943 seconds

We see that LLL or BKZ alone do not always find a shortest vector.


