# Created using:
#  Subdivfit -mf - -fi mechpart.4102.pts -crep 1e-5 -csharp 2e-06 -recon
# Initial mesh:
# Genus: c=1 b=0  v=175 f=358 e=537  genus=3  sharpe=143 cuspv=0
# 4102 points read
# Starting reconstruction sequence
#  crep=1e-05, csharp=2e-06
#  internal xform: F 1  1 0 0  0 1 0  0 0 1  0 0 0  0
# 
# (Timing on cpu=alpha host=opal)
# fgfit_before: v=175 nse=143/537  edis=0.00661521 etot=0.00865121
# fgfit_after : v=175 nse=143/537  edis=0.000179099 etot=0.0022151
#  (_fgfit:                309.82)
# 
# Stoc, crep=1e-05 csharp=2e-06 wcrep=1e-05 wcsharp=2e-06
# stoc_before: v=175 nse=143/537  edis=0.000179099 etot=0.0022151
# it 727, last search: 42 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=125 f=258 e=387  genus=3  sharpe=122 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts       727         0         0         0
#              success        50         0         0         0
#      positive_energy       188         0         0         0
#         bad_dihedral        89         0         0         0
#            bad_sharp       380         0         0         0
#         illegal_move        20         0         0         0
#  last cedis=0.000247864 cetot=0.00174186
# stoc_after : v=125 nse=122/387  edis=0.000247847 etot=0.00174185
#  (_stoc:                1266.27)
# 
# fgfit_before: v=125 nse=122/387  edis=0.000247854 etot=0.00174185
# fgfit_after : v=125 nse=122/387  edis=0.000190165 etot=0.00168417
#  (_fgfit:                 82.13)
# 
# Stoc, crep=1e-05 csharp=2e-06 wcrep=1e-05 wcsharp=2e-06
# stoc_before: v=125 nse=122/387  edis=0.000190166 etot=0.00168417
# it 428, last search: 139 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=125 f=258 e=387  genus=3  sharpe=121 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts         0         0         0       428
#              success         0         0         0         5
#      positive_energy         0         0         0       186
#         bad_dihedral         0         0         0         0
#            bad_sharp         0         0         0       237
#         illegal_move         0         0         0         0
#  last cedis=0.000163021 cetot=0.00165502
# stoc_after : v=125 nse=121/387  edis=0.000164294 etot=0.00165629
#  (_stoc:                 613.10)
# 
# fgfit_before: v=125 nse=121/387  edis=0.000164294 etot=0.00165629
# fgfit_after : v=125 nse=121/387  edis=0.000138871 etot=0.00163087
#  (_fgfit:                 85.13)
# 
# Stoc, crep=1e-05 csharp=2e-06 wcrep=1e-05 wcsharp=2e-06
# stoc_before: v=125 nse=121/387  edis=0.000138871 etot=0.00163087
# it 518, last search: 35 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=122 f=252 e=378  genus=3  sharpe=115 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts       518         0       510       515
#              success         3         0         6         5
#      positive_energy       241         0       270       224
#         bad_dihedral       258         0       216         0
#            bad_sharp         0         0         0       286
#         illegal_move        16         0        18         0
#  last cedis=0.00013745 cetot=0.00158745
# stoc_after : v=122 nse=115/378  edis=0.000137191 etot=0.00158719
#  (_stoc:                3010.58)
# 
# fgfit_before: v=122 nse=115/378  edis=0.000137189 etot=0.00158719
# fgfit_after : v=122 nse=115/378  edis=0.000122075 etot=0.00157207
#  (_fgfit:                 79.33)
#  (Subdivfit:            5446.88)
# Summary of timers (cpu=alpha host=opal):
#  ___submesh:         (1535  )    0.48:4.47     av=     1.29   sum=  1974.12
#  ____gmakeSpatial:   (10    )    4.43:5.67     av=     4.69   sum=    46.95
#  ____gspatialproject:(10    )    1.32:1.58     av=     1.39   sum=    13.87
#  ___gallproject:     (10    )    5.80:7.28     av=     6.11   sum=    61.15
#  ___computegrad:     (268   )    0.48:0.53     av=     0.51   sum=   136.73
#  ___gneighproject:   (463   )    0.73:0.83     av=     0.77   sum=   356.08
#  __fgfit_iter:       (70    )    5.82:9.52     av=     7.38   sum=   516.27
#  _fgfit:             (4     )   79.33:309.82   av=   139.10   sum=   556.42
#  _initial_fit:       (3     )    8.03:10.05    av=     8.71   sum=    26.13
#  ____lmakespatial:   (1525  )    0.00:0.52     av=     0.11   sum=   169.00
#  ____lspatialproject:(1525  )    0.00:0.65     av=     0.07   sum=   100.38
#  ___lallproject:     (1525  )    0.00:1.17     av=     0.18   sum=   270.80
#  ____lneighproject:  (5737  )    0.05:0.72     av=     0.26   sum=  1473.85
#  ____lcombinations:  (5737  )    0.00:0.12     av=     0.03   sum=   173.35
#  _____qrdLLS:        (5737  )    0.00:0.47     av=     0.05   sum=   275.78
#  ____lsolve:         (5737  )    0.00:0.52     av=     0.06   sum=   318.27
#  ___loptimize:       (5737  )    0.07:1.43     av=     0.37   sum=  2115.62
#  __tryecol:          (829   )    1.45:13.02    av=     3.61   sum=  2996.23
#  __lattempt:         (1673  )    0.00:14.20    av=     2.89   sum=  4839.77
#  _stoc:              (3     )  613.10:3010.58  av=  1629.98   sum=  4889.95
#  __tryesha:          (420   )    1.78:6.87     av=     3.17   sum=  1331.63
#  __tryeswa:          (276   )    1.00:5.73     av=     1.85   sum=   510.35
#  Subdivfit:          (1     )        :         av=  5446.88   sum=  5446.88
# Summary of statistics:
# Sospobcells:(177472 )          1:12          av=3.08456     sd=1.4458
# Sssncellsv: (183204 )          1:216         av=5.22496     sd=13.4771
# Sssnelemsv: (183204 )          1:121         av=14.4881     sd=12.2257
# Sospcelln:  (76494  )          1:41          av=7.15642     sd=5.13886
# Sprojquick: (9406000)          0:1           av=0.475498    sd=0.499399
# Sprojnei:   (4933467)          0:6           av=0.0586869   sd=0.251464
# Sprojf:     (4933467)          4:51          av=11.4561     sd=2.92603
# Sprojunexp: (273288 )          0:1           av=0.00454466  sd=0.0672608
# Sneval:     (70     )          6:9           av=6.61429     sd=0.728167
# Stmin:      (70     )  0.0205714:0.096       av=0.0568352   sd=0.0134613
# Ssetfrem:   (1525   )        416:1936        av=1219.92     sd=220.985
# Ssetvrem:   (1525   )        218:949         av=613.313     sd=105.526
# Secolpts:   (829    )        455:3306        av=1511.67     sd=514.268
# Secolmf:    (829    )         41:206         av=101.169     sd=30.492
# Secolmv:    (829    )          5:24          av=10.3752     sd=3.32686
# Secolsmv:   (829    )        359:1720        av=874.534     sd=255.207
# Soptnit:    (1185   )          2:10          av=3.59156     sd=1.98651
# Soptnig:    (76     )         11:12          av=11.0526     sd=0.224781
# Sechange:   (69     )-1.86599e-05:-5.2348e-07 av=-7.91015e-06 sd=5.06318e-06
# Seshapts:   (420    )        632:3306        av=1524.41     sd=486.673
# Seshamf:    (420    )         49:208         av=100.743     sd=30.0637
# Seshamv:    (420    )          6:25          av=10.9381     sd=3.55737
# Seshasmv:   (420    )        431:1736        av=871.331     sd=250.503
# Seswapts:   (276    )        347:1844        av=918.605     sd=356.011
# Seswamf:    (276    )         29:114         av=61.7283     sd=18.1453
# Seswamv:    (276    )          4:4           av=4           sd=0
# Seswasmv:   (276    )        255:1002        av=544.105     sd=159.075
Vertex 1  0.475623 0.46925 0.403842
Vertex 2  0.366657 0.412934 0.448751
Vertex 5  0.405249 0.408255 0.454294
Vertex 8  0.366677 0.366647 0.633334
Vertex 10  0.366669 0.46152 0.408913
Vertex 14  0.419064 0.366672 0.442207
Vertex 20  0.366661 0.366653 0.366652
Vertex 21  0.366711 0.400553 0.514805
Vertex 24  0.366684 0.510203 0.366723
Vertex 33  0.402936 0.366578 0.507073
Vertex 35  0.420072 0.477124 0.415236
Vertex 43  0.366631 0.530948 0.399529
Vertex 47  0.42109 0.435291 0.366577
Vertex 48  0.40088 0.474147 0.366416
Vertex 54  0.478066 0.513684 0.409188
Vertex 58  0.366623 0.427018 0.562425
Vertex 64  0.456254 0.406599 0.366728
Vertex 67  0.461273 0.36664 0.406299
Vertex 68  0.540147 0.366691 0.403051
Vertex 74  0.402797 0.401393 0.529294
Vertex 85  0.366688 0.586664 0.44904
Vertex 87  0.398722 0.515362 0.366803
Vertex 88  0.397956 0.464433 0.0999808
Vertex 94  0.510098 0.482672 0.408043
Vertex 99  0.42101 0.585191 0.441011
Vertex 109  0.412952 0.435095 0.580713
Vertex 110  0.366699 0.476082 0.598337
Vertex 113  0.366623 0.633339 0.36665
Vertex 118  0.39929 0.53456 0.099999
Vertex 129  0.528043 0.526134 0.406491
Vertex 130  0.471116 0.55948 0.405401
Vertex 132  0.425892 0.366732 0.573334
Vertex 133  0.453323 0.465235 0.58293
Vertex 139  0.505358 0.394025 0.366627
Vertex 140  0.510085 0.39458 0.100027
Vertex 147  0.45725 0.633302 0.409932
Vertex 149  0.39812 0.598302 0.496887
Vertex 151  0.413902 0.55546 0.366586
Vertex 160  0.633385 0.366687 0.366663
Vertex 161  0.633288 0.4445 0.410284
Vertex 164  0.589085 0.41569 0.447054
Vertex 165  0.548074 0.515791 0.402129
Vertex 178  0.412956 0.633376 0.448501
Vertex 180  0.560677 0.414693 0.366697
Vertex 194  0.366794 0.593906 0.556293
Vertex 195  0.440976 0.586864 0.366726
Vertex 203  0.586262 0.366685 0.454066
Vertex 209  0.633355 0.522408 0.401625
Vertex 213  0.52868 0.633316 0.401588
Vertex 219  0.503159 0.366618 0.605069
Vertex 220  0.503333 0.506987 0.596041
Vertex 231  0.366615 0.536718 0.591517
Vertex 233  0.395367 0.506061 0.633324
Vertex 234  0.44244 0.412798 0.633343
Vertex 235  0.407205 0.452316 0.633351
Vertex 236  0.487923 0.396536 0.633334
Vertex 255  0.633325 0.397506 0.481369
Vertex 258  0.633315 0.462279 0.366672
Vertex 260  0.633306 0.569109 0.430514
Vertex 267  0.633369 0.366687 0.633338
Vertex 273  0.404123 0.633304 0.550626
Vertex 293  0.497432 0.604488 0.36665
Vertex 311  0.592988 0.406145 0.517064
Vertex 314  0.604584 0.51257 0.3667
Vertex 326  0.461711 0.53124 0.578716
Vertex 331  0.595553 0.455665 0.36662
Vertex 339  0.366671 0.63335 0.63334
Vertex 343  0.394827 0.497863 0.90001
Vertex 359  0.430471 0.578699 0.0999981
Vertex 363  0.444879 0.407299 0.0999812
Vertex 367  0.633305 0.415464 0.55994
Vertex 374  0.589865 0.59296 0.457128
Vertex 375  0.633322 0.633326 0.366652
Vertex 379  0.579038 0.36669 0.567161
Vertex 380  0.566459 0.426794 0.576506
Vertex 381  0.549472 0.450427 0.588117
Vertex 382  0.532876 0.503015 0.5948
Vertex 388  0.422014 0.580603 0.555189
Vertex 401  0.411173 0.556426 0.63334
Vertex 406  0.542556 0.404589 0.633309
Vertex 413  0.548308 0.593307 0.366678
Vertex 425  0.599505 0.366639 0.503903
Vertex 434  0.595114 0.455775 0.10001
Vertex 439  0.575796 0.633335 0.437646
Vertex 469  0.583696 0.438184 0.63335
Vertex 479  0.585262 0.558858 0.366639
Vertex 482  0.486182 0.60577 0.100006
Vertex 495  0.633322 0.63331 0.478976
Vertex 501  0.633205 0.591146 0.470797
Vertex 508  0.633326 0.471254 0.59262
Vertex 524  0.453793 0.594224 0.633332
Vertex 533  0.55387 0.411878 0.899998
Vertex 535  0.500754 0.393898 0.900006
Vertex 558  0.604528 0.510287 0.100003
Vertex 560  0.599351 0.598993 0.531551
Vertex 566  0.633242 0.524219 0.596373
Vertex 569  0.470737 0.633334 0.598367
Vertex 574  0.56186 0.415979 0.0999396
Vertex 611  0.598869 0.633352 0.488147
Vertex 612  0.633465 0.595865 0.513684
Vertex 617  0.633274 0.574642 0.568092
Vertex 618  0.531285 0.530644 0.58616
Vertex 621  0.526447 0.633348 0.591648
Vertex 627  0.506556 0.60511 0.63333
Vertex 634  0.410264 0.447941 0.899987
Vertex 637  0.602523 0.483971 0.633321
Vertex 665  0.633324 0.540594 0.633334
Vertex 668  0.567171 0.633248 0.575039
Vertex 670  0.633334 0.633327 0.633337
Vertex 683  0.599042 0.537156 0.633339
Vertex 684  0.602926 0.483104 0.89999
Vertex 691  0.445379 0.4101 0.900004
Vertex 698  0.550781 0.594429 0.0999755
Vertex 706  0.584907 0.633373 0.534312
Vertex 717  0.561388 0.586214 0.633322
Vertex 718  0.50831 0.603455 0.900019
Vertex 733  0.585223 0.442662 0.90002
Vertex 742  0.58711 0.555276 0.0999914
Vertex 765  0.556397 0.588664 0.899991
Vertex 791  0.409829 0.557202 0.899994
Vertex 813  0.458418 0.596696 0.900002
Vertex 850  0.597563 0.542687 0.900013
Face 6  5 1 2
Face 16  1 10 2
Face 22  14 1 5
Face 32  20 8 2
Face 33  2 8 21
Face 48  2 21 5
Face 52  1 35 10
Face 64  10 43 24
Face 78  5 33 14
Face 82  8 14 33
Face 86  74 5 21
Face 87  10 35 43
Face 104  20 47 64
Face 152  8 58 21
Face 157  85 24 43
Face 161  24 87 48
Face 168  1 67 68
Face 170  14 20 67
Face 185  54 99 43
Face 205  110 109 58
Face 208  24 85 113
Face 215  48 87 88
Face 216  118 88 87
Face 228  67 20 68
Face 234  68 94 1
Face 237  129 54 1
Face 239  54 130 99
Face 243  8 33 132
Face 249  133 74 109
Face 251  43 99 85
Face 263  58 8 110
Face 281  151 118 87
Face 291  68 20 160
Face 305  1 94 129
Face 306  165 129 94
Face 316  147 99 130
Face 325  33 74 132
Face 326  133 132 74
Face 340  147 178 99
Face 341  85 99 149
Face 342  178 149 99
Face 345  180 160 139
Face 346  64 140 139
Face 396  94 161 165
Face 397  209 165 161
Face 448  110 8 231
Face 459  234 8 236
Face 496  164 255 161
Face 513  130 213 147
Face 544  160 180 258
Face 562  8 235 233
Face 567  267 236 8
Face 574  195 113 293
Face 583  195 118 151
Face 605  160 161 255
Face 611  161 258 209
Face 639  231 326 220
Face 641  178 113 273
Face 654  331 258 180
Face 668  231 339 194
Face 721  255 311 367
Face 730  314 258 331
Face 734  129 260 374
Face 739  213 375 113
Face 750  379 267 219
Face 767  231 194 326
Face 773  194 149 388
Face 792  231 8 339
Face 796  235 343 233
Face 797  401 233 343
Face 805  236 267 406
Face 813  194 339 113
Face 836  140 363 88
Face 839  311 203 425
Face 855  331 180 434
Face 862  375 213 439
Face 874  220 382 381
Face 885  326 194 388
Face 887  149 273 388
Face 906  64 363 140
Face 929  406 267 469
Face 940  413 375 479
Face 975  260 375 495
Face 981  479 258 314
Face 982  331 434 314
Face 986  374 260 501
Face 991  213 129 439
Face 992  374 439 129
Face 1010  382 508 381
Face 1014  380 311 379
Face 1015  380 381 311
Face 1036  326 388 273
Face 1068  406 469 533
Face 1072  236 406 535
Face 1085  195 293 482
Face 1099  140 88 434
Face 1103  425 267 379
Face 1117  495 375 439
Face 1126  258 479 375
Face 1134  374 501 560
Face 1149  566 508 382
Face 1155  326 273 569
Face 1157  220 326 569
Face 1178  524 339 401
Face 1244  311 381 367
Face 1257  314 558 479
Face 1278  382 220 618
Face 1279  382 618 566
Face 1280  617 566 618
Face 1287  621 220 569
Face 1288  273 339 569
Face 1297  140 574 180
Face 1316  235 634 343
Face 1327  406 533 535
Face 1352  359 118 195
Face 1389  220 621 618
Face 1390  668 618 621
Face 1421  683 637 665
Face 1434  234 236 535
Face 1470  560 706 611
Face 1471  612 617 560
Face 1477  566 617 665
Face 1485  618 560 617
Face 1488  621 670 668
Face 1498  717 670 627
Face 1499  627 524 718
Face 1521  637 683 684
Face 1557  612 495 617
Face 1562  558 742 479
Face 1572  618 668 560
Face 1591  627 718 717
Face 1609  691 535 634
Face 1626  533 469 733
Face 1638  293 413 482
Face 1652  617 495 670
Face 1654  479 742 413
Face 1682  791 524 401
Face 1721  413 742 698
Face 1744  524 791 813
Face 1772  482 413 698
Face 1785  524 813 718
Face 1802  684 533 733
Face 1811  469 684 733
Face 1881  683 717 850
Face 1895  717 765 850
Face 1927  234 235 8
Face 1929  791 401 343
Face 1930  203 160 425
Face 1931  267 425 160
Face 1940  164 161 94
Face 1952  133 109 220
Face 1953  110 220 109
Face 1977  611 706 495
Face 1978  147 213 113
Face 1981  508 367 381
Face 1987  560 501 612
Face 1988  495 612 501
Face 1992  48 47 24
Face 1996  501 260 495
Face 2011  220 380 219
Face 2012  379 219 380
Face 2014  379 311 425
Face 2015  374 611 439
Face 2016  495 439 611
Face 2019  195 482 359
Face 2024  231 220 110
Face 2039  683 665 717
Face 2040  670 717 665
Face 2043  569 339 621
Face 2044  670 621 339
Face 2047  24 47 20
Face 2053  87 24 151
Face 2059  132 219 8
Face 2060  267 8 219
Face 2063  194 113 85
Face 2067  203 311 164
Face 2068  255 164 311
Face 2075  133 219 132
Face 2078  706 560 668
Face 2086  560 611 374
Face 2096  273 149 178
Face 2097  234 535 691
Face 2107  5 74 33
Face 2109  94 68 164
Face 2110  203 164 68
Face 2118  574 140 434
Face 2123  68 160 203
Face 2129  8 20 14
Face 2138  670 665 617
Face 2140  508 267 367
Face 2141  209 260 165
Face 2142  129 165 260
Face 2143  255 367 160
Face 2144  267 160 367
Face 2146  178 147 113
Face 2148  684 765 533
Face 2153  765 684 850
Face 2158  718 813 765
Face 2160  684 469 637
Face 2162  380 220 381
Face 2165  267 665 469
Face 2166  637 469 665
Face 2174  220 219 133
Face 2175  765 813 343
Face 2176  791 343 813
Face 2177  634 535 343
Face 2184  339 273 113
Face 2185  14 67 1
Face 2187  434 482 742
Face 2188  698 742 482
Face 2189  47 48 88
Face 2195  2 10 20
Face 2196  24 20 10
Face 2199  267 508 665
Face 2200  566 665 508
Face 2203  495 706 670
Face 2204  668 670 706
Face 2207  718 765 717
Face 2212  558 314 434
Face 2214  194 85 149
Face 2219  258 375 209
Face 2220  260 209 375
Face 2226  359 482 118
Face 2235  35 1 43
Face 2236  54 43 1
Face 2237  129 213 54
Face 2238  130 54 213
Face 2239  683 850 684
Face 2241  88 118 434
Face 2242  482 434 118
Face 2244  434 180 574
Face 2246  180 139 140
Face 2247  20 64 160
Face 2248  139 160 64
Face 2250  258 161 160
Face 2251  343 535 765
Face 2252  533 765 535
Face 2255  64 47 363
Face 2256  88 363 47
Face 2258  558 434 742
Face 2259  293 113 413
Face 2260  375 413 113
Face 2261  235 234 634
Face 2262  691 634 234
Face 2263  58 109 21
Face 2264  74 21 109
Face 2265  8 233 339
Face 2266  401 339 233
Face 2269  339 524 670
Face 2270  627 670 524
Face 2271  24 113 151
Face 2272  195 151 113
Edge 791 813 {sharp}
Edge 339 670 {sharp}
Edge 113 339 {sharp}
Edge 113 375 {sharp}
Edge 566 617 {sharp}
Edge 1 35 {sharp}
Edge 1 5 {sharp}
Edge 2 10 {sharp}
Edge 2 21 {sharp}
Edge 569 621 {sharp}
Edge 343 791 {sharp}
Edge 343 634 {sharp}
Edge 683 717 {sharp}
Edge 5 74 {sharp}
Edge 118 359 {sharp}
Edge 684 850 {sharp}
Edge 684 733 {sharp}
Edge 233 401 {sharp}
Edge 233 235 {sharp}
Edge 234 236 {sharp}
Edge 234 235 {sharp}
Edge 8 339 {sharp}
Edge 8 267 {sharp}
Edge 8 20 {sharp}
Edge 236 406 {sharp}
Edge 10 43 {sharp}
Edge 14 67 {sharp}
Edge 14 33 {sharp}
Edge 129 374 {sharp}
Edge 469 637 {sharp}
Edge 132 219 {sharp}
Edge 698 742 {sharp}
Edge 359 482 {sharp}
Edge 20 24 {sharp}
Edge 20 160 {sharp}
Edge 21 58 {sharp}
Edge 24 113 {sharp}
Edge 139 180 {sharp}
Edge 140 363 {sharp}
Edge 140 574 {sharp}
Edge 367 508 {sharp}
Edge 255 367 {sharp}
Edge 482 698 {sharp}
Edge 258 375 {sharp}
Edge 33 132 {sharp}
Edge 260 501 {sharp}
Edge 147 213 {sharp}
Edge 147 178 {sharp}
Edge 374 560 {sharp}
Edge 375 495 {sharp}
Edge 149 388 {sharp}
Edge 151 195 {sharp}
Edge 718 765 {sharp}
Edge 718 813 {sharp}
Edge 379 425 {sharp}
Edge 267 665 {sharp}
Edge 495 670 {sharp}
Edge 43 85 {sharp}
Edge 611 706 {sharp}
Edge 612 617 {sharp}
Edge 273 569 {sharp}
Edge 160 258 {sharp}
Edge 160 267 {sharp}
Edge 47 48 {sharp}
Edge 47 64 {sharp}
Edge 48 87 {sharp}
Edge 161 255 {sharp}
Edge 161 209 {sharp}
Edge 501 612 {sharp}
Edge 164 311 {sharp}
Edge 54 99 {sharp}
Edge 621 668 {sharp}
Edge 508 566 {sharp}
Edge 58 110 {sharp}
Edge 627 717 {sharp}
Edge 401 524 {sharp}
Edge 64 139 {sharp}
Edge 178 273 {sharp}
Edge 406 469 {sharp}
Edge 293 413 {sharp}
Edge 67 68 {sharp}
Edge 180 331 {sharp}
Edge 68 203 {sharp}
Edge 634 691 {sharp}
Edge 637 683 {sharp}
Edge 524 627 {sharp}
Edge 413 479 {sharp}
Edge 74 133 {sharp}
Edge 533 733 {sharp}
Edge 533 535 {sharp}
Edge 194 231 {sharp}
Edge 195 293 {sharp}
Edge 535 691 {sharp}
Edge 85 194 {sharp}
Edge 765 850 {sharp}
Edge 87 151 {sharp}
Edge 314 331 {sharp}
Edge 314 479 {sharp}
Edge 88 363 {sharp}
Edge 88 118 {sharp}
Edge 203 425 {sharp}
Edge 94 164 {sharp}
Edge 434 574 {sharp}
Edge 434 558 {sharp}
Edge 209 260 {sharp}
Edge 99 149 {sharp}
Edge 665 670 {sharp}
Edge 439 611 {sharp}
Edge 326 388 {sharp}
Edge 213 439 {sharp}
Edge 668 706 {sharp}
Edge 558 742 {sharp}
Edge 219 379 {sharp}
Edge 560 618 {sharp}
Edge 110 231 {sharp}
