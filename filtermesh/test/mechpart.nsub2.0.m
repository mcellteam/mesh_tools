# Created using:
#  Subdivfit -mf - -fi mechpart.4102.pts -crep 1e-5 -csharp 2e-06 -recon
# Initial mesh:
# Genus: c=1 b=0  v=184 f=376 e=564  genus=3  sharpe=147 cuspv=0
# 4102 points read
# Starting reconstruction sequence
#  crep=1e-05, csharp=2e-06
#  internal xform: F 1  1 0 0  0 1 0  0 0 1  0 0 0  0
# 
# (Timing on cpu=? host=ficus.salk.edu)
# fgfit_before: v=184 nse=147/564  edis=0.00763097 etot=0.00976497
# fgfit_after : v=184 nse=147/564  edis=0.000269079 etot=0.00240308
#  (_fgfit:                 61.89)
# 
# Stoc, crep=1e-05 csharp=2e-06 wcrep=1e-05 wcsharp=2e-06
# stoc_before: v=184 nse=147/564  edis=0.000269079 etot=0.00240308
# it 687, last search: 23 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=131 f=270 e=405  genus=3  sharpe=130 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts       687         0         0         0
#              success        53         0         0         0
#      positive_energy       195         0         0         0
#         bad_dihedral        81         0         0         0
#            bad_sharp       335         0         0         0
#         illegal_move        23         0         0         0
#  last cedis=0.000388642 cetot=0.00195864
# stoc_after : v=131 nse=130/405  edis=0.000386101 etot=0.0019561
#  (_stoc:                 296.38)
# 
# fgfit_before: v=131 nse=130/405  edis=0.000386101 etot=0.0019561
# fgfit_after : v=131 nse=130/405  edis=0.000317293 etot=0.00188729
#  (_fgfit:                 17.79)
# 
# Stoc, crep=1e-05 csharp=2e-06 wcrep=1e-05 wcsharp=2e-06
# stoc_before: v=131 nse=130/405  edis=0.00031729 etot=0.00188729
# it 623, last search: 74 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=131 f=270 e=405  genus=3  sharpe=124 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts         0         0         0       623
#              success         0         0         0        26
#      positive_energy         0         0         0       275
#         bad_dihedral         0         0         0         0
#            bad_sharp         0         0         0       322
#         illegal_move         0         0         0         0
#  last cedis=0.000219958 cetot=0.00177796
# stoc_after : v=131 nse=124/405  edis=0.000218452 etot=0.00177645
#  (_stoc:                 217.77)
# 
# fgfit_before: v=131 nse=124/405  edis=0.000218455 etot=0.00177645
# fgfit_after : v=131 nse=124/405  edis=0.000195273 etot=0.00175327
#  (_fgfit:                 17.52)
# 
# Stoc, crep=1e-05 csharp=2e-06 wcrep=1e-05 wcsharp=2e-06
# stoc_before: v=131 nse=124/405  edis=0.00019527 etot=0.00175327
# it 480, last search: 7 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=125 f=258 e=387  genus=3  sharpe=120 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts       480         0       474       474
#              success         6         0        18         0
#      positive_energy       232         0       253       208
#         bad_dihedral       230         0       191         0
#            bad_sharp         0         0         0       266
#         illegal_move        12         0        12         0
#  last cedis=0.000175419 cetot=0.00166542
# stoc_after : v=125 nse=120/387  edis=0.000176277 etot=0.00166628
#  (_stoc:                 634.14)
# 
# fgfit_before: v=125 nse=120/387  edis=0.000176278 etot=0.00166628
# fgfit_after : v=125 nse=120/387  edis=0.000154696 etot=0.0016447
#  (_fgfit:                 18.27)
#  (Subdivfit:            1263.95)
# Summary of timers (cpu=? host=ficus.salk.edu):
#  ___submesh:         (1587  )    0.09:0.77     av=     0.25   sum=   402.49
#  ____gmakeSpatial:   (10    )    0.83:1.05     av=     0.88   sum=     8.79
#  ____gspatialproject:(10    )    0.35:0.44     av=     0.37   sum=     3.73
#  ___gallproject:     (10    )    1.20:1.48     av=     1.27   sum=    12.66
#  ___computegrad:     (267   )    0.07:0.09     av=     0.07   sum=    19.48
#  ___gneighproject:   (477   )    0.15:0.19     av=     0.17   sum=    80.15
#  __fgfit_iter:       (70    )    1.21:1.97     av=     1.53   sum=   107.28
#  _fgfit:             (4     )   17.52:61.89    av=    28.87   sum=   115.47
#  _initial_fit:       (3     )    1.71:2.05     av=     1.83   sum=     5.48
#  __lattempt:         (1790  )   -0.01:3.05     av=     0.64   sum=  1137.62
#  ____lmakespatial:   (1577  )    0.00:0.09     av=     0.02   sum=    29.84
#  ____lspatialproject:(1577  )    0.00:0.17     av=     0.02   sum=    24.83
#  ___lallproject:     (1577  )    0.00:0.25     av=     0.04   sum=    55.21
#  ____lneighproject:  (6333  )    0.01:0.14     av=     0.05   sum=   345.75
#  ____lcombinations:  (6333  )    0.00:0.03     av=     0.01   sum=    51.89
#  _____qrdLLS:        (6333  )   -0.01:0.16     av=     0.02   sum=   110.02
#  ____lsolve:         (6333  )    0.00:0.17     av=     0.02   sum=   128.16
#  ___loptimize:       (6333  )    0.02:0.37     av=     0.09   sum=   563.59
#  __tryecol:          (797   )    0.29:2.84     av=     0.83   sum=   661.30
#  _stoc:              (3     )  217.77:634.14   av=   382.76   sum=  1148.28
#  __tryesha:          (509   )    0.34:1.70     av=     0.71   sum=   362.93
#  __tryeswa:          (271   )    0.20:1.27     av=     0.42   sum=   112.56
#  Subdivfit:          (1     )        :         av=  1263.95   sum=  1263.95
# Summary of statistics:
# Sospobcells:(170176 )          1:13          av=3.08913     sd=1.40586
# Sssncellsv: (182435 )          1:216         av=4.80585     sd=12.4308
# Sssnelemsv: (182435 )          1:153         av=14.0912     sd=11.6012
# Sospcelln:  (73408  )          1:53          av=7.16128     sd=5.13504
# Sprojquick: (9693596)          0:1           av=0.46084     sd=0.498464
# Sprojnei:   (5226400)          0:8           av=0.0576647   sd=0.252145
# Sprojf:     (5226400)          4:60          av=11.3407     sd=2.81197
# Sprojunexp: (282198 )          0:1           av=0.00464567  sd=0.0680008
# Sneval:     (70     )          6:8           av=6.81429     sd=0.665786
# Stmin:      (70     )  0.0245868:0.104848    av=0.0588161   sd=0.0147243
# Ssetfrem:   (1577   )        352:2032        av=1178.54     sd=246.387
# Ssetvrem:   (1577   )        176:978         av=595.358     sd=117.941
# Secolpts:   (797    )        472:2895        av=1430.62     sd=443.34
# Secolmf:    (797    )         36:186         av=95.0151     sd=28.2069
# Secolmv:    (797    )          5:21          av=10.0414     sd=2.85707
# Secolsmv:   (797    )        309:1584        av=821.764     sd=238.5
# Soptnit:    (1277   )          2:10          av=3.7502      sd=2.17158
# Soptnig:    (114    )         11:12          av=11.0789     sd=0.270847
# Sechange:   (103    )-2.92948e-05:-5.22439e-07 av=-6.22933e-06 sd=5.04782e-06
# Seshapts:   (509    )        554:2918        av=1347.95     sd=455.451
# Seshamf:    (509    )         45:192         av=96.4833     sd=27.2882
# Seshamv:    (509    )          6:22          av=11.0098     sd=3.11554
# Seshasmv:   (509    )        399:1625        av=833.527     sd=229.138
# Seswapts:   (271    )        360:1627        av=891.089     sd=299.667
# Seswamf:    (271    )         29:99          av=58.8044     sd=15.5059
# Seswamv:    (271    )          4:4           av=4           sd=0
# Seswasmv:   (271    )        259:869         av=518.358     sd=137.483
Vertex 1  0.366669 0.424629 0.434301
Vertex 2  0.434133 0.430277 0.439116
Vertex 7  0.366719 0.396055 0.495517
Vertex 8  0.366661 0.366664 0.366665
Vertex 11  0.366673 0.461861 0.366662
Vertex 14  0.453526 0.366609 0.407096
Vertex 18  0.408076 0.366718 0.460993
Vertex 23  0.366614 0.486528 0.395573
Vertex 27  0.465066 0.462913 0.411176
Vertex 42  0.366669 0.366676 0.633326
Vertex 47  0.447764 0.409826 0.366676
Vertex 48  0.405956 0.445566 0.366697
Vertex 52  0.464132 0.468198 0.407043
Vertex 53  0.402523 0.366668 0.50972
Vertex 58  0.366613 0.422718 0.55693
Vertex 60  0.366692 0.563666 0.422051
Vertex 63  0.394237 0.522352 0.366665
Vertex 65  0.396724 0.46739 0.286529
Vertex 68  0.497234 0.506841 0.400803
Vertex 69  0.527692 0.366732 0.402948
Vertex 72  0.469654 0.538475 0.410125
Vertex 74  0.413127 0.412216 0.562407
Vertex 84  0.366674 0.633315 0.366666
Vertex 94  0.560143 0.442493 0.426295
Vertex 97  0.527173 0.494141 0.401722
Vertex 100  0.44872 0.536355 0.422144
Vertex 102  0.425291 0.366599 0.564495
Vertex 103  0.483542 0.481287 0.59473
Vertex 107  0.496655 0.396132 0.366658
Vertex 110  0.366721 0.466031 0.592748
Vertex 112  0.36669 0.502953 0.633345
Vertex 113  0.366635 0.633395 0.436567
Vertex 118  0.403508 0.545306 0.100027
Vertex 124  0.569067 0.366589 0.433361
Vertex 136  0.416288 0.599287 0.466372
Vertex 147  0.366672 0.633436 0.439884
Vertex 149  0.402835 0.584539 0.468445
Vertex 150  0.36671 0.590628 0.458782
Vertex 151  0.4355 0.586695 0.366699
Vertex 154  0.393566 0.48503 0.0999492
Vertex 167  0.466682 0.633408 0.40421
Vertex 177  0.426608 0.63326 0.439079
Vertex 179  0.633337 0.366644 0.366665
Vertex 180  0.547098 0.407137 0.366677
Vertex 181  0.496166 0.394241 0.0999824
Vertex 191  0.366653 0.633326 0.463301
Vertex 193  0.407276 0.59246 0.515728
Vertex 194  0.366803 0.584152 0.547303
Vertex 204  0.591964 0.366807 0.465547
Vertex 206  0.633297 0.455718 0.404624
Vertex 212  0.554142 0.551762 0.417356
Vertex 213  0.531369 0.63332 0.407593
Vertex 219  0.473304 0.366745 0.596053
Vertex 225  0.579921 0.434715 0.366659
Vertex 231  0.366558 0.530341 0.5984
Vertex 234  0.429993 0.419636 0.633333
Vertex 235  0.402609 0.467469 0.63334
Vertex 236  0.491141 0.392263 0.633336
Vertex 255  0.633313 0.400861 0.466283
Vertex 258  0.633325 0.480174 0.366663
Vertex 259  0.633415 0.535199 0.406465
Vertex 267  0.490376 0.366668 0.633337
Vertex 269  0.535335 0.366647 0.593007
Vertex 285  0.3975 0.506819 0.633292
Vertex 293  0.489337 0.604625 0.366672
Vertex 294  0.549397 0.593315 0.366668
Vertex 305  0.428922 0.416439 0.100022
Vertex 309  0.633334 0.366676 0.633335
Vertex 311  0.597037 0.405871 0.49464
Vertex 316  0.633308 0.572769 0.4374
Vertex 319  0.600547 0.599427 0.469096
Vertex 320  0.633338 0.633327 0.36666
Vertex 322  0.525861 0.477775 0.595243
Vertex 330  0.397076 0.633362 0.53768
Vertex 331  0.605263 0.48436 0.366695
Vertex 333  0.550184 0.407934 0.100016
Vertex 339  0.366651 0.633339 0.633323
Vertex 343  0.41018 0.552811 0.633323
Vertex 365  0.593956 0.366654 0.505521
Vertex 367  0.633341 0.405891 0.534375
Vertex 368  0.580501 0.366657 0.559116
Vertex 376  0.569229 0.63331 0.434766
Vertex 386  0.445783 0.564803 0.586168
Vertex 406  0.561919 0.412428 0.633319
Vertex 410  0.425505 0.423824 0.900005
Vertex 429  0.585391 0.415051 0.558535
Vertex 449  0.46165 0.633344 0.592665
Vertex 462  0.633339 0.633347 0.633344
Vertex 468  0.400476 0.473035 0.899997
Vertex 474  0.485429 0.39339 0.900004
Vertex 477  0.593666 0.550933 0.366673
Vertex 493  0.633359 0.456077 0.59413
Vertex 503  0.598321 0.633363 0.471982
Vertex 523  0.455189 0.596369 0.633363
Vertex 524  0.418844 0.569452 0.9
Vertex 541  0.55616 0.586322 0.0999818
Vertex 558  0.605216 0.516824 0.162905
Vertex 560  0.633307 0.595357 0.484381
Vertex 566  0.633316 0.537031 0.59714
Vertex 568  0.519234 0.633285 0.599359
Vertex 573  0.589891 0.445648 0.0999922
Vertex 586  0.596254 0.45581 0.633331
Vertex 588  0.605367 0.480909 0.80394
Vertex 590  0.541727 0.405255 0.899999
Vertex 611  0.633323 0.633327 0.529848
Vertex 612  0.633387 0.593118 0.541985
Vertex 613  0.5882 0.633353 0.557812
Vertex 618  0.54011 0.544414 0.590338
Vertex 626  0.520674 0.604915 0.633197
Vertex 637  0.604298 0.520449 0.633302
Vertex 652  0.450093 0.593723 0.0999775
Vertex 655  0.605363 0.494625 0.0999983
Vertex 667  0.590515 0.589504 0.546463
Vertex 674  0.480634 0.603157 0.691095
Vertex 696  0.505937 0.605353 0.100024
Vertex 724  0.398081 0.516033 0.900004
Vertex 738  0.571158 0.5706 0.154334
Vertex 754  0.575416 0.572496 0.633381
Vertex 756  0.549015 0.594381 0.89997
Vertex 765  0.604608 0.507798 0.899975
Vertex 768  0.591419 0.550488 0.900138
Vertex 772  0.59071 0.441138 0.900028
Vertex 780  0.590235 0.554328 0.100011
Vertex 843  0.479466 0.608034 0.900019
Vertex 877  0.408064 0.410767 0.472194
Face 10  2 1 7
Face 11  7 1 8
Face 31  1 11 8
Face 35  1 23 11
Face 39  2 14 27
Face 62  7 8 42
Face 69  11 48 47
Face 79  2 27 52
Face 86  53 18 877
Face 98  23 60 11
Face 103  11 63 48
Face 105  47 48 65
Face 113  69 68 27
Face 125  72 23 52
Face 128  7 58 74
Face 145  11 47 8
Face 155  60 84 11
Face 162  48 63 65
Face 168  27 14 69
Face 175  68 69 94
Face 187  68 100 72
Face 193  53 74 102
Face 200  8 47 107
Face 207  42 112 58
Face 216  118 65 63
Face 228  14 8 69
Face 233  68 94 97
Face 243  42 53 102
Face 271  60 149 150
Face 272  60 150 113
Face 278  11 84 63
Face 279  151 63 84
Face 283  154 65 118
Face 315  68 167 100
Face 319  100 136 60
Face 342  177 149 136
Face 344  8 107 179
Face 345  180 179 107
Face 366  113 150 147
Face 367  191 147 150
Face 370  149 194 150
Face 405  212 68 97
Face 429  147 191 177
Face 438  225 179 180
Face 439  107 181 180
Face 452  112 234 235
Face 500  94 206 97
Face 501  259 97 206
Face 510  68 212 213
Face 516  84 177 167
Face 574  151 84 293
Face 580  63 151 118
Face 599  124 204 94
Face 608  206 94 255
Face 611  206 258 259
Face 615  97 259 212
Face 616  316 212 259
Face 636  110 231 103
Face 654  331 258 225
Face 712  255 179 258
Face 717  94 204 311
Face 718  365 311 204
Face 721  255 311 367
Face 765  231 386 103
Face 772  193 330 386
Face 773  194 193 386
Face 795  343 112 285
Face 810  236 410 234
Face 831  154 305 65
Face 839  311 365 368
Face 863  213 376 320
Face 885  231 194 386
Face 908  65 305 47
Face 926  234 410 235
Face 936  236 474 410
Face 970  493 309 367
Face 971  367 429 493
Face 975  316 258 320
Face 977  259 258 316
Face 993  319 503 212
Face 1001  429 368 322
Face 1031  330 449 386
Face 1034  339 462 449
Face 1045  47 181 107
Face 1047  47 305 181
Face 1081  294 477 541
Face 1113  429 322 493
Face 1139  212 503 376
Face 1179  343 524 523
Face 1199  406 309 586
Face 1207  474 236 590
Face 1257  331 558 477
Face 1266  376 503 320
Face 1268  560 612 319
Face 1276  322 566 493
Face 1284  103 568 618
Face 1285  613 618 568
Face 1300  462 523 626
Face 1328  236 406 590
Face 1352  652 118 151
Face 1375  331 225 573
Face 1376  331 573 558
Face 1392  449 462 568
Face 1403  523 524 674
Face 1404  523 674 626
Face 1424  406 586 588
Face 1487  613 667 618
Face 1542  151 293 696
Face 1543  294 696 293
Face 1544  151 696 652
Face 1558  566 612 611
Face 1566  573 655 558
Face 1594  674 756 626
Face 1603  343 724 524
Face 1617  754 637 462
Face 1627  772 406 588
Face 1650  612 566 667
Face 1655  738 541 477
Face 1672  626 756 754
Face 1692  772 765 590
Face 1694  637 754 768
Face 1714  294 541 696
Face 1721  738 780 541
Face 1722  655 780 558
Face 1773  780 738 558
Face 1776  738 477 558
Face 1804  772 590 406
Face 1809  586 637 588
Face 1810  768 588 637
Face 1839  588 768 765
Face 1840  588 765 772
Face 1849  674 843 756
Face 1875  765 768 756
Face 1924  331 477 258
Face 1938  754 756 768
Face 1956  367 311 429
Face 1957  368 429 311
Face 1968  319 212 316
Face 1974  136 100 177
Face 1975  167 177 100
Face 1986  311 255 94
Face 1991  319 612 667
Face 2002  8 179 69
Face 2003  124 69 179
Face 2010  84 320 293
Face 2011  294 293 320
Face 2020  177 191 330
Face 2030  322 368 269
Face 2032  611 462 566
Face 2036  84 113 177
Face 2037  147 177 113
Face 2053  258 477 320
Face 2057  42 58 7
Face 2058  219 269 267
Face 2062  255 258 206
Face 2066  560 611 612
Face 2072  560 319 316
Face 2073  191 339 330
Face 2074  449 330 339
Face 2077  7 877 2
Face 2079  7 74 877
Face 2080  53 877 74
Face 2081  58 112 110
Face 2090  285 112 235
Face 2091  568 462 613
Face 2094  462 611 613
Face 2108  258 179 225
Face 2114  368 365 309
Face 2128  136 149 60
Face 2131  112 42 234
Face 2137  154 118 305
Face 2142  204 124 179
Face 2145  18 8 14
Face 2150  477 294 320
Face 2155  103 322 219
Face 2156  269 219 322
Face 2164  765 524 410
Face 2165  493 566 309
Face 2166  462 309 566
Face 2168  754 462 626
Face 2171  23 72 60
Face 2172  100 60 72
Face 2173  212 376 213
Face 2175  320 611 316
Face 2176  560 316 611
Face 2189  204 179 365
Face 2190  309 365 179
Face 2191  474 590 410
Face 2192  765 410 590
Face 2193  112 343 339
Face 2199  84 60 113
Face 2201  269 368 267
Face 2202  309 267 368
Face 2203  181 333 180
Face 2205  18 53 8
Face 2206  42 8 53
Face 2213  42 267 234
Face 2214  236 234 267
Face 2215  674 524 843
Face 2218  194 231 339
Face 2219  42 102 267
Face 2220  219 267 102
Face 2224  667 566 618
Face 2227  69 124 94
Face 2229  193 194 149
Face 2235  235 468 285
Face 2238  468 235 410
Face 2239  343 285 724
Face 2240  468 724 285
Face 2251  27 68 52
Face 2252  72 52 68
Face 2256  213 167 68
Face 2259  468 410 724
Face 2260  524 724 410
Face 2262  181 305 118
Face 2263  333 118 573
Face 2264  655 573 118
Face 2267  333 181 118
Face 2270  655 118 780
Face 2272  780 118 541
Face 2273  652 696 118
Face 2274  541 118 696
Face 2275  267 309 236
Face 2276  406 236 309
Face 2277  343 523 339
Face 2278  462 339 523
Face 2281  167 213 84
Face 2282  320 84 213
Face 2283  309 462 586
Face 2284  637 586 462
Face 2285  177 330 149
Face 2286  193 149 330
Face 2287  2 52 1
Face 2288  23 1 52
Face 2289  191 150 339
Face 2290  194 339 150
Face 2291  179 255 309
Face 2292  367 309 255
Face 2293  319 667 503
Face 2294  613 503 667
Face 2295  103 386 568
Face 2296  449 568 386
Face 2297  103 618 322
Face 2298  566 322 618
Face 2299  503 613 320
Face 2300  611 320 613
Face 2301  14 2 18
Face 2302  877 18 2
Face 2303  225 180 573
Face 2304  333 573 180
Face 2305  58 110 74
Face 2306  103 74 110
Face 2307  112 339 110
Face 2308  231 110 339
Face 2309  102 74 219
Face 2310  103 219 74
Face 2311  524 765 843
Face 2312  756 843 765
Edge 339 462 {sharp}
Edge 113 147 {sharp}
Edge 566 612 {sharp}
Edge 1 23 {sharp}
Edge 1 7 {sharp}
Edge 568 613 {sharp}
Edge 343 523 {sharp}
Edge 118 652 {sharp}
Edge 118 154 {sharp}
Edge 7 58 {sharp}
Edge 573 655 {sharp}
Edge 8 179 {sharp}
Edge 8 11 {sharp}
Edge 8 42 {sharp}
Edge 234 236 {sharp}
Edge 234 235 {sharp}
Edge 235 285 {sharp}
Edge 462 611 {sharp}
Edge 236 406 {sharp}
Edge 11 84 {sharp}
Edge 124 204 {sharp}
Edge 14 69 {sharp}
Edge 14 18 {sharp}
Edge 468 724 {sharp}
Edge 18 53 {sharp}
Edge 586 637 {sharp}
Edge 474 590 {sharp}
Edge 136 149 {sharp}
Edge 23 60 {sharp}
Edge 590 772 {sharp}
Edge 365 368 {sharp}
Edge 367 493 {sharp}
Edge 255 367 {sharp}
Edge 258 320 {sharp}
Edge 259 316 {sharp}
Edge 147 191 {sharp}
Edge 149 177 {sharp}
Edge 376 503 {sharp}
Edge 150 194 {sharp}
Edge 151 293 {sharp}
Edge 493 566 {sharp}
Edge 267 309 {sharp}
Edge 154 305 {sharp}
Edge 42 267 {sharp}
Edge 42 112 {sharp}
Edge 269 368 {sharp}
Edge 47 107 {sharp}
Edge 47 48 {sharp}
Edge 48 63 {sharp}
Edge 503 613 {sharp}
Edge 618 667 {sharp}
Edge 53 102 {sharp}
Edge 167 213 {sharp}
Edge 167 177 {sharp}
Edge 58 110 {sharp}
Edge 285 343 {sharp}
Edge 60 150 {sharp}
Edge 626 754 {sharp}
Edge 63 151 {sharp}
Edge 177 330 {sharp}
Edge 179 309 {sharp}
Edge 179 258 {sharp}
Edge 406 586 {sharp}
Edge 293 294 {sharp}
Edge 180 225 {sharp}
Edge 294 477 {sharp}
Edge 68 212 {sharp}
Edge 68 94 {sharp}
Edge 181 305 {sharp}
Edge 181 333 {sharp}
Edge 69 124 {sharp}
Edge 523 626 {sharp}
Edge 410 468 {sharp}
Edge 410 474 {sharp}
Edge 637 754 {sharp}
Edge 524 724 {sharp}
Edge 524 843 {sharp}
Edge 74 103 {sharp}
Edge 74 877 {sharp}
Edge 756 843 {sharp}
Edge 756 768 {sharp}
Edge 191 339 {sharp}
Edge 193 194 {sharp}
Edge 193 386 {sharp}
Edge 194 231 {sharp}
Edge 309 462 {sharp}
Edge 84 113 {sharp}
Edge 84 320 {sharp}
Edge 311 429 {sharp}
Edge 311 365 {sharp}
Edge 765 768 {sharp}
Edge 765 772 {sharp}
Edge 652 696 {sharp}
Edge 541 780 {sharp}
Edge 541 696 {sharp}
Edge 655 780 {sharp}
Edge 316 560 {sharp}
Edge 204 365 {sharp}
Edge 319 667 {sharp}
Edge 206 255 {sharp}
Edge 206 259 {sharp}
Edge 320 611 {sharp}
Edge 322 429 {sharp}
Edge 212 319 {sharp}
Edge 100 136 {sharp}
Edge 213 376 {sharp}
Edge 102 219 {sharp}
Edge 103 618 {sharp}
Edge 103 386 {sharp}
Edge 103 322 {sharp}
Edge 330 449 {sharp}
Edge 331 477 {sharp}
Edge 219 269 {sharp}
Edge 333 573 {sharp}
Edge 107 180 {sharp}
Edge 560 612 {sharp}
Edge 449 568 {sharp}
Edge 110 231 {sharp}
Edge 112 339 {sharp}
Edge 225 331 {sharp}
