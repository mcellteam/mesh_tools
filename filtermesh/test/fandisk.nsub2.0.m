# Created using:
#  Subdivfit -fi - -mf - -record fandisk.nsub2.crep1e-4p.rec -crep 1e-4
#   -csharp .2e-4 -reconstruct
# 16475 points read
# Initial mesh:
# Genus: c=1 b=0  v=184 f=364 e=546  genus=0  sharpe=140 cuspv=0
# Starting reconstruction sequence
#  crep=0.0001, csharp=2e-05
#  internal xform: F 1  0.152204 0 0  0 0.152204 0  0 0 0.152204  0.1326 -1.81818 0.508021  0
# 
# (Timing on cpu=alpha host=wally.cs.washington.edu)
# fgfit_before: v=184 nse=140/546  edis=1.07804 etot=1.99317
# fgfit_after : v=184 nse=140/546  edis=0.149265 etot=1.0644
#  (_fgfit:                807.00)
# 
# Stoc, crep=0.0001 csharp=2e-05 wcrep=0.00431666 wcsharp=0.000863331
# stoc_before: v=184 nse=140/546  edis=0.149255 etot=1.06439
# it 682, last search: 81 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=89 f=174 e=261  genus=0  sharpe=92 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts       682         0         0         0
#              success        95         0         0         0
#      positive_energy       166         0         0         0
#         bad_dihedral        93         0         0         0
#            bad_sharp       299         0         0         0
#         illegal_move        29         0         0         0
#  last cedis=0.165531 cetot=0.62914
# stoc_after : v=89 nse=92/261  edis=0.165376 etot=0.628984
#  (_stoc:                2776.87)
# 
# fgfit_before: v=89 nse=92/261  edis=0.165373 etot=0.628982
# fgfit_after : v=89 nse=92/261  edis=0.150403 etot=0.614012
#  (_fgfit:                223.58)
# 
# Stoc, crep=0.0001 csharp=2e-05 wcrep=0.00431666 wcsharp=0.000863331
# stoc_before: v=89 nse=92/261  edis=0.150328 etot=0.613937
# it 348, last search: 46 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=89 f=174 e=261  genus=0  sharpe=97 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts         0         0         0       348
#              success         0         0         0         9
#      positive_energy         0         0         0       147
#         bad_dihedral         0         0         0         0
#            bad_sharp         0         0         0       192
#         illegal_move         0         0         0         0
#  last cedis=0.100298 cetot=0.568223
# stoc_after : v=89 nse=97/261  edis=0.100076 etot=0.568001
#  (_stoc:                 818.92)
# 
# fgfit_before: v=89 nse=97/261  edis=0.100077 etot=0.568002
# fgfit_after : v=89 nse=97/261  edis=0.0917022 etot=0.559628
#  (_fgfit:                211.27)
# 
# Stoc, crep=0.0001 csharp=2e-05 wcrep=0.00431666 wcsharp=0.000863331
# stoc_before: v=89 nse=97/261  edis=0.0917024 etot=0.559628
# it 310, last search: 22 wasted attempts
# New mesh:
# Genus: c=1 b=0  v=87 f=170 e=255  genus=0  sharpe=96 cuspv=0
# Summary of attempts and results:
#                           ecol      espl      eswa      esha
#       total_attempts       310         0       308       308
#              success         2         0        14         0
#      positive_energy       168         0       179       132
#         bad_dihedral       126         0        99         0
#            bad_sharp         0         0         0       176
#         illegal_move        14         0        16         0
#  last cedis=0.0761227 cetot=0.534551
# stoc_after : v=87 nse=96/255  edis=0.0760805 etot=0.534509
#  (_stoc:                2934.52)
# 
# fgfit_before: v=87 nse=96/255  edis=0.0760794 etot=0.534508
# fgfit_after : v=87 nse=96/255  edis=0.0711869 etot=0.529616
#  (_fgfit:                208.67)
#  (Subdivfit:            7982.07)
# Summary of timers (cpu=alpha host=wally.cs.washington.edu):
#  ___submesh:         (1141  )    0.23:3.15     av=     0.61   sum=   695.53
#  ____gmakeSpatial:   (10    )    2.63:4.30     av=     3.00   sum=    30.00
#  ____gspatialproject:(10    )    2.57:3.45     av=     2.76   sum=    27.62
#  ___gallproject:     (10    )    5.22:7.73     av=     5.78   sum=    57.85
#  ___computegrad:     (270   )    1.28:1.50     av=     1.35   sum=   365.70
#  ___gneighproject:   (500   )    1.93:2.37     av=     2.07   sum=  1035.77
#  __fgfit_iter:       (70    )   14.75:26.67    av=    20.27   sum=  1418.62
#  _fgfit:             (4     )  208.67:807.00   av=   362.63   sum=  1450.52
#  _initial_fit:       (3     )    6.57:9.97     av=     7.71   sum=    23.13
#  ____lmakespatial:   (1131  )    0.00:0.40     av=     0.09   sum=    98.37
#  ____lspatialproject:(1131  )    0.02:1.98     av=     0.17   sum=   197.38
#  ___lallproject:     (1131  )    0.02:2.37     av=     0.26   sum=   296.52
#  ____lneighproject:  (4954  )    0.15:3.85     av=     0.75   sum=  3725.02
#  ____lcombinations:  (4954  )    0.02:0.23     av=     0.09   sum=   446.03
#  _____qrdLLS:        (4954  )    0.02:0.57     av=     0.14   sum=   689.32
#  ____lsolve:         (4954  )    0.02:0.63     av=     0.17   sum=   822.83
#  ___loptimize:       (4954  )    0.20:4.52     av=     1.06   sum=  5261.90
#  __tryecol:          (650   )    0.88:25.27    av=     6.88   sum=  4472.07
#  __lattempt:         (1340  )    0.00:28.50    av=     4.84   sum=  6488.33
#  _stoc:              (3     )  818.92:2934.52  av=  2176.77   sum=  6530.30
#  __tryesha:          (288   )    2.08:16.23    av=     4.95   sum=  1424.22
#  __tryeswa:          (193   )    1.25:14.32    av=     3.06   sum=   590.77
#  Subdivfit:          (1     )        :         av=  7982.07   sum=  7982.07
# Summary of statistics:
# Sospobcells:(125696 )          1:25          av=4.4472      sd=2.16926
# Sssncellsv: (761648 )          1:729         av=6.4875      sd=23.8343
# Sssnelemsv: (761648 )          1:94          av=8.47552     sd=7.09962
# Sospcelln:  (114778 )          1:38          av=4.87023     sd=3.42588
# Sprojquick: (38032596)          0:1           av=0.465139    sd=0.498783
# Sprojnei:   (20342150)          0:7           av=0.0554936   sd=0.24906
# Sprojf:     (20342150)          4:61          av=11.144      sd=2.50591
# Sprojunexp: (1049950)          0:1           av=0.0034554   sd=0.0586811
# Sneval:     (70     )          6:9           av=7.14286     sd=0.571946
# Stmin:      (70     )  0.0018147:0.0144793   av=0.0049611   sd=0.00263704
# Ssetfrem:   (1131   )        288:1840        av=788.216     sd=191.062
# Ssetvrem:   (1131   )        152:957         av=401.437     sd=95.7759
# Secolpts:   (650    )       1412:13103       av=6733.58     sd=2214.68
# Secolmf:    (650    )         30:156         av=72.9646     sd=18.6734
# Secolmv:    (650    )          4:16          av=9.24308     sd=2.18962
# Secolsmv:   (650    )        265:1352        av=629.438     sd=160.692
# Soptnit:    (925    )          2:10          av=4.00757     sd=2.61136
# Soptnig:    (133    )         11:12          av=11.0226     sd=0.149046
# Sechange:   (120    ) -0.0144287:-0.000101271 av=-0.00421697 sd=0.00285219
# Seshapts:   (288    )       3532:11745       av=7130.12     sd=1692.38
# Seshamf:    (288    )         37:115         av=69.9167     sd=16.469
# Seshamv:    (288    )          6:15          av=9.75347     sd=1.93581
# Seshasmv:   (288    )        323:979         av=601.108     sd=140.699
# Seswapts:   (193    )       2053:8903        av=4942.85     sd=1543.44
# Seswamf:    (193    )         25:79          av=47.3005     sd=10.0738
# Seswamv:    (193    )          4:4           av=4           sd=0
# Seswasmv:   (193    )        215:706         av=414.705     sd=89.5502
Vertex 1  -3.00848e-05 15.4037 -1.38369
Vertex 8  0.414365 15.3993 -1.04481
Vertex 22  -5.76485e-05 15.4317 3.89295e-05
Vertex 47  0.48471 15.4579 -1.17114
Vertex 108  0.725949 15.4901 -3.10427e-05
Vertex 133  1.33479 15.7402 -0.809726
Vertex 160  1.38141 15.8522 -4.89972e-06
Vertex 209  1.97417 16.3312 -0.762236
Vertex 230  1.75725 16.1209 2.95574e-06
Vertex 278  3.50375 15.5058 -1.03916
Vertex 279  3.1716 15.3929 -1.32318
Vertex 284  2.72381 15.3029 -2.07461
Vertex 285  3.11489 15.7474 -0.834588
Vertex 306  2.81556 15.5032 -1.10869
Vertex 309  1.57133 15.3405 -1.2774
Vertex 473  2.70026 15.2003 -2.68019
Vertex 567  2.30435 12.6165 -2.2247
Vertex 569  1.40771 12.5753 -2.21747
Vertex 578  2.43993 12.753 -2.24872
Vertex 587  1.30866 12.7697 -2.2516
Vertex 591  3.1807 14.8952 -2.62684
Vertex 592  3.1832 14.0828 -2.4832
Vertex 614  2.69498 14.8951 -2.62633
Vertex 618  2.69844 14.0828 -2.48305
Vertex 623  2.42596 13.203 -2.32809
Vertex 646  2.69824 13.2365 -2.33407
Vertex 756  1.34739 13.5442 -2.38823
Vertex 849  1.27108 13.8752 -2.44661
Vertex 923  1.00863 14.2609 -2.51459
Vertex 1009  0.524061 14.4872 -2.5544
Vertex 1041  -8.00072e-05 14.4534 -2.54821
Vertex 1078  0.000124945 15.2037 -2.68079
Vertex 1138  -9.38589e-05 14.6327 -1.52452
Vertex 1199  3.76884e-05 14.6914 -0.000118889
Vertex 1319  4.82875 17.8329 -0.538184
Vertex 1321  4.82799 14.833 0.000165101
Vertex 1323  4.82769 17.8595 0.000216875
Vertex 1325  4.82769 15.6409 -0.169318
Vertex 1427  4.82831 15.5737 -0.58079
Vertex 1434  4.82788 15.2665 -0.521349
Vertex 1445  4.82783 15.1812 -1.00426
Vertex 1489  4.82817 14.37 -0.853046
Vertex 1499  4.82794 14.4521 -0.388764
Vertex 1540  4.82772 13.6268 0.000123933
Vertex 1583  4.82791 13.6009 -0.236819
Vertex 1593  4.19392 13.6133 -0.27746
Vertex 1603  3.12234 13.5326 6.3798e-05
Vertex 1626  3.27962 13.5101 -0.776473
Vertex 1664  2.42571 13.3952 -6.59291e-05
Vertex 1713  2.42877 13.3885 -0.985046
Vertex 1770  2.74662 13.366 -1.6422
Vertex 1856  2.7411 15.0078 -1.98815
Vertex 1865  3.18518 15.152 -1.17066
Vertex 1898  3.24334 15.0207 -1.91521
Vertex 1910  4.01444 15.2527 -0.5988
Vertex 1935  4.08795 15.5505 -0.63638
Vertex 1971  1.46173 12.8172 0.00111366
Vertex 1980  2.1851 12.7928 -0.000587371
Vertex 1982  1.28428 12.9534 -0.000244333
Vertex 2007  0.436298 14.6957 5.50683e-05
Vertex 2021  1.37894 13.986 0.000259938
Vertex 2059  0.780883 14.6108 -5.47941e-05
Vertex 2075  1.04872 14.3971 -3.13156e-05
Vertex 2145  2.4509 12.8393 0.000295352
Vertex 2535  2.32866 16.8368 0.000319488
Vertex 2777  3.95322 17.7941 0.000467548
Vertex 2983  2.95695 17.3642 -0.000271925
Vertex 3611  2.32891 16.8549 -0.642398
Vertex 3683  3.522 16.6535 -0.690141
Vertex 3701  3.15272 17.4632 -0.551551
Vertex 3778  3.80346 17.6941 -0.508094
Vertex 3817  4.4476 17.8449 -0.464979
Vertex 3821  4.25901 17.5565 -0.519747
Vertex 3879  4.43327 15.6459 -0.155032
Vertex 4035  3.4938 15.5961 -0.455572
Vertex 4193  2.9297 15.4779 -0.990426
Vertex 4267  2.42732 12.7528 -1.02137
Vertex 4287  1.46717 12.7909 -0.900875
Vertex 4950  0.585826 14.661 -1.48515
Vertex 5000  1.08078 14.3708 -1.26053
Vertex 5930  2.74559 14.1985 -1.8264
Vertex 5944  3.19883 14.3395 -1.02877
Vertex 5948  3.89261 14.4198 -0.571407
Vertex 5957  4.50141 14.4509 -0.396016
Vertex 6285  3.93052 14.3467 -0.986008
Vertex 6384  3.23711 14.201 -1.81212
Vertex 6503  4.01715 15.167 -1.08536
Face 6524  133 108 160
Face 6537  230 133 160
Face 6600  47 1 8
Face 6672  230 209 133
Face 6921  47 8 133
Face 6995  284 309 279
Face 7133  47 309 1
Face 7313  133 285 306
Face 8911  1499 1321 1540
Face 9301  1770 623 646
Face 9355  1540 1593 1583
Face 9427  1626 1603 1713
Face 9468  1713 1603 1664
Face 9504  1713 1770 1626
Face 9616  614 284 1856
Face 9634  1865 284 279
Face 9706  278 1865 279
Face 9741  1910 1427 1434
Face 9792  284 1865 1856
Face 10361  2021 1664 2075
Face 10425  1603 1321 1664
Face 10778  1540 1321 1603
Face 11361  1664 1321 2535
Face 14252  4267 1980 4287
Face 14329  569 567 4287
Face 14737  4287 1980 1971
Face 14744  1971 1982 4287
Face 15226  923 5000 4950
Face 15744  756 2021 849
Face 18179  1445 1489 6503
Face 18858  278 4193 4035
Face 18967  1323 3817 2777
Face 19091  3817 3778 2777
Face 19153  2777 3778 3701
Face 19215  3701 2983 2777
Face 19254  3611 2983 3701
Face 19284  3611 2535 2983
Face 20061  4193 306 285
Face 20098  2007 4950 2059
Face 20127  5948 1593 1626
Face 20397  1078 1009 1041
Face 21110  306 4193 278
Face 21481  1980 4267 2145
Face 21509  1325 3821 1319
Face 21516  1910 1935 1427
Face 21647  646 5930 1770
Face 21845  5948 5957 1593
Face 22159  1319 3821 3817
Face 22257  6503 6285 1898
Face 22258  6384 1898 6285
Face 22308  1489 6285 6503
Face 22419  1770 1713 623
Face 22581  4035 4193 285
Face 22584  1499 1583 1593
Face 22669  1593 5957 1499
Face 22716  2535 230 1664
Face 22734  3817 1323 1319
Face 22742  1593 1603 1626
Face 22843  5000 2075 2059
Face 22851  6384 591 1898
Face 22875  646 618 5930
Face 22899  5944 1626 5930
Face 22900  1770 5930 1626
Face 22906  2145 4267 1664
Face 22918  2535 3611 230
Face 22960  4287 587 569
Face 22985  2075 1664 2059
Face 23015  3611 209 230
Face 23044  3683 4035 285
Face 23046  1427 1325 1434
Face 23051  1971 1980 1982
Face 23058  278 279 306
Face 23079  1319 1323 1325
Face 23080  1321 1325 1323
Face 23094  3821 1325 3879
Face 23115  1499 1540 1583
Face 23125  2021 5000 849
Face 23126  923 849 5000
Face 23131  1540 1603 1593
Face 23156  6384 592 591
Face 23158  923 4950 1009
Face 23167  592 618 591
Face 23173  473 284 614
Face 23193  1982 1980 2021
Face 23222  1982 587 4287
Face 23226  1434 1499 1445
Face 23227  1489 1445 1499
Face 23237  133 209 285
Face 23238  5948 1626 5944
Face 23243  1009 1078 923
Face 23244  567 569 578
Face 23256  646 623 618
Face 23312  2021 2075 5000
Face 23313  1713 1664 623
Face 23321  6503 1898 1910
Face 23326  578 4267 567
Face 23327  4287 567 4267
Face 23345  108 2007 160
Face 23368  1138 1199 1
Face 23372  1664 230 2059
Face 23376  1321 1323 2535
Face 23397  1321 1499 1325
Face 23398  1434 1325 1499
Face 23419  309 284 1
Face 23536  587 1982 756
Face 23549  1856 1865 614
Face 23643  1434 1445 1910
Face 23644  6503 1910 1445
Face 23652  4950 5000 2059
Face 23654  569 587 578
Face 23656  3683 3821 4035
Face 23657  3879 4035 3821
Face 23696  587 756 578
Face 23697  623 578 756
Face 23698  8 108 133
Face 23700  3879 1325 4035
Face 23748  849 923 614
Face 23755  1078 473 923
Face 23761  473 614 923
Face 23810  618 592 5930
Face 23812  592 6384 5930
Face 23856  5944 5930 6384
Face 23912  5944 6384 5948
Face 23915  6384 6285 5948
Face 23921  6285 5957 5948
Face 23926  1 1199 22
Face 23961  2777 2983 1323
Face 23962  2535 1323 2983
Face 23997  1489 1499 6285
Face 24009  1499 5957 6285
Face 24011  278 4035 1935
Face 24012  1935 4035 1427
Face 24017  1910 1898 1865
Face 24019  1325 1427 4035
Face 24026  591 614 1898
Face 24028  623 756 618
Face 24037  1898 614 1865
Face 24071  849 614 756
Face 24075  591 618 614
Face 24076  756 614 618
Face 24078  3821 3683 3778
Face 24081  3701 3778 3611
Face 24097  3778 3683 3611
Face 24098  3817 3821 3778
Face 24166  1664 2021 2145
Face 24176  1982 2021 756
Face 24187  2145 2021 1980
Face 24189  279 309 306
Face 24190  230 160 2059
Face 24191  2007 2059 160
Face 24194  1138 4950 1199
Face 24195  2007 1199 4950
Face 24196  285 209 3683
Face 24197  3611 3683 209
Face 24198  1 284 1078
Face 24199  473 1078 284
Face 24200  1009 4950 1041
Face 24201  1138 1041 4950
Face 24202  8 1 108
Face 24203  22 108 1
Face 24204  47 133 309
Face 24205  306 309 133
Face 24206  278 1935 1865
Face 24207  1910 1865 1935
Face 24210  1078 1041 1
Face 24211  1138 1 1041
Face 24212  623 1664 578
Face 24213  4267 578 1664
Face 24214  108 22 2007
Face 24215  1199 2007 22
Edge 1583 1593 {sharp}
Edge 1 1078 {sharp}
Edge 1 22 {sharp}
Edge 567 4287 {sharp}
Edge 567 578 {sharp}
Edge 567 569 {sharp}
Edge 569 587 {sharp}
Edge 230 2535 {sharp}
Edge 1138 1199 {sharp}
Edge 1593 1626 {sharp}
Edge 578 623 {sharp}
Edge 923 1009 {sharp}
Edge 133 209 {sharp}
Edge 1489 1499 {sharp}
Edge 1489 6285 {sharp}
Edge 1603 1664 {sharp}
Edge 473 1078 {sharp}
Edge 473 614 {sharp}
Edge 22 1199 {sharp}
Edge 22 108 {sharp}
Edge 587 756 {sharp}
Edge 1041 1138 {sharp}
Edge 1041 1078 {sharp}
Edge 2059 2075 {sharp}
Edge 591 614 {sharp}
Edge 591 1898 {sharp}
Edge 591 592 {sharp}
Edge 592 6384 {sharp}
Edge 592 618 {sharp}
Edge 1499 5957 {sharp}
Edge 1499 1583 {sharp}
Edge 3879 4035 {sharp}
Edge 1626 1770 {sharp}
Edge 160 230 {sharp}
Edge 47 133 {sharp}
Edge 1856 1865 {sharp}
Edge 3778 3817 {sharp}
Edge 2535 2983 {sharp}
Edge 614 1856 {sharp}
Edge 1971 4287 {sharp}
Edge 1971 1982 {sharp}
Edge 1971 1980 {sharp}
Edge 278 1935 {sharp}
Edge 278 279 {sharp}
Edge 618 646 {sharp}
Edge 618 5930 {sharp}
Edge 279 309 {sharp}
Edge 279 284 {sharp}
Edge 5930 5944 {sharp}
Edge 1865 1910 {sharp}
Edge 849 923 {sharp}
Edge 284 473 {sharp}
Edge 623 646 {sharp}
Edge 623 1713 {sharp}
Edge 1980 2145 {sharp}
Edge 285 306 {sharp}
Edge 285 3683 {sharp}
Edge 1982 2021 {sharp}
Edge 2777 2983 {sharp}
Edge 3683 3821 {sharp}
Edge 3683 3778 {sharp}
Edge 5944 5948 {sharp}
Edge 1199 2007 {sharp}
Edge 6285 6384 {sharp}
Edge 1540 1603 {sharp}
Edge 1540 1583 {sharp}
Edge 1427 1434 {sharp}
Edge 1427 1935 {sharp}
Edge 5948 5957 {sharp}
Edge 1319 3821 {sharp}
Edge 1319 3817 {sharp}
Edge 1319 1323 {sharp}
Edge 1319 1325 {sharp}
Edge 756 849 {sharp}
Edge 1321 1323 {sharp}
Edge 1321 1540 {sharp}
Edge 1434 1910 {sharp}
Edge 1434 1445 {sharp}
Edge 4035 4193 {sharp}
Edge 1323 2777 {sharp}
Edge 306 4193 {sharp}
Edge 646 1770 {sharp}
Edge 1664 2145 {sharp}
Edge 1664 1713 {sharp}
Edge 1325 3879 {sharp}
Edge 1325 1427 {sharp}
Edge 3701 3778 {sharp}
Edge 2007 2059 {sharp}
Edge 1445 6503 {sharp}
Edge 1445 1489 {sharp}
Edge 1898 6503 {sharp}
Edge 209 3611 {sharp}
Edge 2021 2075 {sharp}
Edge 1009 1041 {sharp}
Edge 3611 3701 {sharp}
Edge 108 160 {sharp}
