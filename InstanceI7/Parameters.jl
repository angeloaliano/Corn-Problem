using JuMP, Printf, LinearAlgebra, Random, Gaston, Gurobi, JLD2, FileIO

origem = 34
silos = 67
rodo = 77
ferro = 37
hidro = 30
portos = 11
dest_in = 26
dest_ex = 2
estados=26

T = 1:1
safra=2
O = 1:origem
D_i = origem+silos+rodo+ferro+hidro+portos+1:origem+silos+rodo+ferro+hidro+portos+dest_in
D_e = origem+silos+rodo+ferro+hidro+portos+dest_in+1:origem+silos+rodo+ferro+hidro+portos+dest_in+dest_ex
Silos = origem+1:origem+silos
V = origem+silos+1:origem+silos+rodo+ferro+hidro+portos+dest_in+dest_ex
n = origem + silos + rodo + ferro + hidro + portos + dest_in + dest_ex

L = [
-11.3	-62.6	1
-8.9	-70.2	2
-3.8	-64.4	3
1.4	    -61.2	4
-4.7	-53.1	5
-10.3	-48.1	6
-5.3	-45.3	7
-7.8	-42.6	8
-5.2	-39.4	9
-5.8	-36.4	10
-7.2	-36.6	11
-8.6	-37.9	12
-9.7	-36.5	13
-10.6	-37.4	14
-12.7	-41.7	15
-18.8	-44.5	16
-19.7	-40.6	17
-22.5	-43.5	18
-22.0	-48.9	19
-23.3	-53.1	20
-25.2	-53.5	20
-24.2	-50.4	20
-25.3	-49.9	20
-26.1	-51.5	20
-27.1	-50.1	21
-29.5	-53.3	22
-10.4	-58.9	23
-11.0	-52.5	23
-13.6	-58.5	23
-13.7	-53.5	23
-16.5	-55.4	23
-20.8	-54.5	24
-16.1	-49.5	25
-15.8	-47.7	26
-11.3	-61.4	1
-8.6	-64.1	1
-10.0	-67.8	2
-3.1	-60.0	3
2.8	    -60.7	4
-1.3	-48.4	5
-11.8	-49.5	6
-5.5	-47.4	7
-2.6	-44.2	7
-6.8	-43.0	8
-2.9	-41.7	8
-7.1	-41.5	8
-5.1	-42.8	8
-7.2	-39.3	9
-3.9    -38.6	9
-4.9	-38.0	9
-6.4	-38.9	9
-6.4	-39.3	9
-5.6	-39.4	9
-3.7	-40.4	9
-5.2	-40.7	9
-5.6	-36.9	10
-6.3	-36.5	10
-5.2	-37.4	10
-6.5	-37.1	10
-5.8	-35.2	10
-6.0	-37.8	10
-7.3	-35.9	11
-7.1	-34.9	11
-7.9	-37.1	11
-7.0	-37.3	11
-8.4	-37.1	12
-8.1	-34.9	12
-9.7	-35.7	13
-9.4	-36.7	13
-10.7	-37.5	14
-12.5	-40.3	15
-10.8	-38.6	15
-11.3	-41.8	15
-16.7	-43.9	16
-20.9	-47.0	16
-21.1	-45.1	16
-19.0	-48.3	16
-21.6	-45.5	16
-19.7	-47.9	16
-19.7	-46.2	16
-20.3	-40.3	17
-20.9	-41.1	17
-19.5	-40.6	17
-22.9	-43.2	18
-23.0	-49.5	19
-22.2	-49.7	19
-23.3	-51.3	20
-25.1	-50.2	20
-23.3	-51.4	20
-27.2	-51.5	21
-29.9	-51.2	22
-16.4	-54.6	23
-12.5	-55.7	23
-20.5	-54.6	24
-17.8	-50.9	25
-16.5	-50.4	25
-16.7	-49.2	25
-15.6	-50.0	25
-17.5	-49.4	25
-17.8	-50.2	25
-15.7	-48.0	26
-32.0	-52.1	22
-30.0	-51.2	22
-27.8	-50.3	21
-25.5	-48.5	20
-24.0	-46.3	19
-22.9	-43.2	18
-20.2	-42.1	16
-18.9	-42.0	16
-16.5	-41.5	16
-13.8	-40.1	15
-11.0	-38.8	15
-8.1	-39.1	12
-4.9	-38.0	9
-3.7	-38.5	9
-29.9	-50.3	22
-28.2	-48.7	21
-26.9	-48.7	21
-26.2	-48.6	21
-20.3	-40.3	17
-16.4	-39.6	15
-12.1	-38.4	15
-10.9	-37.1	14
-9.7	-35.7	13
-8.1	-34.9	12
-7.1	-34.9	11
-5.8	-35.2	10
-10.0	-67.8	2
-8.6	-64.1	1
-11.3	-61.4	1
-13.5	-57.9	23
-15.6	-56.1	23
-17.9	-51.7	25
-20.0	-48.9	16
-22.6	-47.4	19
-7.5	-63.0	3
-7.2	-59.9	3
-4.2	-56.0	5
-3.1	-52.1	5
-5.4	-49.1	5
-7.5	-46.0	7
-7.1	-41.5	8
-7.0	-37.3	11
-25.0	-53.5	20
-22.2	-54.8	24
-20.5	-54.6	24
-18.5	-54.7	24
-11.9	-55.5	23
-7.0	-55.4	5
-2.4	-54.7	5
-30.0	-52.9	22
-27.0	-51.7	22
-24.5	-50.4	20
-21.7	-49.7	19
-16.7	-49.2	25
-14.4	-49.2	25
-11.7	-49.1	6
-8.8	-48.5	6
-29.7	-53.8	22
-26.2	-52.7	20
-23.1	-52.5	20
-20.9	-51.4	19
-15.9	-52.3	23
-13.0	-51.6	23
-10.6	-51.5	23
-7.8	-49.9	5
2.8	    -60.7	4
1.1	    -60.4	4
-3.1	-60.0	3
-19.9	-43.9	16
-16.7	-43.9	16
-13.3	-44.4	15
-10.9	-45.2	15
-8.1	-43.7	8
-5.2	-44.5	7
-2.6	-44.4	7
-20.2	-56.4	24
-19.6	-47.0	16
-32.0	-52.1	22
-29.9	-54.8	22
-28.3	-54.3	22
-30.0	-51.2	22
-26.2	-48.6	21
-25.5	-48.5	20
-23.3	-51.2	20
-25.1	-50.2	20
-25.4	-51.5	20
-21.8	-52.1	19
-23.0	-49.9	19
-23.3	-47.8	19
-24.0	-46.3	19
-21.8	-48.2	19
-20.2	-51.0	19
-16.7	-49.2	25
-19.7	-47.9	16
-17.4	-44.9	16
-22.5	-44.5	18
-18.9	-42.0	16
-20.3	-40.3	17
-16.4	-54.6	23
-18.8	-52.6	24
-17.8	-50.9	25
-14.4	-49.2	25
-10.7	-48.4	6
-19.0	-57.7	24
-20.8	-51.7	24
-2.6	-44.4	7
-5.1	-42.8	8
-5.2	-40.7	9
-3.7	-40.4	9
-3.7	-38.5	9
-5.4	-49.1	5
-5.0	-47.5	7
-3.7	-45.4	7
-6.6	-47.4	7
-12.0	-48.5	6
-10.7	-48.4	6
-5.4	-49.1	5
-1.5	-48.7	5
-15.9	-52.3	23
-14.7	-52.4	23
-8.3	-49.3	5
-9.4	-40.5	15
-12.2	-43.2	15
-17.3	-44.9	16
-8.6	-64.1	1
-7.5	-63.0	3
-3.1	-58.4	3
-25.4	-54.6	20
-23.9	-54.3	24
-21.8	-52.1	19
-20.2	-51.0	19
-19.0	-50.5	25
-22.3	-48.8	19
-32.0	-52.1	22
-30.0	-51.2	22
-30.0	-52.9	22
-6.5	-64.4	3
-3.1	-60.0	3
-0.4	-65.0	3
-2.5	-66.1	3
-6.4	-68.2	3
2.8  	-60.7	4
-4.2	-56.0	5
-2.4	-54.7	5
-24.0	-46.3	19
-25.5	-48.5	20
-32.0	-52.1	22
-26.2	-48.6	21
-2.6	-44.4	7
-20.3	-40.3	17
-3.1	-58.4	3
-2.4	-54.7	5
-1.5	-48.7	5
-28.2	-48.7	21
-26.9	-48.7	21
-10.9	-62.0	1
-10.0	-67.8	2
-3.1	-60.0	3
2.8	   -60.7	4
-5.4	-49.1	5
-7.1	-48.0	6
-2.6	-44.2	7
-5.1	-42.8	8
-3.5	-39.7	9
-6.5	-37.1	10
-7.3	-35.9	11
-8.2	-36.0	12
-9.8	-36.7	13
-10.2	-37.4	14
-14.7	-40.8	15
-19.0	-48.3	16
-20.9	-41.1	17
-21.2	-41.9	18
-22.8	-47.1	19
-25.0	-53.5	20
-26.8	-49.2	21
-28.2	-52.5	22
-11.9	-55.5	23
-22.2	-54.8	24
-17.8	-50.9	25
-15.7	-48.0	26
31.2	121.7	0
51.9	4.4	   0
]

#lat=L[:,1]
#lon=L[:,2]

#coord_x = 6371*sind.(lat).*cosd.(lon)
#coord_y = 6371*sind.(lat).*sind.(lon)
#Lxy = [coord_x coord_y]

#coord_x = 6371*sind.(lat).*cosd.(lon)
#coord_y = 6371*sind.(lat).*sind.(lon)
Lxy = 60 * 1.852 * L[:, 1:2]
Estado = Int.(L[:, 3])


α_safra = [
    2.485 4.277
    2.307 1.118
    2.58 2.267
    3.455 2.5
    3.29 3.073
    3.325 3.752
    2.205 1.532
    1.443 1.84
    0.242 3.544
    0.382 1.0
    0.185 1.0
    0.228 0.071
    1.512 0.6
    1.0 1.591
    4.291 1.497
    6.388 2.555
    2.826 3.85
    2.277 4.995
    6.454 4.244
    3.726 3.877
    6.624 3.917
    8.96 5.278
    8.836 5.031
    7.213 5.744
    7.022 1.0
    6.406 1.0
    9.386 3.549
    4.786 4.102
    3.815 3.265
    4.571 4.683
    3.576 5.065
    8.369 4.598
    7.858 3.582
    6.74 2.399
]


#scatter(L[origem+silos+1:origem+silos+rodo,1],L[origem+silos+1:origem+silos+rodo,2],size="400,400",linewidth="1",
# pointtype="fsquare",plotstyle="linespoints",
# linecolor = "'blue'",Axes(key = :off))

#p=[]
#for i in origem+silos+1:origem+silos+rodo
# push!(p,plot!([L[i,1]],[L[i,2]].+0.5,supp=["$i"],w = "labels"))
#end
#p[end]
#plot!(L[origem+1:origem+silos,1],L[origem+1:origem+silos,2],
#size="400,400",linewidth="1",
#pointtype="fcircle", plotstyle="points",linecolor = "'red'")
#plot!(L[n-destino+1:n,1],L[n-destino+1:n,2],size="400,400",linewidth="3", pointtype="fcircle", plotstyle="points",linecolor = "'violet'")
#plot!(L[origem+1:n-destino-1,1],L[origem+1:n-destino-1,2],size="400,400",linewidth="3", pointtype="fcircle", plotstyle="points",linecolor = "'green'")


β_area = [
    40006 126812
    35014 1373
    3658 262
    4530 600
    171016 27445
    56875 104383
    258433 75676
    461359 22379
    475578 68
    36866 0
    83798 0
    112494 110020
    29475 35
    0 172285
    356168 264668
    848268 370994
    12608 697
    2390 427
    417115 449397
    2955 154055
    58462 784080
    76399 869855
    154095 335417
    130775 17990
    360341 0
    740510 0
    12647 1673690
    16957 2655162
    4757 541506
    6679 57730
    2184 49455
    14500 551338
    238063 1336478
    28600 41563
]



Cap_arm = [
    5700
    7100
    600
    3800
    9700
    20600
    32000
    34000
    4500
    4300
    4100
    3200
    12300
    4900
    30200
    3500
    3000
    3500
    3200
    4800
    3300
    3200
    3500
    3200
    3000
    13800
    3500
    7000
    8900
    3500
    2800
    1800
    26600
    2800
    3200
    3200
    2030
    4400
    23400
    9900
    10300
    10000
    239800
    24600
    25800
    10000
    41600
    27500
    74600
    26800
    33900
    31000
    25800
    400000
    55300
    7300
    16630
    92100
    45900
    47900
    48810
    14000
    6700
    9400
    14000
    42000
    63800
]


Custo_plantio = [
    2912.82 1493.71
    2912.82 1432.975
    2912.82 1432.975
    2912.82 1432.975
    2912.82 1432.975
    2912.82 1372.24
    2326.93 1771.64
    2221.87 1752.538
    2236.06 1752.538
    2236.06 1752.538
    2236.06 1752.538
    2236.06 1752.538
    2236.06 1752.538
    2236.06 1752.538
    2159.38 1752.538
    2825.21 2092.79
    2825.21 2092.79
    2825.21 2092.79
    2825.21 2092.79
    3094.18 1982.643
    3094.18 1982.643
    3094.18 1982.643
    3094.18 1982.643
    3094.18 1982.643
    4564.24 1982.643
    2213.98 1982.643
    2682.066 1817.86
    2682.066 1587.43
    2682.066 1587.43
    2682.066 1587.43
    2682.066 1587.43
    2682.066 1587.43
    2682.066 1915.04
    2682.066 1773.443
]


custo_secagem = [
    27.16666667
    27.16666667
    27.31944444
    27.31944444
    27.45833333
    27.31944444
    27.33333333
    27.33333333
    27.33333333
    27.5
    27.5
    27.5
    27.5
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.33333333
    27.33333333
    27.33333333
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.08333333
    27.08333333
    27.08333333
    27.16666667
    27.25
    26.83333333
    26.83333333
    27.08333333
    27.16666667
    27.16666667
    27.16666667
    27.16666667
    27.16666667
    27.16666667
    27.02777778
]


custo_silo_bolsa = 7.5

custo_estoque = [
    27.16666667
    27.16666667
    27.31944444
    27.31944444
    27.45833333
    27.31944444
    27.33333333
    27.33333333
    27.33333333
    27.5
    27.5
    27.5
    27.5
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.38888889
    27.33333333
    27.33333333
    27.33333333
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.25
    27.08333333
    27.08333333
    27.08333333
    27.16666667
    27.25
    26.83333333
    26.83333333
    27.08333333
    27.16666667
    27.16666667
    27.16666667
    27.16666667
    27.16666667
    27.16666667
    27.02777778
]


consumo_suino = [
    0.203497
    0.203497
    0.203497
    0.203497
    0.203497
    0.203497
    0.203497
    0.203497
    0.156448
    0.203497
    0.203497
    0.208730
    0.203497
    0.203497
    0.203497
    0.201918
    0.200625
    0.203497
    0.180340
    0.204212
    0.227406
    0.226303
    0.231952
    0.183907
    0.216631
    0.203497
]

consumo_ave_corte = [
    0.002441
    0.002441
    0.002441
    0.002441
    0.002441
    0.002441
    0.002441
    0.002441
    0.002527
    0.002441
    0.002441
    0.002718
    0.002441
    0.002441
    0.002441
    0.002873
    0.001745
    0.002441
    0.002987
    0.002296
    0.002355
    0.002481
    0.002466
    0.002701
    0.001698
    0.002441
]

consumo_ave_ovos = [
    0.026778
    0.026778
    0.026778
    0.026778
    0.026778
    0.026778
    0.026778
    0.026778
    0.025748
    0.026778
    0.026778
    0.025763
    0.026778
    0.026778
    0.026778
    0.027918
    0.025685
    0.026778
    0.02792
    0.028791
    0.026778
    0.026778
    0.027353
    0.027353
    0.024475
    0.026778
]

rebanho_suino = [
    224176
    139150
    65507
    49073
    636859
    308422
    1241502
    793301
    1249739
    273518
    179258
    630065
    127441
    138877
    1126310
    5090238
    230748
    64492
    1367512
    7092317
    7099184
    5927862
    1267038
    2538530
    1988478
    168394
]

rebanho_ave_corte = [
    1534660
    2248382
    1857591
    185731
    23090044
    11358917
    7033375
    7625775
    17764366
    1967536
    7467547
    22657693
    6265868
    5423880
    39186162
    100280606
    15478226
    11530347
    149150884
    308694152
    124842367
    115934790
    21371058
    53224506
    56896284
    14428459
]

rebanho_ave_ovos = [
    1493953
    636671
    2599226
    327186
    3288844
    2249177
    2366897
    2188621
    10434505
    2582480
    2176070
    11631875
    1808178
    1566378
    5546645
    20421989
    16163733
    660434
    47942947
    23174302
    15303715
    18776182
    3318814
    10347908
    11824308
    1238664
]

ind_origem_silos = [
    [1, 2],
    [3],
    [4],
    [5],
    [6],
    [7],
    [8, 9],
    [10, 11, 12, 13],
    [14, 15, 16, 17, 18, 19, 20, 21],
    [22, 23, 24, 25, 26, 27],
    [28, 29, 30, 31],
    [32, 33],
    [34, 35],
    [36],
    [37, 38, 39],
    [40, 41, 42, 43, 44, 45, 46],
    [47, 48, 49],
    [50],
    [51, 52],
    [53, 54, 55],
    [53, 54, 55],
    [53, 54, 55],
    [53, 54, 55],
    [53, 54, 55],
    [56],
    [57],
    [58, 59],
    [58, 59],
    [58, 59],
    [58, 59],
    [58, 59],
    [60],
    [61, 62, 63, 64, 65, 66],
    [67],
]


C_porto = 1.1*[
    11284129
    2799188
    244026
    1009920
    771647
    2276916
    892814
    735402
    1290593
    525596
    3506
]

custo_frete_rodo = [
    0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.1467 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.1496 0.1468 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.1465 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.2338 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.1532 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1786 0.1358 0.1358 0.1484 0.3117 0.1358 0.1358
    0.1469 0.1450 0.1358 0.1358 0.1475 0.1455 0.1358 0.1457 0.1483 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1737 0.1358 0.1358 0.1488 0.1777 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.3449 0.1358 0.1358 0.1454 0.1722 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1696 0.1358 0.1358 0.1447 0.1494 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1629 0.1358 0.1358 0.1358 0.1496 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1630 0.1457 0.1358 0.1358 0.1488 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1479 0.1358 0.1358 0.1449 0.1501 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1501 0.1463 0.1358 0.1358 0.1472 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1564 0.1358 0.1358 0.1358 0.1471 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1503 0.1358 0.1358 0.1358 0.1465 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1494 0.1458 0.1358 0.1358 0.1470 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1529 0.1455 0.1358 0.1358 0.1473 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1581 0.1449 0.1358 0.1358 0.1478 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1712 0.1358 0.1358 0.1358 0.1485 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1694 0.1358 0.1358 0.1358 0.1488 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1494 0.1358 0.1358 0.1358 0.1456 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1485 0.1358 0.1358 0.1358 0.1452 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1498 0.1358 0.1358 0.1358 0.1460 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1488 0.1449 0.1358 0.1358 0.1457 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1481 0.1358 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1494 0.1450 0.1358 0.1358 0.1463 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1477 0.1452 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1473 0.1448 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1485 0.1462 0.1358 0.1358 0.1455 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1491 0.1453 0.1358 0.1358 0.1456 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1484 0.1464 0.1358 0.1358 0.1455 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1469 0.1455 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1468 0.1468 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1473 0.1471 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1473 0.1483 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1453 0.1358 0.1358 0.1358 0.1474 0.1500 0.1358 0.1358 0.1358 0.1358 0.1358
    0.1358 0.1358 0.1358 0.1358 0.1479 0.1487 0.1358 0.1358 0.1449 0.1358 0.1358
    0.1453 0.1358 0.1358 0.1358 0.1472 0.1486 0.1358 0.1358 0.1446 0.1358 0.1358
    0.1502 0.1478 0.1358 0.1472 0.1358 0.1592 0.1358 0.1358 0.1358 0.1461 0.1470
    0.1674 0.1637 0.1465 0.1573 0.1358 0.1551 0.1358 0.1358 0.1358 0.1500 0.1544
    0.1559 0.1589 0.1461 0.1524 0.1358 0.1674 0.1358 0.1358 0.1358 0.1496 0.1503
    0.1708 0.1503 0.1450 0.1497 0.1358 0.1502 0.1358 0.1358 0.1358 0.1487 0.1495
    0.1820 0.1627 0.1464 0.1563 0.1358 0.1633 0.1358 0.1358 0.1358 0.1500 0.1541
    0.1780 0.1536 0.1456 0.1503 0.1358 0.1499 0.1358 0.1358 0.1358 0.1492 0.1501
    0.1694 0.1498 0.1358 0.1493 0.1358 0.1629 0.1358 0.1358 0.1358 0.1482 0.1491
    0.1503 0.1480 0.1358 0.1474 0.1358 0.3502 0.1358 0.1358 0.1358 0.1463 0.1472
    0.1556 0.1487 0.1358 0.1482 0.1358 0.2944 0.1358 0.1358 0.1358 0.1470 0.1480
    0.1501 0.1478 0.1358 0.1472 0.1358 0.3028 0.1358 0.1358 0.1358 0.1460 0.1470
    0.1811 0.1551 0.1457 0.1505 0.1358 0.2931 0.1358 0.1358 0.1358 0.1493 0.1503
    0.1669 0.1784 0.1480 0.1751 0.1358 0.1811 0.1358 0.1358 0.1358 0.1622 0.1713
    0.1355 0.1746 0.1476 0.1694 0.1358 0.1485 0.1358 0.1358 0.1358 0.1565 0.1672
    0.1725 0.1337 0.1486 0.1772 0.1358 0.1476 0.1358 0.1358 0.1358 0.1643 0.1750
    0.1789 0.2581 0.1502 0.2231 0.1358 0.1479 0.1358 0.1358 0.1358 0.1373 0.2090
    0.1715 0.1817 0.1486 0.1764 0.1358 0.1475 0.1358 0.1358 0.1358 0.1636 0.1743
    0.1624 0.1479 0.1672 0.1696 0.1358 0.1461 0.1358 0.1358 0.1358 0.1576 0.1904
    0.1499 0.1666 0.2090 0.1735 0.1358 0.1358 0.1358 0.1358 0.1358 0.1895 0.1799
    0.1481 0.1472 0.1358 0.1467 0.1358 0.1452 0.1358 0.1447 0.1358 0.1456 0.1465
    0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.1358 0.1483 0.1358 0.1358 0.1358
    0.1500 0.1499 0.1467 0.1494 0.1358 0.1450 0.1358 0.1358 0.1358 0.1483 0.1492
    0.1503 0.1490 0.1358 0.1485 0.1358 0.1482 0.1358 0.1358 0.1358 0.1473 0.1483
    0.1500 0.1482 0.1358 0.1477 0.1358 0.1476 0.1358 0.1358 0.1447 0.1466 0.1476
    0.1501 0.1486 0.1358 0.1482 0.1358 0.1482 0.1358 0.1358 0.1451 0.1470 0.1480
    0.1528 0.1477 0.1358 0.1472 0.1449 0.1475 0.1358 0.1358 0.1456 0.1461 0.1470
    0.1552 0.1491 0.1358 0.1486 0.1358 0.1486 0.1358 0.1358 0.1358 0.1475 0.1484
    0.1542 0.1490 0.1358 0.1485 0.1358 0.1485 0.1358 0.1358 0.1358 0.1474 0.1483
    0.1499 0.1478 0.1358 0.1473 0.1358 0.1491 0.1358 0.1358 0.1452 0.1462 0.1471
]

#--------------------------------------------------------------------

Aliquota_ICMS = [
    0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.17 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.17 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.17 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.18 0.07 0.12 0.12 0.12 0.12 0.12 0.07 0.07 0.07 0.07 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.17 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.00
    0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.12 0.07 0.20 0.12 0.12 0.12 0.12 0.07 0.07 0.07 0.07 0.00
    0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.12 0.07 0.12 0.18 0.12 0.12 0.12 0.07 0.07 0.07 0.07 0.00
    0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.12 0.07 0.12 0.12 0.18 0.12 0.12 0.07 0.07 0.07 0.07 0.00
    0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.12 0.07 0.12 0.12 0.12 0.17 0.12 0.07 0.07 0.07 0.07 0.00
    0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.07 0.12 0.07 0.12 0.12 0.12 0.12 0.18 0.07 0.07 0.07 0.07 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.17 0.12 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.17 0.12 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.17 0.12 0.00
    0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.18 0.00
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
]

J=[]
for e in 1:estados
    push!(J,findall(L[:,3] .== e))
end


#include("Vetores_dados.jl")

#
n_silo=2 #Pensar em mudar para 1 - Número de arestas do silo ao rodo
n_linha=2 #Pensar em mudar para 1 - Número de arestas do rodo ao destino interno
n_ferro_rodo=2 #Pensar em mudar para 1 - Número de arestas de cada ferro ao rodo
n_hidro_rodo=2 #Pensar em mudar para 1 - Número de arestas de cada ferro ao hidro
n_porto=2 #Pensar em mudar para 1 - Número de arestas que ligam o rodo ao porto
#
c = round.([norm(Lxy[i,:] - Lxy[j,:]) for i in 1:n,j in 1:n])
#
α=zeros(origem,safra,length(T))
for e in 1:origem, s in 1:safra, t in T
    α[e,s,t] = α_safra[e,s]
end

Estados_seca = [19,20,21,22,23,24,25,26]
α[Estados_seca,:,:] = 0.8*α[Estados_seca,:,:]
#
β=zeros(origem,safra,length(T))
for e in 1:origem, s in 1:safra, t in T
    β[e,s,t] = β_area[e,s]
end
#
δ=zeros(n,length(T))
for j in Silos, t in T
    δ[j,t] = Cap_arm[j-origem]
end
#https://www.econstor.eu/bitstream/10419/121587/1/797027513.pdf
Cap_min_trem=5000
Cap_max_trem=2.5e6
Cap_min_navio=5000
Cap_max_navio=5e6
Cap_max_rodo=18e6

ε=0*ones(n,n,length(T))
for j in setdiff(1:n,1:origem), k in 1:n, t in T
    if j >= origem+1 && j<=  origem+silos+rodo
        ε[j,k,t] = Cap_max_rodo
    end
    if j >= origem+silos+rodo+1 && j<=  origem+silos+rodo+ferro
        ε[j,k,t] = Cap_max_trem
    end
    if j >= origem+silos+rodo+ferro+1 && j<=  origem+silos+rodo+ferro+hidro
        ε[j,k,t] = Cap_max_navio
    end
    if j>= origem+silos+rodo+ferro+hidro+1 && k <= origem+silos+rodo+ferro+hidro
        ε[j,k,t] = Cap_max_rodo
    end
end
#
γ1 = zeros(origem,safra,length(T))
for e in 1:origem, s in 1:safra, t in T
    γ1[e,s,t] = Custo_plantio[e,s]
end
#
γ2 = zeros(n,length(T))
γ3 = zeros(n,length(T))
γ4 = zeros(n,length(T))
for j in Silos, t in T
    γ2[j,t] = custo_secagem[j-origem]
    γ3[j,t] = custo_silo_bolsa
    γ4[j,t] = custo_estoque[j-origem]
end
#
I0=zeros(n)
for j in Silos
    I0[j] = 0
end
#
γ5=zeros(n,n,length(T))
for j in origem+silos+rodo+ferro+hidro+1:origem+silos+rodo+ferro+hidro+portos, k in 1:n, t in T
    γ5[j,k,t] = 1060
end
#
σ=zeros(n,length(T))
for k in 1:length(D_i), t in T
    σ[k+origem+silos+rodo+ferro+hidro+portos,t] = consumo_suino[k]*rebanho_suino[k] + consumo_ave_corte[k]*rebanho_ave_corte[k]
    +consumo_ave_ovos[k]*rebanho_ave_ovos[k]
end
for k in 1:length(D_e), t in T
    σ[k+origem+silos+rodo+ferro+hidro+portos+dest_in,t] = 21833737/2
end
σ[279,1]=σ[279,1]+5087.3
σ[281,1]=σ[281,1]+731.4
σ[276,1]=σ[276,1]+38.4

#
Cap_porto=zeros(n)
Cap_porto[origem+silos+rodo+ferro+hidro+1:origem+silos+rodo+ferro+hidro+portos] .= C_porto
#
#frete_rodo=0.14390
#frete_ferro=0.07180
#frete_hidro=0.04500
##http://www.ppe.ufrj.br/images/publica%C3%A7%C3%B5es/mestrado/Raphael_Benirschke_Terra.pdf
frete_rodo=0.193229*1.6204
frete_ferro=0.072*1.6204
frete_hidro=0.03181*1.6204
frete_hidro_ex=0.5*0.03181*1.6204
#
Preco_ICMS = 705
Custo_ICMS = zeros(estados, estados)
for e = 1:estados, elinha = 1:estados
    Custo_ICMS[e, elinha] = Preco_ICMS * Aliquota_ICMS[e,elinha]
end

#Custo_ICMS = zeros(n, estados)
#for i = 1:n, j = 1:estados
#    Custo_ICMS[i, j] = Preco_ICMS * Aliquota_ICMS[Int(L[i,3]),j]
#end

λ = [
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.6
0.8
0.8
0.8
0.8
0.77
0.77
0.77
0.77
0.77
1
0.6
0.9
0.9
0.9
0.9
0.9
0.9
0.9
0.9
]

λ[[19,20,21,22,23,24,25,26]] .= 0.8*λ[[19,20,21,22,23,24,25,26]]

perda = 0.02

#--------------------------------------------------------------------------------
include("Matriz_adjacencias.jl")
#--------------------------------------------------------------------------------
Intermodais = zeros(n,n)
Intermodais[origem+silos+1:origem+silos+rodo+ferro+hidro,origem+silos+1:origem+silos+rodo+ferro+hidro] .= A_modais
for i in origem+silos+1:origem+silos+rodo, j in origem+silos+1:origem+silos+rodo
    Intermodais[i,j] = 0
end
for i in origem+silos+rodo+1:origem+silos+rodo+ferro, j in origem+silos+rodo+1:origem+silos+rodo+ferro
    Intermodais[i,j] = 0
end
for i in origem+silos+rodo+ferro+1:origem+silos+rodo+ferro+hidro, j in origem+silos+rodo+ferro+1:origem+silos+rodo+ferro+hidro
    Intermodais[i,j] = 0
end
#https://ontl.epl.gov.br/aplicacoes/simulador-de-custo-de-transbordo/
cT = zeros(n,n)
for i in 1:n, j in 1:n
    if i>= origem+silos+1 && i<=origem+silos+rodo && j>= origem+silos+rodo+1 && j<=origem+silos+rodo+ferro
        cT[i,j] = 9.74
    end
    if i>= origem+silos+1 && i<=origem+silos+rodo && j>= origem+silos+rodo+ferro+1 && j<=origem+silos+rodo+ferro+hidro
        cT[i,j] = 9.95
    end
    if i>= origem+silos+rodo+1 && i<=origem+silos+rodo+ferro && j>= origem+silos+1 && j<=origem+silos+rodo
        cT[i,j] = 8.89
    end
    if i>= origem+silos+rodo+1 && i<=origem+silos+rodo+ferro &&  j>= origem+silos+rodo+ferro+1 && j<=origem+silos+rodo+ferro+hidro
        cT[i,j] = 9.73
    end
    if  i>= origem+silos+rodo+ferro+1 && i<=origem+silos+rodo+ferro+hidro &&  j>= origem+silos+1 && j<=origem+silos+rodo
        cT[i,j] = 10.67
    end
    if  i>= origem+silos+rodo+ferro+1 && i<=origem+silos+rodo+ferro+hidro && j>= origem+silos+rodo+1 && j<=origem+silos+rodo+ferro
        cT[i,j] = 11.82
    end
end
#----------------------------------------------------------------------------
#http://www.ppe.ufrj.br/images/publica%C3%A7%C3%B5es/mestrado/Raphael_Benirschke_Terra.pdf

τ = zeros(n, n)
for i = 1:n, j = 1:n
    if i >= 1 && i <= origem + silos + rodo
        #3km/l, 1l de diesel = 3.13 kgCO2, 57 ton
        τ[i, j] = (15.4*3.13/1000)/1000 #(0.5*3.13/57)/1000
    end
    if i >= origem + silos + rodo + 1 && i <= origem + silos + rodo + ferro && j >= origem + silos + rodo + 1 && j <= origem + silos + rodo + ferro
        #9l/km, 1l de diesel = 3.13 kgCO2, 5000 ton
        τ[i, j] = (5.7*3.13/1000)/1000 #(9*3.13/3000)/1000
    end
    if i >= origem + silos + rodo + ferro + 1 &&  i <= origem + silos + rodo + ferro + hidro && j >= origem + silos + rodo + ferro + 1 &&  j <= origem + silos + rodo + ferro + hidro
        τ[i, j] = (4.1*3.13/1000)/1000 #0.01252 / 1000 #0.00533227/1000
    end
    if i >= origem + silos + rodo + ferro + hidro + 1 &&  i <= origem + silos + rodo + ferro + hidro + portos && j >= origem + silos + rodo + ferro + hidro + portos + estados+1
        #2l/1000 t por km
        τ[i, j] = (4.1*3.13/1000)/1000 #0.01252 / 1000 #0.00533227/1000#0.4*0.01252 / 1000
    end
    if i >= origem + silos + rodo + ferro + hidro + 1 &&  i <= origem + silos + rodo + ferro + hidro + portos && j <= origem + silos + rodo + ferro + hidro + portos + estados
        τ[i, j] = (15.4*3.13/1000)/1000 #(0.5*3.13/57)/1000
    end
    if Intermodais[i,j] == 1
        τ[i,j] = (15.4*3.13/1000)/1000
    end
end
μ = zeros(n, n)
for i = 1:n, j = 1:n
    μ[i, j] = τ[i, j] * c[i, j]
end
#
u=zeros(n,n)
l=zeros(n,n)
for i in 1:n, j in 1:n
    u[i,j] = maximum(ε)
    l[i,j] = 5000
end
dist_max_imp = 500
#
