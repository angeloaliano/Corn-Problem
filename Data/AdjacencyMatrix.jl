A_origem_silos=zeros(origem,silos)
for i in 1:origem
    for j in ind_origem_silos[i]
        A_origem_silos[i,j] = 1
    end
end
#
A_silos_rodo=zeros(silos,rodo)
for i in 1:silos
    ind_per=sortperm(c[i+origem,origem+silos+1:origem+silos+rodo])[1:n_silo]
    A_silos_rodo[i,ind_per] .= 1
end

A_silos_ferro=zeros(silos,ferro)
for i in 1:silos
    ind_per=findall(c[i+origem,origem+silos+rodo+1:origem+silos+rodo+ferro] .== 0)
    A_silos_ferro[i,ind_per] .= 1
end

A_silos_hidro=zeros(silos,hidro)
for i in 1:silos
    ind_per=findall(c[i+origem,origem+silos+rodo+ferro+1:origem+silos+rodo+ferro+hidro] .== 0)
    A_silos_hidro[i,ind_per] .= 1
end

A_silos_modais = Int.([A_silos_rodo A_silos_ferro A_silos_hidro])

#-----------------------------------------------------------------------
ind_rodo_rodo=[
(102,103)
(103,104)
(104,105)
(105,106)
(106,107)
(107,108)
(108,109)
(109,110)
(110,111)
(111,112)
(112,113)
(113,114)
(114,115)
(116,117)
(117,118)
(118,119)
(119,105)
(105,106)
(106,107)
(107,120)
(120,121)
(121,122)
(122,123)
(123,124)
(124,125)
(125,126)
(126,127)
(128,129)
(129,130)
(130,131)
(131,132)
(132,133)
(133,134)
(134,135)
(135,106)
(136,137)
(137,138)
(138,139)
(139,140)
(140,141)
(141,142)
(142,143)
(143,126)
(144,145)
(145,146)
(146,147)
(147,132)
(132,148)
(148,149)
(149,150)
(151,152)
(152,153)
(153,154)
(154,134)
(134,155)
(155,156)
(156,157)
(157,158)
(158,140)
(159,160)
(160,161)
(161,162)
(162,133)
(133,163)
(163,164)
(164,165)
(165,166)
(166,140)
(167,168)
(168,169)
(169,136)
(136,129)
(170,171)
(171,172)
(172,173)
(173,174)
(174,175)
(175,176)
(177,146)
(146,162)
(162,134)
(134,178)
(178,170)
(170,108)
(108,120)
(153,105)
(104,117)
]

A_rodo_rodo=zeros(rodo,rodo)
for k in 1:length(ind_rodo_rodo)
    A_rodo_rodo[ind_rodo_rodo[k][1] - origem - silos,ind_rodo_rodo[k][2] - origem - silos] = 1
    A_rodo_rodo[ind_rodo_rodo[k][2] - origem - silos,ind_rodo_rodo[k][1] - origem - silos] = 1
end
##############################
#Novo - conectando as federais com menos de 150km
fed_fed=[]
for i in 1:rodo
    dist = c[origem+silos+1:origem+silos+rodo,origem+silos+i]
    ind_per = findall(dist .< 150) .+ origem .+ silos
    for j in 1:length(ind_per)
        if i+origem+silos != ind_per[j]
            push!(fed_fed,(i+origem+silos,ind_per[j]))
        end
    end
end
for k in 1:length(fed_fed)
    A_rodo_rodo[fed_fed[k][1] - origem - silos,fed_fed[k][2] - origem - silos] = 1
    A_rodo_rodo[fed_fed[k][2] - origem - silos,fed_fed[k][1] - origem - silos] = 1
end
A_rodo_rodo=Int.(A_rodo_rodo)
#-----------------------------------------------------------------------
ind_ferro_ferro=[
(179,180)
(180,182)
(181,180)
(182,183)
(185,186)
(187,186)
(186,184)
(186,183)
(188,189)
(189,185)
(190,186)
(190,191)
(191,192)
(192,193)
(194,195)
(195,196)
(196,197)
(197,191)
(194,195)
(195,191)
(195,197)
(197,191)
(196,197)
(195,198)
(198,199)
(197,199)
(200,201)
(201,193)
(193,202)
(202,203)
(203,204)
(205,206)
(206,190)
(207,208)
(208,209)
(209,210)
(210,211)
(212,213)
(213,214)
(214,207)
(204,215)
(215,213)
]
#-------------------------------------------------------------------------
A_ferro_ferro=zeros(ferro,ferro)
for k in 1:length(ind_ferro_ferro)
    A_ferro_ferro[ind_ferro_ferro[k][1]-origem-silos-rodo,ind_ferro_ferro[k][2]-origem-silos-rodo] = 1
    A_ferro_ferro[ind_ferro_ferro[k][2]-origem-silos-rodo,ind_ferro_ferro[k][1]-origem-silos-rodo] = 1
end
A_ferro_ferro=Int.(A_ferro_ferro)
#-------------------------------------------------------------------------
ind_hidro_hidro=[
(216,217)
(217,218)
(218,219)
(220,221)
(221,222)
(222,218)
(223,224)
(224,225)
(226,227)
(227,228)
(229,230)
(230,231)
(231,232)
(232,233)
(234,232)
(235,236)
(236,237)
(238,239)
(241,239)
(242,241)
(240,239)
(243,239)
(239,228)
(228,245)
(244,245)
]

A_hidro_hidro=zeros(hidro,hidro)
for k in 1:length(ind_hidro_hidro)
    A_hidro_hidro[ind_hidro_hidro[k][1]-origem-silos-rodo-ferro,ind_hidro_hidro[k][2]-origem-silos-rodo-ferro] = 1
    A_hidro_hidro[ind_hidro_hidro[k][2]-origem-silos-rodo-ferro,ind_hidro_hidro[k][1]-origem-silos-rodo-ferro] = 1
end
A_hidro_hidro=Int.(A_hidro_hidro)
#-------------------------------------------------------------------------
A_rodo_ferro=zeros(rodo,ferro)
for j in 1:ferro
    ind_per = sortperm(c[origem+silos+1:origem+silos+rodo,j+origem+silos+rodo])[1:n_ferro_rodo]
    A_rodo_ferro[ind_per,j] .= 1
end
A_rodo_ferro=Int.(A_rodo_ferro)
A_ferro_rodo=A_rodo_ferro'
#-------------------------------------------------------------------------
A_rodo_hidro=zeros(rodo,hidro)
for j in 1:hidro
    ind_per = sortperm(c[origem+silos+1:origem+silos+rodo,j+origem+silos+rodo+ferro])[1:n_hidro_rodo]
    A_rodo_hidro[ind_per,j] .= 1
end
A_rodo_hidro=Int.(A_rodo_hidro)
A_hidro_rodo=A_rodo_hidro'
#-------------------------------------------------------------------------
A_ferro_hidro=zeros(ferro,hidro)
for j in 1:hidro
    ind_per = findall(c[origem+silos+rodo+1:origem+silos+rodo+ferro,j+origem+silos+rodo+ferro] .< 10)
    A_ferro_hidro[ind_per,j] .= 1
end
A_ferro_hidro=Int.(A_ferro_hidro)
A_hidro_ferro=A_ferro_hidro'
#-------------------------------------------------------------------------
A_rodo_porto=zeros(rodo,portos)
for j in 1:portos
    ind_per = sortperm(c[origem+silos+1:origem+silos+rodo,j+origem+silos+rodo+ferro+hidro])[1:n_porto]
    A_rodo_porto[ind_per,j] .= 1
end
A_rodo_porto=Int.(A_rodo_porto)
#-------------------------------------------------------------------------
A_ferro_porto=zeros(ferro,portos)
for j in 1:portos
    ind_per = findall(c[origem+silos+rodo+1:origem+silos+rodo+ferro,j+origem+silos+rodo+ferro+hidro] .< 10)
    if length(ind_per) > 0
        A_ferro_porto[ind_per,j] .= 1
    end
end
A_ferro_porto=Int.(A_ferro_porto)
#-------------------------------------------------------------------------
A_hidro_porto=zeros(hidro,portos)
for j in 1:portos
    ind_per = findall(c[origem+silos+rodo+ferro+1:origem+silos+rodo+ferro+hidro,j+origem+silos+rodo+ferro+hidro] .< 10)
    if length(ind_per) > 0
        A_hidro_porto[ind_per,j] .= 1
    end
end
A_hidro_porto=Int.(A_hidro_porto)
#-------------------------------------------------------------------------
A_linha_1 = zeros(rodo,dest_in)
for j in 1:dest_in
    ind_per = sortperm(c[origem+silos+1:origem+silos+rodo,j+origem+silos+rodo+ferro+hidro+portos])[1:n_linha]
    A_linha_1[ind_per,j] .= 1
end
#
A_linha_2 = zeros(ferro,dest_in)
for j in 1:dest_in
    ind_per = findall(c[origem+silos+rodo+1:origem+silos+rodo+ferro,j+origem+silos+rodo+ferro+hidro+portos] .< 10)
    if length(ind_per) > 0
        A_linha_2[ind_per,j] .= 1
    end
end
#
A_linha_3 = zeros(hidro,dest_in)
for j in 1:dest_in
    ind_per = findall(c[origem+silos+rodo+ferro+1:origem+silos+rodo+ferro+hidro,j+origem+silos+rodo+ferro+hidro+portos] .< 10)
    if length(ind_per) > 0
        A_linha_3[ind_per,j] .= 1
    end
end
#
A_porto_exp = ones(portos,dest_ex)
#
A_modais = [A_rodo_rodo A_rodo_ferro A_rodo_hidro;
            A_ferro_rodo A_ferro_ferro A_ferro_hidro;
            A_hidro_rodo A_hidro_ferro A_hidro_hidro]
#
A_modais_porto = [A_rodo_porto;A_ferro_porto;A_hidro_porto]
#A_porto_modais = A_modais_porto'
#
A_linha = [A_linha_1;A_linha_2;A_linha_3]
#
A = zeros(n,n)
A[1:origem,origem+1:origem+silos] .= A_origem_silos
A[origem+1:origem+silos,origem+silos+1:origem+silos+rodo+ferro+hidro] .= A_silos_modais
A[origem+silos+1:origem+silos+rodo+ferro+hidro,origem+silos+1:origem+silos+rodo+ferro+hidro] .= A_modais
A[origem+silos+1:origem+silos+rodo+ferro+hidro,origem+silos+rodo+ferro+hidro+1:origem+silos+rodo+ferro+hidro+portos] .= A_modais_porto
#A[origem+silos+rodo+ferro+hidro+1:origem+silos+rodo+ferro+hidro+portos,origem+silos+1:origem+silos+rodo+ferro+hidro] .= A_porto_modais
A[origem+silos+1:origem+silos+rodo+ferro+hidro,origem+silos+rodo+ferro+hidro+portos+1:origem+silos+rodo+ferro+hidro+portos+dest_in] .= A_linha
A[origem+silos+rodo+ferro+hidro+1:origem+silos+rodo+ferro+hidro+portos,origem+silos+rodo+ferro+hidro+portos+dest_in+1:origem+silos+rodo+ferro+hidro+portos+dest_in+dest_ex] .= A_porto_exp
#
AlinhaO=zeros(n,n)
AlinhaO[1:origem,origem+1:origem+silos] .= A_origem_silos
#------------------------------------------------------------

#Nos_iguais=[]
#for i in origem+silos+1:origem+silos+rodo
#    nos_iguais_i=[]
#    for j in origem+silos+1:origem+silos+rodo+ferro+hidro
#        if L[i,1] == L[j,1] && L[i,2] == L[j,2]
#            push!(nos_iguais_i,j)
#        end
#    end
#    push!(Nos_iguais,nos_iguais_i)
    #Nos_iguais=[Nos_iguais;[nos_iguais_i]]
#end
#
