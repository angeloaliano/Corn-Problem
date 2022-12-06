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

include("Vetores_dados.jl")

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
