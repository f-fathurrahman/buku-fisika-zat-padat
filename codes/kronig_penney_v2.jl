using Printf

import PyPlot
const plt = PyPlot

include("root_bisection.jl")

struct KronigPenneyParams
    a::Float64
    b::Float64
    U_0::Float64
end

# Default values
function KronigPenneyParams( ; a=10.0, b=1.0, U_0=1.0 )
    return KronigPenneyParams(a, b, U_0)
end

function kronig_penney_model(
    E::Float64, k::Float64;
    kp_params=KronigPenneyParams()
)
    a = kp_params.a
    b = kp_params.b
    U_0 = kp_params.U_0
    #
    κ = sqrt(2*(U_0 - E))
    K = sqrt(2*E)
    res = (κ^2 - K^2)/(2*κ*K) * sinh(κ*b) * sin(K*a) + cosh(κ*b)*cos(K*a) - cos(k*(a+b))
    return res
end

function scan_for_root(f, E0, EN; ΔE=0.1)
    res = []
    Ninterval = round( Int64, (EN-E0)/ΔE )
    println("Ninterval = ", Ninterval)
    for i in 1:Ninterval
        a = E0 + (i-1)*ΔE
        b = E0 + i*ΔE
        @printf("Scan a = %f, b = %f\n", a, b)
        fa = f(a)
        fb = f(b)
        if fa*fb < 0.0
            #@printf("Akar ada di antara a dan b\n")
            #println("fa = ", f(a))
            #println("fb = ", f(b))
            push!(res, (a,b))
        end
    end
    return res
end

function do_plot(f, Egrid, filesave)
    plt.clf()
    NptsPlot = length(Egrid)
    fE = zeros(Float64,NptsPlot)
    for i in 1:NptsPlot
        fE[i] = f(Egrid[i])
    end
    plt.plot(Egrid, fE)
    plt.grid(true)
    plt.savefig(filesave)
end

function find_roots()
    kp_params = KronigPenneyParams(a=5.0, b=5.0, U_0=0.1)
    a = kp_params.a
    b = kp_params.b
    U_0 = kp_params.U_0

    k1 = -pi/(a+b)
    k2 = pi/(a+b)
    Nk = 51
    Δk = (k2-k1)/(Nk-1)

    plt.clf()

    Nroots_plt = 3

    Egrid = range(1e-10, U_0, length=1000)
    for ik in 1:Nk
        k = k1 + (ik-1)*Δk
        println("k = ", k)
        f(E) = kronig_penney_model(E, k, kp_params=kp_params)
        do_plot(f, Egrid, "IMG_kp_k_"*string(ik)*".png")


        #res = scan_for_root(f, 1e-10, 0.9*U_0, ΔE=1e-3*U_0)
        #Nroots = length(res)
        ##println("Nroots = ", Nroots)
        ##if Nroots < Nroots_plt
        ##    error("Too few bands found for k = ", k)
        ##end
        #Ebands = []
        #for p in res
        ##for ist in 1:Nroots_plt
        ##    p = res[ist]
        #    E1 = p[1]
        #    E2 = p[2]
        #    xr = root_bisection(f, E1, E2)
        #    push!(Ebands, xr)
        #    plt.plot( [k], [xr], marker="o" )
        #end
        #println("Ebands = ", Ebands)
    end
    #plt.grid(true)
    #plt.savefig("IMG_bands_kp_1.pdf")
end
find_roots()