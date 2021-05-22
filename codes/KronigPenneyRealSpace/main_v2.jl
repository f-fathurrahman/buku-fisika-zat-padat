using Printf
using LinearAlgebra
using LaTeXStrings

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)
#plt.rc("font", size=16)

include("init_FD1d_p_grid.jl")
include("build_D1_matrix_p_3pt.jl")
include("build_D2_matrix_p_3pt.jl")

function pot_kronig_penney( x; V0 = 1.0, L=1.0 )
    if (x >= L/4) && (x < 3*L/4)
        return V0
    else
        return 0.0
    end
end

function pot_mathieu(x; V0=1.0, L=1.0 )
    return V0*( 1 + cos(2*pi*x/L) )
end

function build_Ham_matrix( D1, D2, Vpot, k::Float64 )
    N = size(Vpot,1)
    Ham = -0.5*( D2 + 2*im*k*D1 - k^2*diagm(0 => ones(N)) ) + diagm( 0 => Vpot )
    return Hermitian(Ham)
end

function main()
    # Initialize the grid points
    xmin =  0.0
    xmax = 10.0
    L = xmax - xmin
    N = 101
    x, h = init_FD1d_p_grid(xmin, xmax, N)
    # Build derivative matrices
    D1 = build_D1_matrix_p_3pt(N, h)
    D2 = build_D2_matrix_p_3pt(N, h)
    # Potential
    Vpot = pot_kronig_penney.(x, L=L, V0=0.05)
    #Vpot = pot_mathieu.(x, L=L, V0=0.1)

    plot_title = "N="*string(N)
    plt.clf()

    Nk = 51*4
    k = range(-3*pi/L, 3*pi/L, length=Nk)
    Nstates = 3
    ebands = zeros(Nk, Nstates)
    for ik in 1:Nk
        # Hamiltonian
        Ham = build_Ham_matrix(D1, D2, Vpot, k[ik])
        # Solve the eigenproblem
        evals = eigvals( Ham )
        for ist in 1:Nstates
            ebands[ik,ist] = evals[ist]
        end
    end

    plt.figure(figsize=(7,10))
    for ist in 1:Nstates
        labelstr = "band-"*string(ist)
        plt.plot(k, ebands[:,ist], label=labelstr, marker="o", markersize=0.5)
    end
    E_free = 0.5*k.^2 .+ minimum(ebands[:,1])
    plt.plot(k, E_free, label="free-electron", linestyle="--")
    plt.legend()
    plt.grid()
    #plt.tight_layout()
    plt.title("Kronig-Penney Model")
    plt.xlabel("k")
    plt.ylabel("E (Ha)")
    plt.xlim(k[1], k[end])
    plt.xticks([-3*pi/L, -2*pi/L, -pi/L, 0, pi/L, 2*pi/L, 3*pi/L],
        [L"$-3\pi/L$", L"$-2\pi/L$", L"$-\pi/L$", "0", L"$\pi/L$", L"$2\pi/L$", L"$3\pi/L$"])
    plt.savefig("IMG_bands_kronigpenney_"*string(N)*"_ext.pdf")

end

main()