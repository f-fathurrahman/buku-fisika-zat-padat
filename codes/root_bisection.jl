function root_bisection(f, x1, x2, TOL=1e-10, NiterMax=100, verbose=true)
    f1 = f(x1)
    f2 = f(x2)

    if f1*f2 > 0
        verbose && @printf("f1 = %18.10f", f1)
        verbose && @printf("f2 = %18.10f", f2)
        verbose && error("f1 dan f2 memiliki tanda yang sama")
    end
    
    xr = Inf
    for i in 1:NiterMax
        xr = 0.5*(x1 + x2)
        fxr = f(xr)
        if abs(fxr) <= TOL
            verbose && println("Iterasi konvergen: akar ditemukan")
            return xr
        end        
        verbose && @printf("Iter = %5d %18.10f %15.5e\n", i, xr, abs(fxr))
        if f1*fxr < 0.0
            x2 = xr
            f2 = fxr
        else
            x1 = xr
            f1 = fxr
        end
    end
    
    @printf("WARNING: Konvergensi tidak diperoleh setelah %d iterasi\n", NiterMax)
    @printf("WARNING: Nilai tebakan akhir akan dikembalikan\n")
    return xr
end