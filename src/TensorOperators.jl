import WignerSymbols: wigner3j, wigner6j


# Check the triangle condition
δtri(j1,j2,j3) = (j3 <= j1+j2) && (j1 <= j2+j3) && (j2 <= j3+j1) && isinteger(j1+j2+j3)
δsum(m1,m2,m3) = (m1+m2+m3 == 0)


# Compute Wigner symbols
function wigner3j_(j1,j2,j3,m1,m2,m3)
    if δtri(j1,j2,j3) && δsum(m1,m2,m3)
        try
            wigner3j(Float64,j1,j2,j3,m1,m2,m3)
        catch
            0.0
        end
    else
        0.0
    end
end

function wigner6j_(j1,j2,j3,m1,m2,m3)
    try
        wigner6j(Float64,j1,j2,j3,m1,m2,m3)
    catch
        0.0
    end
end


# Just let this one use the standard WignerSybmols 6j call.
function wigner9j_(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    val = 0.0
    kmin = maximum([abs(j1-j9), abs(j4-j8), abs(j2-j6)])
    kmax = minimum([abs(j1+j9), abs(j4+j8), abs(j2+j6)])
    if kmax >= kmin
        val += sum(
            (-1)^(2k) * (2k+1) *
            wigner6j(Float64,j1,j4,j7,j8,j9,k) * 
            wigner6j(Float64,j2,j5,j8,j4,k,j6) *
            wigner6j(Float64,j3,j6,j9,k,j1,j2) for k in kmin:kmax)
    end
    return val
end