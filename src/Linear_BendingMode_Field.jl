using Parameters
using CompositeStructs
using InfiniteArrays


export 
    LinearCaseA_Bend_Field, LinearCaseB_Bend_Field,
    Rotation, SpinRotation, HyperfineFermi, HyperfineIS,
    lDoubling

Base.@kwdef struct LinearCaseA_Bend_Field <: AbstractCaseA #LinearBasisState
    Λ::Rational{Int64}
    l::Rational{Int64}
    S::Rational{Int64}
    Σ::Rational{Int64}
    J::Rational{Int64}
    Ω::Rational{Int64}
    I::Rational{Int64}
    F::Rational{Int64}
    M::Rational{Int64}

    constraints = (
        Σ = -S:S,
        J = abs(Σ+l):∞,
        Ω = max(-J, Λ+Σ):min(J, Λ+Σ),
        F = abs(J-I):abs(J+I),
        M = -F:F
    )
end


Base.@kwdef struct LinearCaseB_Bend_Field <: AbstractCaseB #LinearBasisState
    Λ::Rational{Int64}
    l::Rational{Int64}
    N::Rational{Int64}
    S::Rational{Int64}
    J::Rational{Int64}
    I::Rational{Int64}
    F::Rational{Int64}
    M::Rational{Int64}

    constraints = (
        N = abs(Λ+l):∞,
        J = abs(N-S):abs(N+S),
        F = abs(J-I):abs(J+I),
        M = -F:F
    )
end


function unpack(state::LinearCaseA_Bend_Field)
    return (state.Λ, state.l, state.S, state.Σ, state.J, state.Ω, state.I, state.F, state.M)
end

function unpack(state::LinearCaseB_Bend_Field)
    return (state.Λ, state.l, state.N, state.S, state.J, state.I, state.F, state.M)
end


function Rotation(state::LinearCaseA_Bend_Field, state′::LinearCaseA_Bend_Field)
    """
    Rotational hamiltonian.
    """

    Λ, l, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, l′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    ME = 0.0
    if δ(l,l′) && δ(J,J′) && δ(F,F′) && δ(M,M′)
        if δ(Σ,Σ′)
            ME += J*(J+1)-(l+Σ)^2 + S*(S+1) - Σ^2
        elseif δ(Σ,Σ′-1)
            ME += -1.0*sqrt(J*(J+1) - (l+Σ)*(l+Σ+1))*sqrt(S*(S+1)-Σ*(Σ+1))
        elseif δ(Σ,Σ′+1)
            ME += -1.0*sqrt(J*(J+1)-(l+Σ)*(l+Σ-1))*sqrt(S*(S+1)-Σ*(Σ-1))
        end
    end
    return ME 
end


function Rotation(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    """
    Rotational hamiltonian.
    """

    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = 0.0
    if δ(l,l′) && δ(N,N′) && δ(F,F′) && δ(M,M′) && δ(J,J′)
        ME = N*(N+1) - l^2
    end
    return ME
end



function SpinRotation(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    # See Hirota, eq. (2.3.35). There's a simpler form if Λ=0, but this will reduce to it when appropriate.
    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = (-1)^(J + S + N) * (-1)^(N - l) *
        sqrt( S * (S + 1) * (2S + 1) * (2N + 1) * (2N′ + 1) ) *
        wigner6j_(N′, S ,J, S, N, 1) *
        sum( sqrt(2k + 1) *
            (
                (-1)^k * 
                sqrt( N′ * (N′ + 1) * (2N′ + 1) ) * 
                wigner6j_(1, 1, k, N, N′, N′) +
                sqrt( N * (N + 1) * (2N + 1) ) *
                wigner6j_(1, 1, k, N′, N, N)
            ) *
                wigner3j_(N, k, N′, -l, q, l′) * T[q+2, k+1]
            for k in 0:2, q in -1:1
        )

    return  ME * δ(J, J′) * δ(F, F′) * δ(M, M′) * δ(S,S′) * δ(I,I′)
end

function HyperfineFermi(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    """
    Fermi contact hyperfine term.
    """

    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = 0.0
    if δ(l,l′) && δ(N,N′) && δ(F,F′) && δ(M,M′)
        ME = (-1)^(N+S+J′+J+I+F+1) * sqrt((2J′+1)*(2J+1)*S*(S+1)*(2S+1)*I*(I+1)*(2I+1)) *
            wigner6j_(I,J′,F,J,I,1) * wigner6j_(S,J′,N,J,S,1)
    end
    return ME
end

function Hyperfine_Dipolar_c(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    """
    Hyperfine IzSz - (I.S)/3 term. Note this can be defined differently in different papers,
    so check that the right linear combination of constants is being used.
    """

    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME1 = 0.0
    ME2 = 0.0
    if δ(l,l′) && δ(M,M′) && δ(F,F′)
        ME1 = (-1)^(J+I+F+1+N′-l′) * sqrt(30)/2 * sqrt((2J′+1)*(2J+1)*(2N′+1)*(2N+1)) * 
            wigner6j_(I,J′,F,J,I,1) * wigner3j_(N,2,N′,-l′,0,l′) *
            wigner9j_(N′,N,2,S,S,1,J′,J,1)
    end
    if δ(N,N′) && δ(F,F′) && δ(M,M′) && δ(l,l′)
        ME2 = -1/3 * (-1)^(N+S+J′+J+I+F+1) * sqrt((2J′+1)*(2J+1)*S*(S+1)*(2S+1)*I*(I+1)*(2I+1)) * 
            wigner6j_(I,J′,F,J,I,1) * wigner6j_(S,J′,N,J,S,1)
    end
    return ME1+ME2
end

function lDoubling(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    """
    Hyperfine IzSz - (I.S)/3 term. Note this can be defined differently in different papers,
    so check that the right linear combination of constants is being used.
    """

    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = 0.0
    if δ(N,N′) && δ(M,M′) && δ(J,J′) && δ(F,F′)
        Δl = l-l′
        ME = 0.5*N*(N+1) * ( δ(Δl,-2)*(-1)^(-abs(l)) + δ(Δl,2)*(-1)^(abs(l)) )
    end
    return ME
end

function Zeeman(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    # Hirota, equation (2.5.16) and (2.5.19)
    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return (-1)^(F - M) *
        wigner3j_(F, 1, F′, -M, 0, M) *
        (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * 
        wigner6j_(J′, F′, I, F, J, 1) *
        (-1)^(N + S + J + 1) * 
        sqrt( (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1) ) *
        wigner6j_(S, J′, N, J, S, 1) * 
        δ(Λ, Λ′) * δ(N, N′) * δ(M, M′) * δ(l,l′)
end

function Stark(state::LinearCaseB_Bend_Field, state′::LinearCaseB_Bend_Field)
    # Hirota, equation (2.5.35)
    Λ, l, N, S, J, I, F, M = unpack(state)
    Λ′, l′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return -(-1)^(F - M) * 
        wigner3j_(F, 1, F′, -M, 0, M′) * 
        (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) *
        wigner6j_(J, F, I, F′, J′, 1) *
        (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) ) *
        wigner6j_(N, J, S, J′, N′, 1) *
        (-1)^(N - Λ) * sqrt( (2N + 1) * (2N′ + 1) ) *
        wigner3j_(N, 1, N′, -l, 0, l′)
end
