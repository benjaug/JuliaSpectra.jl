using Parameters

export 
    Params_Linear_2Sigma, Params_Linear_2Pi, Params_SymTop_2A1, Params_SymTop_2E,
    MolParams


const sol = 29979.2458 # converts cm-1 to MHz.

# Start a dictionary to hold all molecule choices.
MolParams = Dict{String, Any}()

""" Define custom types to hold the parameters. """ 

@with_kw struct Params_Linear_2Sigma
    T0 = 0.0 # origin
    B = 0.0 # rotational constant
    D = 0.0 # rotational distortion
    γ = 0.0 # spin-rotation
    γD = 0.0 # spin-rotation distortion
    bF = 0.0 # Fermi contact hyperfine
    c = 0.0 # dipolar hyperfine

    gS = 2.0023 # spin g-factor
    gL = 1.0 # Orbital g-factor
    gl′ = 0.0 # parity-dependent g-factor
    gl = 0.0 

    μa = 0.0 # electric dipole moment along symmetry axis
end

@with_kw struct Params_Linear_2Pi
    T0 = 0.0 # origin
    B = 0.0 # rotational constant
    D = 0.0 # rotational distortion
    p2q = 0.0 # Lambda-doubling
    q = 0.0 # Lambda-doubling
    aSO = 0.0 # spin-orbit constant
    a = 0.0 # Orbital hyperfine

    gS = 2.0023 # spin g-factor
    gL = 1.0 # Orbital g-factor
    gl′ = 0.0 # parity-dependent g-factor
    gl = 0.0 

    μa = 0.0 # electric dipole moment along symmetry axis
end

@with_kw struct Params_SymTop_2A1
    T0 = 0.0 # origin
    A = 0.0 # rotational constant A 
    B = 0.0 # rotational constant B = C
    DN = 0.0 # rotational distortion on N 
    DNK = 0.0 # rotational distortion constant 
    DK = 0.0 # rotational distortion
    ϵaa = 0.0 # spin-rotation constant, a axis 
    ϵbc = 0.0 # spin-rotation constant, b ( =c) axis 

    aF = 0.0 # Fermi contact hyperfine
    Taa = 0.0 # dipolar hyperfine tensor component
    TbbmTcc = 0.0 # dipolar hyperfine tensor Tbb - Tcc component.

    μa = 0.0 # dipole moment along a-axis.
    gS = 2.0023 
end 

@with_kw struct Params_SymTop_2E
    T0 = 0.0 # origin
    A = 0.0 # rotational constant A
    B = 0.0 # rotational constant B
    aζed = 0.0 # spin-orbit including Coriolis and quenching
    ϵaa = 0.0 # spin-rotation on a-axis
    ϵbc = 0.0 # spin-rotation on b=c axis.
    Aζt = 0.0 # Coriolis
    ϵ1 = 0.0 # one of the Jahn-Teller parameters (like p in Λ-doubling)
    h1 = 0.0 # another Jahn-Teller parameter (like q in Λ-doubling)

    μa = 0.0 # dipole moment along a-axis.
    gS = 2.0023
end



# Now define the molecules

# CaF 
X = Params_Linear_2Sigma(
    B= 10303.988/sol, # Childs 1981, 10.1016/0022-2852(81)90288-5
    D = 4.969e-7, # Devlin 2015, 10.1016/j.jms.2015.07.009
    γ = 39.65891/sol, # Childs 1981, 10.1016/0022-2852(81)90288-5
    bF = 109.1839/sol, # Childs 1981, 10.1016/0022-2852(81)90288-5
    c = 40.1190/sol, # Childs 1981, 10.1016/0022-2852(81)90288-5
    μa = 3.07 # Ernst 1989, 10.1103/PhysRevA.39.1575
    )

A = Params_Linear_2Pi(
    T0 = 16529.10177, # Devlin 2015, 10.1016/j.jms.2015.07.009
    B = 0.347395, # Devlin 2015, 10.1016/j.jms.2015.07.009
    aSO = 72.61743, # Devlin 2015, 10.1016/j.jms.2015.07.009
    p2q = -0.0452, # Devlin 2015, 10.1016/j.jms.2015.07.009
    q = -0.0002916, # Devlin 2015, 10.1016/j.jms.2015.07.009
    gl′ = -0.0611, # Devlin 2015, 10.1016/j.jms.2015.07.009
    μa = 2.45 # Ernst 1989, 10.1103/PhysRevA.39.1575
    )

B = Params_Linear_2Sigma(
    B = 0.341259, # Devlin 2015, 10.1016/j.jms.2015.07.009
    γ = -0.045994, # Devlin 2015, 10.1016/j.jms.2015.07.009
    T0 = 18833.12751, # Devlin 2015, 10.1016/j.jms.2015.07.009
    gl = 0.0696, # Devlin 2015, 10.1016/j.jms.2015.07.009
    gS = 1.9977, # Devlin 2015, 10.1016/j.jms.2015.07.009
    μa = 2.05 # Raouafi 2001, 10.1063/1.1405118
    )

CaF_Params = Dict("X" => X, "A" => A, "B" => B)
MolParams["CaF"] = CaF_Params

# SrF 
X = Params_Linear_2Sigma(
    B= 7487.60/sol, # Barry thesis
    D = 0.0075/sol, # Barry thesis
    γ = 75.02249/sol, # Barry thesis
    bF = 97.6670/sol, # Barry thesis
    c = 29.846/sol, # Barry thesis
    μa = 3.4963 # Barry thesis
    )

A = Params_Linear_2Pi(
    T0 = 15216.34287, # Barry thesis
    B = 0.2528335, # Barry thesis
    aSO = 281.46138, # Barry thesis
    p2q = -0.133002, # Barry thesis
    q = 0.0, # 
    gl′ = 0.0, # 
    μa = 2.064 # Barry thesis
    )

B = Params_Linear_2Sigma(
    B = 0.24961, # Barry thesis
    γ = -0.134, # Barry thesis
    T0 = 17267.4465, # Barry thesis
    gl = 0.0, # Barry thesis
    gS = 2.002, # Barry thesis
    μa = 0.91 # Barry thesis
    )

SrF_Params = Dict("X" => X, "A" => A, "B" => B)
MolParams["SrF"] = SrF_Params

# YbOH
X = Params_Linear_2Sigma(
    B = 0.245116257, # Steimle 2019, 10.1103/PhysRevA.100.052509
    D = 2.029e-7, # Steimle 2019, 10.1103/PhysRevA.100.052509
    γ = -0.002707, # Steimle 2019, 10.1103/PhysRevA.100.052509
    γD = 1.59e-7, # Steimle 2019, 10.1103/PhysRevA.100.052509
    μa = 1.9, # Steimle 2019, 10.1103/PhysRevA.100.052509
    gl = 0.0055, # Steimle 2019, 10.1103/PhysRevA.100.052509
    bF = 4.8/sol, # Nakhate 2018, 10.1016/j.cplett.2018.11.030
    c = 2.46/sol # Nakhate 2018, 10.1016/j.cplett.2018.11.030
    )

A = Params_Linear_2Pi(
    aSO = 1350.0, # Steimle 2019, 10.1103/PhysRevA.100.052509
    B = 0.253052, # Steimle 2019, 10.1103/PhysRevA.100.052509
    p2q = -0.43807, # Steimle 2019, 10.1103/PhysRevA.100.052509
    T0 = 17998.5875, # Steimle 2019, 10.1103/PhysRevA.100.052509
    μa = 0.43, # Steimle 2019, 10.1103/PhysRevA.100.052509
    gl′ = -0.865 # Steimle 2019, 10.1103/PhysRevA.100.052509
    )

YbOH_Params = Dict("X" => X, "A" => A)
MolParams["YbOH"] = YbOH_Params

# SrOH
X = Params_Linear_2Sigma(
    B = 0.249199814, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    D = 2.17437e-7, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    γ = 2.42748e-3, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    bF = 1.713/sol, # Fletcher 1993, 10.1063/1.464218
    c = 1.673/sol, # Fletcher 1993, 10.1063/1.464218
    μa = 1.9 # Steimle 1992, 10.1063/1.462007
    )

A = Params_Linear_2Pi(
    B = 0.2537833, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    aSO = 263.58741, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    p2q = -0.143662, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    q = -1.528e-4, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    T0 = 14674.04063, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    gl′= -0.269, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    μa = 0.59 # Steimle 1992, 10.1063/1.462007
    )

B = Params_Linear_2Sigma(
    B = 0.2522066, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    γ = -0.142583, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    T0 = 16377.49826, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    gl = 0.306, # Nguyen 2018, 10.1016/j.jms.2018.02.007
    μa = 0.396 # Steimle 1992, 10.1063/1.462007
    )

SrOH_Params = Dict("X" => X, "A" => A, "B" => B)
MolParams["SrOH"] = SrOH_Params

# CaOH
X = Params_Linear_2Sigma(
    B = 0.334354, # Steimle 1992, 10.1063/1.462007
    γ = 34.7593/sol, # Scurlock 1993, 10.1006/jmsp.1993.1133
    bF = 2.602/sol, # Scurlock 1993, 10.1006/jmsp.1993.1133
    c = 2.053/sol, # Scurlock 1993, 10.1006/jmsp.1993.1133
    μa = 1.465 # Steimle 1992, 10.1063/1.462007
    )

A = Params_Linear_2Pi(
    B = 0.341428, # Steimle 1992, 10.1063/1.462007
    aSO = 66.795, # Steimle 1992, 10.1063/1.462007
    p2q = -0.042796 + 2*(-3.569e-4), # Steimle 1992, 10.1063/1.462007
    q = -3.569e-4, # Steimle 1992, 10.1063/1.462007
    T0 = 15998.128, # Li 1992, 10.1063/1.463322
    μa = 0.836, # Steimle 1992, 10.1063/1.462007
    gl′ = -0.063 # Augenbraun 2022, 10.1103/PhysRevLett.127.263002
    )

B = Params_Linear_2Sigma(
    B = 0.339385, # Steimle 1992, 10.1063/1.462007
    γ = -0.043418, # Steimle 1992, 10.1063/1.462007
    T0 = 18022.268, # Bernath 1985, 10.1086/162800
    μa = 0.744 # Steimle 1992, 10.1063/1.462007
    )

CaOH_Params = Dict("X" => X, "A" => A, "B" => B)
MolParams["CaOH"] = CaOH_Params


# CaCH3
X = Params_SymTop_2A1(
    A = 5.44831, # Marr 1996, 10.1063/1.472265
    B = 0.25238487, # Marr 1996, 10.1063/1.472265
    DN = 3.544852e-7, # Marr 1996, 10.1063/1.472265
    DNK = 1.995314e-5, # Marr 1996, 10.1063/1.472265
    DK = 7.03e-5, # Marr 1996, 10.1063/1.472265
    ϵbc = 0.00185105, # Marr 1996, 10.1063/1.472265
    μa = 2.62 # Marr 1996, 10.1063/1.472265
)

A = Params_SymTop_2E(
    T0 = 14743.3822, # Marr 1996, 10.1063/1.472265
    A = 5.3855, # Marr 1996, 10.1063/1.472265
    B = 0.254270, # Marr 1996, 10.1063/1.472265
    aζed = 72.7092, # Marr 1996, 10.1063/1.472265
    ϵaa = 0.0129, # Marr 1996, 10.1063/1.472265
    ϵbc = 0.0196, # Marr 1996, 10.1063/1.472265
    Aζt = 5.360, # Marr 1996, 10.1063/1.472265
    ϵ1 = -0.0251, # Marr 1996, 10.1063/1.472265
    h1 = 1.1e-4, # Marr 1996, 10.1063/1.472265
    μa = 1.69  # Marr 1996, 10.1063/1.472265
)

CaCH3_Params = Dict("X" => X, "A" => A)
MolParams["CaCH3"] = CaCH3_Params

# CaOCH3
X = Params_SymTop_2A1(
    A = 5.448303, # Crozet 2002, 10.1006/jmsp.2002.8536
    B = 0.11626491, # Crozet 2002, 10.1006/jmsp.2002.8536
    DN = 2.6885e-8, # Crozet 2002, 10.1006/jmsp.2002.8536
    DNK = 2.358e-6, # Crozet 2002, 10.1006/jmsp.2002.8536
    DK = 7.03e-5, # Crozet 2002, 10.1006/jmsp.2002.8536
    ϵbc = 4.15304e-4, # Crozet 2002, 10.1006/jmsp.2002.8536
    μa = 1.58, # Namiki 1998, 10.1063/1.477146

    aF = -0.421/sol, # Namiki 1998, 10.1063/1.477146
    Taa = 1.070/sol, # Namiki 1998, 10.1063/1.477146
    TbbmTcc = 0.292/sol # Namiki 1998, 10.1063/1.477146
)

A = Params_SymTop_2E(
    T0 = 15925.1232, # Crozet 2002, 10.1006/jmsp.2002.8536
    A = 5.43997, # Crozet 2002, 10.1006/jmsp.2002.8536
    B = 0.1178840, # Crozet 2002, 10.1006/jmsp.2002.8536
    aζed = 66.97448, # Crozet 2002, 10.1006/jmsp.2002.8536
    ϵaa = 3.58e-3, # Crozet 2002, 10.1006/jmsp.2002.8536
    ϵbc = 3.20e-3, # Crozet 2002, 10.1006/jmsp.2002.8536
    Aζt = 5.43730, # Crozet 2002, 10.1006/jmsp.2002.8536
    ϵ1 = -8.208e-3, # Crozet 2002, 10.1006/jmsp.2002.8536
    h1 = 1.50e-4, # Crozet 2002, 10.1006/jmsp.2002.8536
    μa = 1.0 # NOT MEASURED?? Assuming around 1 Debye here... 
)

CaOCH3_Params = Dict("X" => X, "A" => A)
MolParams["CaOCH3"] = CaOCH3_Params

