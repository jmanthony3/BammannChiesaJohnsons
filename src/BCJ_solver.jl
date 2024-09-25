struct BCJ_metal
    θ::Float64                      # temperature
    strain_rate::Float64            # strain rate
    strain_total::Float64           # total strain
    incnum::Int64                   # number of strain increments
    istate::Int64                   # load type (1: uniaxial tension; 2: torsion)
    params::Dict{String, Float64}   # material constants
end

struct BCJ_metal_currentconfiguration
    incnum1::Int64
    μ::Float64                      # shear modulus
    S::Matrix{Float64}              # deviatoric stress tensor
    ϵₚ::Matrix{Float64}             # plastic strain
    ϵₜₒₜₐₗ::Matrix{Float64}         # total strain
    Δϵ::Vector{Float64}             # strain increment
    Sₜᵣ::Vector{Float64}            # deviatoric stress tensor (trial)
    ξ::Vector{Float64}              # overstress tensor (S - 2/3*alpha)
    Δt::Float64                     # timestep
    strain_rate_effective::Float64  # strain rate (effective)
    V::Float64                      # temperature dependent something
    Y::Float64                      # temperature dependent something
    f::Float64                      # temperature dependent something
    rd::Float64                     # dynamic recovery of isotropic hardening
    rs::Float64                     # static recovery of isotropic hardening
    h::Float64                      # isotropic hardening
    Rd::Float64                     # dynamic recovery of kinematic hardening
    Rs::Float64                     # static recovery of kinematic hardening
    H::Float64                      # kinematic hardening
    κ::Vector{Float64}              # isotropic hardening tensor
    αₜᵣ::Vector{Float64}            # kinematic hardening tensor (trial)
    α::Matrix{Float64}              # kinematic hardening tensor
    β::Float64                      # yield function
    Tot::Vector{Float64}            # kappa + alpha + beta
end

function BCJ_metal_currentconfiguration_init(BCJ::BCJ_metal)::BCJ_metal_currentconfiguration
    θ               = BCJ.θ
    strain_rate     = BCJ.strain_rate
    strain_total    = BCJ.strain_total
    incnum          = BCJ.incnum
    istate          = BCJ.istate
    params          = BCJ.params
    incnum1         = incnum+1
    # breakout params into easy variables
    G       = params["bulk_mod"]
    μ       = params["shear_mod"] * 2.
    C1      = params["C01"]
    C2      = params["C02"]
    C3      = params["C03"]
    C4      = params["C04"]
    C5      = params["C05"]
    C6      = params["C06"]
    C7      = params["C07"]
    C8      = params["C08"]
    C9      = params["C09"]
    C10     = params["C10"]
    C11     = params["C11"]
    C12     = params["C12"]
    C13     = params["C13"]
    C14     = params["C14"]
    C15     = params["C15"]
    C16     = params["C16"]
    C17     = params["C17"]
    C18     = params["C18"]
    C19     = params["C19"]
    C20     = params["C20"]


    # array declarations
    # * tenXirs: # = [#_11, #_22, #_33, #_12, #_23, #_13]
    ## OSVs
    S   = zeros(Float64, (6, incnum1))  # deviatoric Stress
    ϵₚ  = zeros(Float64, (6, incnum1))  # plastic Strain
    ϵₜₒₜₐₗ  = zeros(Float64, (6, incnum1))  # total Strain
    ## ISVs
    α= zeros(Float64, (6, incnum1))  # alpha = kinematic hardening
    κ = zeros(Float64,     incnum1 )  # kappa = iXitropic hardening
    Tot = zeros(Float64,     incnum1 )  # kappa + alpha + beta
    ## holding values
    Δϵ  = zeros(Float64,  6)            # strain increment
    Sₜᵣ = zeros(Float64,  6)            # trial stress  (deviatoric)
    αₜᵣ= zeros(Float64,  6)            # trial kinematic
    ξ  = zeros(Float64,  6)            # overstress (S - 2/3*alpha)

    # initialize variables (non-zeros)
    α[:, 1] .= 0.0000001
    S[:, 1]    .= 0.0
    ϵₚ[:, 1]   .= 0.0
    ϵₜₒₜₐₗ[:, 1]   .= 0.0
    κ[1]      = 0.0
    Tot[1]      = 0.0


    # state evaluation - loading type
    if (istate == 1)            # uniaxial tension
        totale_incnum = strain_total / incnum
        Δϵ     .= [totale_incnum, -0.499totale_incnum, -0.499totale_incnum, 0., 0., 0.]
        Δt      = totale_incnum/strain_rate      # timestep
        strain_rate_effective     = strain_rate                    # effective strain rate
    elseif (istate == 2)        # torsion
        # convert equivalent strain to trueshear strain
        strain_total *= 0.5 * √(3.)
        Δϵ     .= [0., 0., 0., strain_total / incnum, 0., 0.]
        # equivalent strain rate to true shear strain rate
        Δt      = Δϵ[3] / strain_rate            # timestep
        strain_rate_effective     = 2strain_rate / √(3.)        # effective strain rate
    end


    # deviatoric strain increment & effective strain rate
    DE_ave = sum(Δϵ[1:3]) / 3.
    # effective strain rate - manually calculated
    # strain_rate_effective  = sqrt((Δϵ[1]^2 + Δϵ[2]^2 + Δϵ[3]^2 \
    #             +(Δϵ[4]^2 + Δϵ[5]^2 + Δϵ[6]^2)*2.)*(2./3.))/Δt


    # temperature dependent constants
    V   = C1 * exp( -C2 / θ )
    Y   = C3 * exp(  C4 / θ )
    f   = C5 * exp( -C6 / θ )

    β= Y + V * asinh( strain_rate_effective / f )
    # β= Y
    # print('β:  ', β)

    rd  = C7 * exp( -C8 / θ )
    h   = C9 -    ( C10 * θ )
    rs  = C11* exp( -C12/ θ )

    Rd  = C13* exp( -C14/ θ )
    H   = C15-    ( C16 * θ )
    Rs  = C17* exp( -C18/ θ )

    Y  *= (C19 < 0.) ? (1.) : (0.5 * ( 1.0 + tanh(max(0., C19 * ( C20 - θ )))))
    return BCJ_metal_currentconfiguration(
        incnum1, μ, S, ϵₚ, ϵₜₒₜₐₗ, Δϵ, Sₜᵣ, ξ, Δt, strain_rate_effective,
        V, Y, f,
        rd, rs, h,
        Rd, Rs, H,
        κ, αₜᵣ, α, β, Tot)
end


"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function solve!(BCJ::BCJ_metal_currentconfiguration)
    incnum1 = BCJ.incnum1
    μ = BCJ.μ
    S = BCJ.S
    ϵₚ, ϵₜₒₜₐₗ, Δϵ = BCJ.ϵₚ, BCJ.ϵₜₒₜₐₗ, BCJ.Δϵ
    Sₜᵣ = BCJ.Sₜᵣ
    ξ = BCJ.ξ
    Δt = BCJ.Δt
    strain_rate_effective = BCJ.strain_rate_effective
    V, Y, f = BCJ.V, BCJ.Y, BCJ.f
    rd, rs, h = BCJ.rd, BCJ.rs, BCJ.h
    Rd, Rs, H = BCJ.Rd, BCJ.Rs, BCJ.H
    κ, αₜᵣ, α, β, Tot = BCJ.κ, BCJ.αₜᵣ, BCJ.α, BCJ.β, BCJ.Tot
    # timestep calculations
    for i ∈ range(2, incnum1)
        Al_mag = sum(α[1:3, i-1] .^ 2.) + 2sum(α[4:6, i-1] .^ 2.)
        # Al_mag = sqrt( Al_mag * 3./2.)       # match cho
        Al_mag = √( Al_mag * 2. / 3.)       # match vumat20


        # trial guesses: ISVs (from recovery) and stress
        ## recovery for alpha
        reco = Δt * (rd*strain_rate_effective + rs) * Al_mag       # original
        # reco = min(reco,1.0)
        αₜᵣ .= α[:, i-1] .* (1 - reco)
        ## recovery for kappa
        Reco = Δt * (Rd*strain_rate_effective + Rs) * κ[i-1]
        # Reco = min(Reco,1.0)
        Katr = κ[i-1] * (1 - Reco)

        ## trial stress guess
        Sₜᵣ    .= S[:, i-1] + μ .* Δϵ           # trial stress
        ξ     .= Sₜᵣ - (2. / 3.) .* αₜᵣ       # trial overstress original
        # ξ     .= Sₜᵣ - sqrt(2. / 3.) .* αₜᵣ   # trial overstress FIT
        Xi_mag  = √(sum(ξ[1:3] .^ 2.) + 2sum(ξ[4:6] .^ 2.))



        # ----------------------------------- #
        ###   ---   YIELD CRITERION   ---   ###
        # ----------------------------------- #
        Crit = Xi_mag - sqrt(2. / 3.)*(Katr + β)         # same as vumat20
        # Crit = Xi_mag - (Katr + β) #changed to FIT
        if Crit <= 0.           # elastic
            # trial guesses are correct
            κ[i]      = Katr
            α[:, i] .= αₜᵣ
            S[:, i]    .= Sₜᵣ
            ϵₚ[:, i]   .= ϵₚ[:, i-1]
            ϵₜₒₜₐₗ[:, i]   .= ϵₜₒₜₐₗ[:, i-1] + Δϵ
        else                    # plastic
            # Radial Return
            γ           = (Crit) / (μ + 2. / 3. *(h+H))     # original
            κ[i]      = Katr + H * sqrt(2. / 3.) * γ  # original
            S[:, i]    .= Sₜᵣ - (μ * γ) .* (ξ ./ Xi_mag)
            α[:, i] .= αₜᵣ + (h * γ) .* (ξ ./ Xi_mag)
            ϵₚ[:, i]   .= ϵₚ[:, i-1] + Δϵ - μ .* (S[:, i] - S[:, i-1])
            ϵₜₒₜₐₗ[:, i]   .= ϵₜₒₜₐₗ[:, i-1] + Δϵ
        end
        Tot[i] = β + α[1, i] + κ[i]
    end
    # return ϵₜₒₜₐₗ, S, α, κ, Tot # return OSVs and ISVs
    return nothing
end