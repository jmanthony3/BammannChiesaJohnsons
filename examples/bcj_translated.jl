###############################################################################
# Daniel Kenney
# Summer 2023 

# Point simulator for calibrating BCJ model to experimental data
# Consistent with vumat20.f - january 1995 version (Bammann - Chisea Model)

#   ***   Some formula changes may be changed to match EPIC results   ***
# This model does not incorporate damage, or the damage parameter m
# Damage must be calibrated by simulating triaxiality
#
###############################################################################
###############################################################################
###############################################################################
# Outline of code:
#
# Declare model constant
# Declare arrrays
# Evaluate loading state
# Calculate temperature dependent variables
# Loop timesteps:
#   - trial guess
#   - yield critereon
#       - elastic = trial guess is correct
#       - plastic = radial return correction
# 

###############################################################################


# One step through
"""
Function to get a full stress-strain curve (and ISV values)

params = material constants

istate: 1 = tension, 2 = torsion

**no damage in this model**
"""
function BCJ(params, Temp, strainrate, totale, incnum, istate)
    incnum1 = incnum+1
    # Breakout params into easy variables
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
    
    bulk_mod= params["bulk_mod"]
    shr_mod = params["shear_mod"] *2. #/5.           # *2 to match vumat20, /5. to get fit to vumat20


# --------------------------------------
# Array declarations:
# --------------------------------------
#   tenXirs: # = [#_11, #_22, #_33, #_12, #_23, #_13]

    # OSVs
    S   = zeros((6,incnum1)   , float)    # Deviatoric Stress
    PE  = zeros((6,incnum1)   , float)    # Plastic Strain
    TE  = zeros((6,incnum1)   , float)    # Total Strain

    # ISVs
    Alph= zeros((6,incnum1)   , float)    #alpha = kinematic hardening
    Kap = zeros( incnum1      , float)    #kappa = iXitropic hardening
    Tot = zeros( incnum1      , float)    #kappa + alpha + beta

    # Holding values
    DE  = zeros(        6    , float)    # strain increment
    Str = zeros(        6    , float)    # trial stress  (deviatoric)
    Altr= zeros(        6    , float)    # trial kinematic
    Xi  = zeros(        6    , float)    # overstress (S - 2/3*alpha)

# --------------------------------------
# Initialize variables (non-zeros)
    for k in range(6)
        Alph[k][1]  = 0.0000001
        S[k][1]     = 0.0
        PE[k][1]    = 0.0
        TE[k][1]    = 0.0
    end
    Kap[1]      = 0.0
    Tot[1]      = 0.0


### ====================================================================== ####
#         State evaluation - Loading type
    if (istate == 1)           # Uniaxial tension:
        DE[1]   = totale/incnum
        DE[2]   = -0.499*DE[1]
        DE[3]   = -0.499*DE[1]
        DE[4]   = 0.
        DE[5]   = 0.
        DE[6]   = 0.
        dt      = DE[1]/strainrate                  # Timestep
        ddd     = strainrate                        # Effective strain rate
    elseif (istate == 2)         # Torsion
        # Convert equivalent strain to trueshear strain
        totale  = totale * 0.5 * sqrt(3.)
        DE[1]   = 0.
        DE[2]   = 0.
        DE[3]   = 0.
        DE[4]   = totale/incnum
        DE[5]   = 0.
        DE[6]   = 0.
        # equivalent strain rate to true shear strain rate
        dt      = DE[3]/strainrate                  # Timestep
        ddd     = strainrate*2. / sqrt(3.)         # Effective strain rate
    end

    # Deviatoric strain increment & effective strain rate
    DE_ave = (DE[1]+DE[2]+DE[3])/3.
    # Effective strain rate - manually calculated
    # ddd  = sqrt((DE[1]^2 + DE[2]^2 + DE[3]^2 \
    #             +(DE[4]^2 + DE[5]^2 + DE[6]^2)*2.)*(2./3.))/dt




# --------------------------------------
# Temperature dependent constants
    V   = C1 * exp( -C2 / Temp )
    Y   = C3 * exp(  C4 / Temp )
    f   = C5 * exp( -C6 / Temp )

    Beta= Y + V*asinh(ddd/f)
    # Beta= Y
    # print('Beta:  ', Beta)

    rd  = C7 * exp( -C8 / Temp )
    h   = C9  - (C10*Temp)
    rs  = C11* exp( -C12/ Temp )

    Rd  = C13* exp( -C14/ Temp )
    H   = C15 - (C16*Temp)    
    Rs  = C17* exp( -C18/ Temp )

    adj = 1.0
    if (C19 > 0.0)
      adj = 0.5*( 1.0 + tanh(max(0,C19*(C20-Temp))))
    end
    Y = adj * Y


### ====================================================================== ####
### ===================     Timestep Calculations     ==================== ####
### ====================================================================== ####


    for i in range(2,incnum1)
        Al_mag = Alph[1][i-1]^2 + Alph[2][i-1]^2 + Alph[4][i-1]^2 + (
            Alph[4][i-1]^2 + Alph[5][i-1]^2 + Alph[6][i-1]^2)*2.
        # Al_mag = sqrt( Al_mag * 3./2.)       # match cho
        Al_mag = sqrt( Al_mag * 2. / 3.)       # match vumat20


        # Trial guesses: ISVs (from recovery) and stress

        #Recovery for alpha:
        reco = dt *(rd*ddd + rs)* Al_mag       # original
        # reco = min(reco,1.0)
        for k in range(6)
            Altr[k] = Alph[k][i-1] - reco*Alph[k][i-1]
        end

        #Recovery for kappa:
        Reco = dt *(Rd*ddd + Rs)* Kap[i-1]
        # Reco = min(Reco,1.0)
        Katr = Kap[i-1]  - Reco*Kap[i-1]

        # Trial stress guess:
        for k in range(6)
            Str[k] = S[k][i-1] + shr_mod*DE[k]        #Trial Stress
            Xi[k] = Str[k] - (2./3.)*Altr[k]            #trial overstress original
            # Xi[k] = Str[k] - sqrt(2. / 3.)*Altr[k]            #trial overstress FIT
        end
        Xi_mag2 = Xi[0]^2 + Xi[1]^2 + Xi[2]^2 + (Xi[3]^2 + Xi[4]^2 + Xi[5]^2)*2.0
        Xi_mag = sqrt(Xi_mag2)



        # ----------------------------------- #
        ###   ---   YIELD CRITERION   ---   ###
        # ----------------------------------- #
        # print('Xi_mag   :' ,Xi_mag)
        # print('Beta     :' , Beta)

        Crit = Xi_mag - sqrt(2. / 3.)*(Katr + Beta)    #same as vumat20
        # Crit = Xi_mag - (Katr + Beta) #changed to FIT
        if Crit <= 0.0
            # print('ELASTIC')
        # -----------------------------------
        # ELASTIC
        # -----------------------------------
            # Trial guesses are correct
            Kap[i] = Katr
            for k in range(6)
                Alph[k][i]  = Altr[k]
                S[k][i]     = Str[k]
                PE[k][i]    = PE[k][i-1]
                TE[k][i]    = TE[k][i-1] + DE[k]
            end


        else
            # print('PLASTIC')
        # -----------------------------------
        # PLASTIC
        # -----------------------------------
        #   Radial Return 

            Gamma = (Crit) / (shr_mod + 2. / 3. *(h+H))            # original
            Kap[i]          = Katr + H*sqrt(2. / 3.)*Gamma     # original

            for k in range(6)
                S[k][i]     = Str[k] - shr_mod*Gamma*Xi[k] / Xi_mag
                Alph[k][i]  = Altr[k] + h*Gamma*Xi[k] / Xi_mag
                TE[k][i]    = TE[k][i-1] + DE[k]
                PE[k][i]    = PE[k][i-1] + DE[k] - shr_mod*(S[k][i]-S[k][i-1])
            end
        end
        Tot[i] = Beta + Alph[1][i] + Kap[i]
    end
    # ----------------------------------- 

    # for i in range(len(S)):          #Artificial Scaling for troubleshooting
    #     for k in range(len(S[i])):
    #         S[i][k] = sqrt(3/2)*S[i][k]


    # return OSVs and ISVs

    return [TE, S, Alph, Kap, Tot]
end