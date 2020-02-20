using Printf
using Plots
# HW 1 Test
ω        = 1  # Relaxation Parameter
max_iter = 10000 # Maximum number of iterations

# Ok how are we going to do this...
# Starting simply?
# -----------------------------------------------
# Problem with a vertical, 1-D, 10-level model
# of equally spaced cubes

# Grid Parameters (Convention will be 1 = bottom)
levels    = collect(1000:-100:0)     #
kmax      = length(levels);          # maximum index of k
z_f       = ones(Int8,1,kmax)*100    # Height of each cell
z_c       = ones(Int8,1,kmax+1)*100  # Height of each cell # Distance between cell midpoints

# Source/Sink Term Options -----------------------------
z_att     =  400 # Attenuation depth for source
ocn_trns  = 0.43 # Transmitted portion through ocn surface
S0        = 6.7 # Constant multiplying source term
cp0       = 3850 # J(kg*C)
rho       = 1025 # kg/m^3

# Setting Boundary Conditions --------------------------
# Currently Hard coded to work with 1-D F_diff, Temperature
# Where F_diff = -κ(ΔT/Δz)
BC_top    = 2    # Top BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)
BC_bot    = 1    # Bot BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)

# Type 1 (Dirichlet) - Specify temperature/value directly
Tb        = 0    # Temperature at bottom
Tt        = 20    # Temperature at top
z_b       = 2    # Distance between midpts of cell k=1 and k=0
z_t       = 2    # Distance between midpts of cell k=n and k=n+1

# Type 2 (Neumann)  - Specify flux/derivative directly
Ft        = S0 /(cp0*rho*z_f[end])  # Flux at the top
Fb        = 0  /(cp0*rho*z_f[end])  # Flux at the bottom



# Initial Guess at the form of the solution ------------
# https://docs.julialang.org/en/v1/manual/metaprogramming/
ex        = :(x^m)
Tz        = ones(Float64,max_iter+1,kmax)# Evaluate the expression with a guess
# Set some counters and their limits
fmax      = kmax+1; # maximum index f fluxes

# Preallocate for temperature --------------------------
#Tz        = zeros(1,length(levels))

#Define Diffusivity
# κ        = ones(Int8,1,kmax+1)*10^-1
mld       = 500
mldout    = length(findall(x->x>mld,levels))
κ         = hcat(ones(Int8,1,mldout+1)*10^-4,ones(Int8,1,kmax-mldout)*10^-1)

# Compute source term
# Compute Source/Sink term for this layer
#S_all   =  ones(Int8,1,kmax+1)#ocn_trns * S0 * exp.(-1 * levels ./ z_att)

S_all    =  (ocn_trns * S0 * exp.(-1 * levels ./ z_att))./(z_f.*cp0.*rho)'
#Preallocate Matrices
B     = ones(Float64,1,kmax)
C     = zeros(Float64,kmax,kmax)

# Loop for each cell, and compute source/sink term
cnt = 1
while cnt < max_iter+1
    for k = 1:kmax
        # Get Level Value
        lev = levels[k]

        # Get source term for that level
        S = S_all[k]

        # Compute Coefficients
        Ck0 = -κ[k]     / (z_f[k] * z_c[k]  )
        Ck2 = -κ[k+1]   / (z_f[k] * z_c[k+1])
        #Ck1 = κ[k+1]/(z_f[k]*z_c[k+1])+κ[k]/(z_f[k]*z_c[k])
        Ck1 = (Ck0 + Ck2) * -1 # Note this might only work in certain cases

        # Special Options for Bottom Boundary - B - B - B - B - B - B - B - B
        if k == 1

            # Remove bottom coefficient
            Ck0 = 0

            # Dirichlet Bottom BC
            if BC_bot == 1
                # Compute Source Term with BC (Prescribed Temp)
                S   = S + Tb * ( (z_b * κ[k]) / (z_f[k] * z_c[k]) )

                # Ck2 remains the same, Ck1 dep. on BCs
                Ck1 = (     κ[k+1]   / (z_f[k] * z_c[k+1]) ) +
                      ( 2 * κ[k]     / (z_f[k] * z_c[k]  ) ) # Prescribed Temp
            # Neumann Bottom Bx
            elseif BC_bot == 2
                # Compute Source Term with BC (Prescribed Flux)
                S   = S + Fb / z_f[k]

                # Ck2 remains the same, Ck1 dep. on BCs
                Ck1 = κ[k+1] / (z_f[k] * z_c[k+1])
            end

        # Special Option for Top Boundary - T - T - T - T - T - T - T - T - T
        elseif k == kmax

            # Remove top coefficient
            Ck2 = 0

            # Dirichlet Top BX
            if BC_top == 1

                # Compute Source Term with BC (Prescribed Temp)
                S   = S + Tt * ( (z_t * κ[k+1]) / (z_f[k] * z_c[k+1]) )

                # Ck0 remains the same
                Ck1 = ( 2 * κ[k+1]   / (z_f[k] * z_c[k+1]) ) +
                      (     κ[k]     / (z_f[k] * z_c[k]  ) ) # Prescribed Temp

            # Neumann Top Bx
            elseif BC_top == 2
                # Compute Source Term with BC (Prescribed Flux)
                S   = S - Ft / z_f[k]

                # Ck0 remains the same
                Ck1 = κ[k] / (z_f[k] * z_c[k])# This depends on the type of BC

            end


        end

        # Find the temp (depending on index)
        if k == kmax
            T1 = 1/Ck1 * (S - Ck0*Tz[cnt,k-1])

            # assign to matrix
            C[k,k-1]   = Ck0
            C[k,k]     = Ck1

        elseif k == 1
            T1 = 1/Ck1 * (S - Ck2*Tz[cnt,k+1])

            # assign to matrix
            C[k,k]   = Ck1
            C[k,k+1] = Ck2

        else

            # Assign to matrix
            global C[k,k-1] = Ck0
            global C[k,k]   = Ck1
            global C[k,k+1] = Ck2
            T1 = 1/Ck1 * (S - Ck0*Tz[cnt,k-1] - Ck2*Tz[cnt,k+1])
        end

        # Store the temp
        global Tz[cnt+1,k] = T1

        #@printf("%i Calculated Temp for level %i (%i) is %f \n",cnt,k,levels[k],T1)

        # Store new values of S
        global B[k] = S

    end
    @printf("Iteration %i Completed\n",cnt)

    global cnt += 1
end

# Compute Residuals (Broadcasting B)
B_ori = B
C_ori = C
r = C*Tz'.-B'
r_all = sum(r.^2,dims=1)
r_all = r_all[2:end]'

# Compute result by inverting matrix
Tz_inv = B * inv(C)
p = plot(Tz_inv',levels,
        xlabel="Temp",
        ylabel="Depth",
        yticks = 0:100:1000,
        yflip=true,
        fmt=png)

plot!(Tz[end,:],levels)

# #Element by Element Check,
#
# r_all2 = zeros(Float64,1,max_iter)
# for m = 1:max_iter
#
#     A = C
#
#     # Skip first index for original guess of T-profile
#     x = Tz[m+1,:]
#     r = sum(A*x-B')
#
#     global r_all2[1,m] = r
# end

# # Make A Plot of the results
# plot(Tz',levels,
#     xlabel="Temp",
#     ylabel="Depth",
#     yticks = 0:100:1000,
#     yflip=true,
#     fmt=png)


# # Make A Plot of the results
# plot(B',levels,
#     xlabel="Source Term",
#     ylabel="Depth",
#     yticks = 0:100:1000,
#     yflip=true)

# Plot the Residuals
# plot([1:1:max_iter],r_all',#broadcast(abs,r_all'),
#     xlabel="Iterations",
#     ylabel="Error")

# So at this point, we have T^1, first guess
# we can continue to iterate within a for loop
# until the residual is satisfactory
