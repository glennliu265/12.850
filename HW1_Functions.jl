using Plots
using Printf
using LinearAlgebra
#using Seaborn

# -------------------------------------------------------
# User Edits
# -------------------------------------------------------
# Grid Parameters (Convention will be 1 = bottom)
levels    = collect(1000:-100:0)     #
#levels    = hcat([1000:100:0])
kmax      = length(levels)         # maximum index of k
z_f       = ones(Int8,1,kmax)*100   # Height of each cell
z_c       = ones(Int8,1,kmax+1)*100  # Height of each cell # Distance between cell midpoints


# Source/Sink Options --------------------------
z_att     =  400  # Attenuation depth for source
ocn_trns  = 0.43   # Transmitted portion through ocn surface
S0        =  125  # Constant multiplying source term
cp0       = 3850   # J(kg*C)
rho       = 1025   # kg/m^3

# Eddy Diffusivity Options --------------------------
mld       =  300  # Mixed Layer Depth
κ_mld     = 10^0 # Eddy Diffusivity in mixed-layer
κ_int     = 10^-2 #Eddy Diffusivity in the interior

# Iteration Parameters ------------------------------
tol       = 1e-8
x_g       = ones(kmax)*5#collect(0:10/(kmax-1):10)
ω         = 0.7
max_iter  = 1000
method    = 1

# Setting Boundary Conditions --------------------------
# Currently Hard coded to work with 1-D F_diff, Temperature
# Where F_diff = -κ(ΔT/Δz)
BC_top    = 2    # Top BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)
BC_bot    = 1   # Bot BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)

val_top   = S0/(cp0*rho*mld) # For this case, I use a flux
val_bot   = 5   # In this case, I just prescribe a temperature at the boundary

z_t       = 2
z_b       = 2

"""
FD_make_κ -------------------------------------------
Create vertical profile of eddy diffusivity (κ) based on
2 layers (mixed layer with higher κ and interior with lower
κ). κ0 is included in the interior layers (extra entry for
κ_mld)

Eventually: Implement "transition" layer rather of
diffusivity values

Inputs
    1) mld     = Mixed Layer Depth
    2) levels  = Depth levels
    3) kmax    = Number of cells
    4) κ_mld   = Eddy Diffusivity in mixed-layer
    5) κ_int   = Eddy Diffusivity in the interior
Outputs
    1) κ       = Vertical diffusivity profile

"""
function FD_make_κ(mld,levels,kmax,κ_mld,κ_int)
    # Find # of levels deeper than mld
    #mldout    = argmin(levels.-mld)
    mldout    = length(findall(x->x>mld,levels))
    # Create κ profile [κ0 as part of interior, κ_int]
    κ         = hcat(ones(Float64,1,mldout + 1)    * κ_int,
                     ones(Float64,1,kmax - mldout) * κ_mld)
   return κ
end
"""
FD_calc_I -------------------------------------------
Calculate the source/sink term for incoming shortwave radiation
Inputs
    1) levels = depth levels
    2) z_f    = height of each cell
    3) z_att  = attentuation depth
    4) S0     = Initial incoming solar radiation
    5) cp0  (Default set to 3850)
    6) rho  (Default set to 1025)
Outputs
    1) S      = Source/Sink Term
"""
function FD_calc_I(levels,z_f,z_att,S0,ocn_trns,rho,cp0,mld)
    #Q = (ocn_trns * S0 * exp.(-1 * levels ./ z_att))
    #S = Q ./ (z_f .* cp0 .* rho)'
    S = (ocn_trns * S0 * exp.(-1 * levels ./ z_att)) ./ (mld .* cp0 .* rho)'
    return S
end
"""
FD_calc_coeff -------------------------------------------
Compute both C (Matrix of Coefficients) and corresponding
forcing terms B
Inputs:
    ||~ Box Geometry and Indices ~||
    kmax   = maximum iteration size
    z_f    = Height of each grid cell
    z_c    = Midpoint distance between cells

    ||~ Initial Profiles         ~||
    κ      = Diffusivity Profile
    S  = Initial forcing/sink term (initial)

    ||~ Boundary Conditions      ~||
    BC_top = 1 for Dirichlet, 2 for Neumann
    val_top = Input value for top
    BC_bot = 1 for Dirichlet, 2 for Neumann
    val_bot = Input value for bot
    z_t     = Height to top cell
    z_b     = Height to bottom cell
Outputs:
    1) C - Matrix of Coefficients [1:3,kmax]
        1 = C(k,k-1); 2 = C(k,k); 3 = C(k,k+1)
    2) B - Column vector of Source/Sink Terms, with
        first and last entries modified for BCs
"""
function FD_calc_coeff(kmax,z_f,z_c, # Cell geometry/indices
    κ,S,                         # κ and Source/Sink Terms
    BC_top,val_top,                  # Bottom BCs
    BC_bot,val_bot,                 # Top BCs
    z_t,z_b)                        # Top and Bottom heights

    # Preallocate
    C = zeros(Float64,3,kmax)
    B = zeros(Float64,kmax)

    # Options for bottom boundary ------------------------------------
    k = 1

    # Dirichlet Bottom
    if BC_bot == 1
        # Compute Source Term with BC (Prescribed Temp)
        B[1]   = S[k] + val_bot *
                 ( (z_t * κ[k+1]) / (z_f[k] * z_c[k+1]) )

        # Ck0 remains the same
        C[2,1] = ( 2 * κ[k+1]   / (z_f[k] * z_c[k+1]) ) +
                 (     κ[k]     / (z_f[k] * z_c[k]  ) ) # Prescribed Temp
    # Newmann Bottom
    elseif BC_bot == 2
        # Compute Source Term with BC (Prescribed Flux)
        global B[1]   = S[k] + val_bot / z_f[k]

        # Ck2 remains the same, Ck1 dep. on BCs
        C[2,1] = κ[k+1] / (z_f[k] * z_c[k+1])
    end

    # Set C(k,k-1) to 0, and set C(k,k+1) to the usual
    C[1,1] = 0
    C[3,1] = -κ[k+1]   / (z_f[k] * z_c[k+1])

    # Options for top boundary ----------------------------------------
    k = kmax

    # Dirichlet Top
    if BC_top == 1

        # Compute Source Term with BC (Prescribed Temp)
        B[kmax]   = S[k] + val_top * ( (z_t * κ[k+1]) / (z_f[k] * z_c[k+1]) )

        # Calculate C(k,k)
        C[2,kmax] = ( 2 * κ[k+1]   / (z_f[k] * z_c[k+1]) ) +
                 (     κ[k]     / (z_f[k] * z_c[k]  ) ) # Prescribed Temp

    # Neumann Top
    elseif BC_top == 2
        # Compute Source Term with BC (Prescribed Flux)
        global B[kmax]   = S[k] - val_top / z_f[k] # CHANGE THIS

        # Calculate C(k,k)
        C[2,kmax] = κ[k] / (z_f[k] * z_c[k])# This depends on the type of BC

    end

    # Set C(k,k+1) to 0, and set C(k,k-1) to the usual
    C[1,kmax] = -κ[k]     / (z_f[k] * z_c[k]  )
    C[3,kmax] = 0

    # Options for interior ----------------------------------------------------
    for k = 2:kmax-1
        B[k] = S[k]
        # Compute Coefficients
        C[1,k] = -κ[k]     / (z_f[k] * z_c[k]  )
        C[3,k] = -κ[k+1]   / (z_f[k] * z_c[k+1])
        C[2,k] = (C[1,k] + C[3,k]) * -1 # Note this might only work in certain cases
    end

    return C,B
end

"""
FD_calc_T -------------------------------------------
For a steady state, 1-D Diffusion Problem
Solves Ax = B using one of the following iterative methods
(1) Jacobi (2) Gauss-Siedel (3) Successive Overrelaxation,
where A is the sparse matrix of coefficients given by C,
x is the temperature profile and B is the source term.

Begins with an initial guess, x and iterates until the residual
(sum of(Ax-B)^2) is below the tolerance level

Input
    1) C       - coefficients
    2) B       - source term
    3) x_g     - array of initial guesses
    4) tol     - tolerance level for residual
    5) Method  - (1) Jacobi, (2) Gauss-Siedel, (3) SOR
    6) max_iter - maximum number of iterations
    7) ω       - weights for successive overrelaxation
"""

function FD_calc_T(C,B,x_g,tol,Method,max_iter,ω)

    #Preallocate [iteration x length]
    Tz = zeros(Float64,max_iter,length(x_g))
    r  = zeros(Float64,max_iter)

    itcnt = 0 # Count of Iterations
    m     = 1 # Iteration Storage Number
    err   = tol+1 # initial value

    while err > tol || m <= max_iter
        r_m  = zeros(Float64,1,kmax)
        for k = 1:kmax

            # Get temperature values
            # Use Initial Guess for first iteration
            if itcnt == 0
                # 1st Cell
                if k == 1
                    T0 = 0
                    T2 = x_g[k+1]
                # Last Cell
                elseif k == kmax
                    T0 = x_g[k-1]
                    T2 = 0
                # Interior Cells
                else
                    T0 = x_g[k-1]
                    T2 = x_g[k+1]
                end
            else

                # First Cell
                if k == 1
                    T0 = 0
                    T2 = Tz[m-1,k+1]
                # Last Cell
                elseif k == kmax
                    # For Jacobi Method, use previous values
                    if method == 1
                        T0 = Tz[m-1,k-1]
                    # For GS and SOR, use new value
                    else
                        T0 = Tz[m  ,k-1]
                    end
                    T2 = 0
                # Interior Cells
                else
                    # For Jacobi Method, use previous values
                    if method == 1
                        T0 = Tz[m-1,k-1]
                    # For GS and SOR, use new value
                    else
                        T0 = Tz[m  ,k-1]
                    end
                    T2 = Tz[m-1,k+1]
                end
            end

            # Get Coefficients and Source Term
            # C[Coeff#,level]
            C0 = C[1,k]
            C1 = C[2,k]
            C2 = C[3,k]
            B1 = B[k]


            # Compute Temperature
            global Tz[m,k] = 1/C1 * (B1 - C2*T2 - C0*T0)

            # Apply weights for sucessive overrelaxation
            if method == 3 && itcnt > 0
                global Tz[m,k] = ω * Tz[m,k] + (1-ω) * Tz[m-1,k]
            end

            r_m[k] = (C0*T0 + C1*Tz[m,k] + C2*T2 - B1).^2
        end

        # Compute residual
        r[m] = sum(r_m,dims=2)[1]
        err = r[m]


        #@printf("Now on iteration %i with residual %.5E \n",m,err)
        if mod(m,5000)==0
            @printf("Now on iteration %i\n",m)
        end
        itcnt += 1
        m += 1
    end

    # Restrict to iteration
    r  =  r[1:itcnt]
    Tz = Tz[1:itcnt,:]

    return Tz, itcnt
end



#Make κ
κ = FD_make_κ(mld,levels,kmax,κ_mld,κ_int)

# Calculate I
S   = FD_calc_I(levels,z_f,z_att,S0,ocn_trns,rho,cp0)

# Calculate Coefficients
C,B_new = FD_calc_coeff(kmax,z_f,z_c,κ,S,BC_top,val_top,BC_bot,val_bot,z_t,z_b)


# Calculate T
Tz_new,itcnt = FD_calc_T(C,B_new,x_g,tol,Method,max_iter,ω)

du = C[3,1:end-1]
dl = C[1,2:end]
d  = C[2,:]
C_new = Tridiagonal(dl,d,du)
Tz_inv = B_new' * inv(C_new)

p = plot(Tz_inv[end,:],levels,
        xlabel="Temp",
        ylabel="Depth",
        yticks = 0:100:1000,
        yflip=true,
        fmt=png)

# plot(S[1:500],levels[1:500],
#         xlabel="Temp",
#         ylabel="Depth",
#         yticks = 0:100:1000,
#         yflip=true,
#         fmt=png)

plot!(Tz_new[end,:],levels)#,
#         xlabel="Temp",
#         ylabel="Depth",
#         yticks = 0:100:1000,
#         yflip=true,
#         fmt=png)
# p = plot(Tz_new[1:100:1000,:]',levels,
#         xlabel="Temp",
#         ylabel="Depth",
#         yticks = 0:100:1000,
#         yflip=true,
#         fmt=png)
