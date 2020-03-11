using Plots
using Printf
using LinearAlgebra
#using PyPlot
#using Seaborn

# -------------------------------------------------------
# User Edits
# -------------------------------------------------------
# Grid Parameters (Convention will be 1 = bottom)
levels    = collect(1000:-100:0)     #
#lowerlev = collect(1000:-100:300)
#upperlev = collect(390:-10:0)
#levels = vcat(lowerlev,upperlev)
δz        = levels[2]-levels[1]*-1  # Cell Height (Constant for now)
kmax      = length(levels)          # maximum index of k
z_f       = ones(Int8,1,kmax)*δz    # Height of each cell
z_c       = ones(Int8,1,kmax+1)*δz  # Height of each cell # Distance between cell midpoints
z_c0      = δz                      # Distance to bottom midpoint
#z_f =hcat(ones(Int8,1,length(lowerlev)-1)*100,ones(Int8,1,length(upperlev)+1)*10)
#z_c = hcat(ones(Int8,1,length(lowerlev))*100,ones(Int8,1,length(upperlev))*10)

# Source/Sink Options --------------------------
z_att     =  400   # Attenuation depth for source
ocn_trns  = 0.43   # Transmitted portion through ocn surface
S0        =  125   # Constant multiplying source term
cp0       = 3850   # J(kg*C)
rho       = 1025   # kg/m^3

# Eddy Diffusivity Options --------------------------
mld       =  300  # Mixed Layer Depth
κ_mld     = 10^-1 # Eddy Diffusivity in mixed-layer
κ_int     = 10^-7 #Eddy Diffusivity in the interior
κ0        = κ_int

# Iteration Parameters ------------------------------
tol       = 1e-12
x_g       = collect(5:5/(kmax-1):10)#ones(kmax)*5
ω         = 1.9
max_iter  = 10379
printint  = 1e6
method    = 3

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
    BC_top  = 1 for Dirichlet, 2 for Neumann
    val_top = Input value for top
    BC_bot  = 1 for Dirichlet, 2 for Neumann
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
    κ,S,                             # κ and Source/Sink Terms
    BC_top,val_top,                  # Bottom BCs
    BC_bot,val_bot,                  # Top BCs
    z_t,z_b)                         # Top and Bottom heights

    # Preallocate
    C = zeros(Float64,3,kmax)
    B = zeros(Float64,kmax)
    A = zeros(Float64,kmax,kmax)

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

    # Output A Matrix (optional)
    du = C[3,1:end-1]
    dl = C[1,2:end]
    d  = C[2,:]
    A  = Tridiagonal(dl,d,du)

    return C,B,A
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

Note that input 8 is temporary (need to figure out correct
way to compute residuals)

Input
    1) C       - coefficients
    2) B       - source term
    3) x_g     - array of initial guesses
    4) tol     - tolerance level for residual
    5) Method  - (1) Jacobi, (2) Gauss-Siedel, (3) SOR
    6) max_iter - maximum number of iterations
    7) ω       - weights for successive overrelaxation
    8) x_inv   - Inverted solution for error calculation
"""

function FD_calc_T(C,B,x_g,tol,Method,max_iter,ω,x_inv,printint,A_in)
    allstart = time()

    # Preallocate Arrays [1 x # of cells]
    Tz  = zeros(Float64,2,length(x_g)) # Array for current profile (itr=m)

    # Input Initial Guess(2nd row of Tz)
    Tz[2,:] = x_g

    r   = Float64[]#zeros(Float64,max_iter+1)

    #Tz_all = zeros(Float64,max_iter,length(x_g))

    # Set up counters and indexing
    itcnt = 0     # Count of Iterations
    m     = 0     # Iteration Storage Number
    err   = tol+1 # initial value, pre error calculation (need to set up for the diff between guess and actual)

    while  err > tol #m <= max_iter #|| err > tol
        start = time()
        #@printf("Starting iter %i...",m)
        r_m  = zeros(Float64,1,kmax+1)
        for k = 1:kmax

            # First Cell
            if k == 1
                T0 = 0
                T2 = Tz[2,k+1] # Take T2 from last iteration (or init guess)

            # Last Cell
            elseif k == kmax
                # For Jacobi Method, use previous values (Row 2 of Tz)
                if Method == 1
                    T0 = Tz[2,k-1]

                # For GS and SOR, use new value (Row 1 of Tz)
                else
                    T0 = Tz[1,k-1]
                end

                T2 = 0

            # Interior Cells
            else
                # For Jacobi Method, use previous values (row 2 Tz)
                if Method == 1
                    T0 = Tz[2,k-1]
                # For GS and SOR, use new value (row 1 Tz)
                else
                    T0 = Tz[1,k-1]
                end

                # Use old T2 value for all cases
                T2 = Tz[2,k+1]
            end

            # Get Coefficients and Source Term
            # C[Coeff#,level]
            C0 = C[1,k]
            C1 = C[2,k]
            C2 = C[3,k]
            B1 = B[k]

            # Compute Temperature
            T1 = 1/C1 * (B1 - C2*T2 - C0*T0)
            r_m[k] = C0*T0+C1*T1+C2*T2 - B1

            # Apply weights for sucessive overrelaxation
            if Method == 3 && itcnt > 0
                global Tz[1,k] = ω * T1 + (1-ω) * Tz[2,k]
            else
                global Tz[1,k] = T1
            end

            #Residuals
            #r_m[k] = (C0*T0 + C1*Tz[m,k] + C2*T2 - B1).^2
            if mod(m,printint)==0
                #@printf("\tIT: %i k=%i with err %e\n",m,k,r_m[k])
            end
        end

        itcnt += 1
        m += 1

        # Cora's Method
        #r   = A_in*Tz[1,:]-B
        #err = norm(r)
        #err = sqrt(sum(r_m.^2)/length(r))

        # Tridiag method
        du = C[3,1:end-1] # Upper diagonal
        dl = C[1,2:end]   # Lower diagonal
        d  = C[2,:]       # Diagonal

        C_tri = Tridiagonal(dl,d,du)

        # Calculate residual
        err = norm(C_tri*Tz[1,:] - B)

        # Compute residual but differencing from last iteration
        #err = maximum(broadcast(sqrt,abs2.((Tz[1,:]-Tz[2,:])./Tz[1,:])))
        #err = sqrt(sum(norm(Tz[1,:]-Tz[2,:])))#./Tz[1,:])))
        #err = norm(Tz[1,:]-Tz[2,:])/length(Tz)
        #r[m] = sum(abs.((Tz[1,:]-Tz[2,:])))#./Tz)))
        #err = maximum(abs.((Tz[1,:]-Tz[2,:])))#./Tz)))

        #r[m] = sum(r_m,dims=2)[1]
        #x_itr = Tz[1,:]
        #x_err = sqrt.((x_itr' .- x_inv).^2)
        #r[m] =sum(x_err,dims=2)[1]
        #err = sum(x_err,dims=2)[1]

        # Push calculated error to residual variable
        push!(r,err)

        # Save new Tz profile to Tz0
        Tz[2,:] = Tz[1,:] # Note that index 2 is the prev

        # Save Tz profile to corresponding index on storage
        #Tz_all[m,:] = Tz[1,:]

        #@printf("Now on iteration %i with residual %.5E \n",m,err)

        elapsed = time() - allstart
        if mod(m,printint)==0
            @printf("Now on iteration %i with resid %.2E in %fs\n",m,err,elapsed)
            #@printf("\tTz ")
            #print(Tz)
            #@printf("\n\tTz0 ")
            #print(Tz0)

            #@printf("Now on iteration %i with err %f\n",m,err)
        end

    end



    # Restrict to iteration
    #r  =  r[1:itcnt]
    #Tz = Tz[1:itcnt,:]
    elapsed = time()-allstart
    @printf("Completed in %i iterations with resid %.2E. Only took %fs\n",itcnt,err,elapsed)
    return Tz, itcnt,r#, Tz_all
end
"""
FD_inv_sol -------------------------------------------
Find the solution of Ax=B by inverting A, which contains
the coefficients of a tridiagonal matrix.
The first dim of a indicates...
    1 = subdiagonal (lower)
    2 = diagonal
    3 = superdiagonal (upper)

Inputs
    1) A       = Coefficients (coeff#,k)
    2) B       = R.H.S. (1xk)

Outputs
    1) A_tri   = Tridiagonal Matrix


"""
function FD_inv_sol(A,B)
    # Create Tridiagonal Matrix
    du = A[3,1:end-1] # Upper diagonal
    dl = A[1,2:end]   # Lower diagonal
    d  = A[2,:]       # Diagonal
    A_tri = Tridiagonal(dl,d,du)
    x = B' * inv(A_tri)
   return A_tri,x
end

#err = B .- Tridiagonal(C[1,2:end],C[2,:],C[3,1:end-1])*Tz'

#Make κ
κ                  = FD_make_κ(mld,levels,kmax,κ_mld,κ_int)

# Calculate I
S                  = FD_calc_I(levels,z_f,z_att,S0,ocn_trns,rho,cp0,mld)

# Calculate Coefficients
C,B_new,A_in       = FD_calc_coeff(kmax,z_f,z_c,κ,S,BC_top,val_top,BC_bot,val_bot,z_t,z_b,z_c0,κ0)

# Get Solution via inverting matrix
C_new,Tz_inv       = FD_inv_sol(C,B_new)

# Calculate T
Tz_new,itcnt,resid = FD_calc_T(C,B_new,x_g,tol,method,max_iter,ω,Tz_inv,printint,A_in)

# ----------------------
# Julia Error
#r  = Ax-b
#err= norm(r)/sqrt(r)    # Sum of the squares

# ----------------------
#
# Plotting Comparison between inverted solution and iterated
#gr()
plot(Tz_inv[end,:],levels,
        title="Temperature Profiles",
        xlabel="Temperature(°C)",
        ylabel="Depth (m)",
        yticks = 0:100:1000,
        yflip=true,
        fmt=png,
        lw=2.5,
        linestyle=:dot,
        linecolor=:black,
        labels="Solution (Inversion)",
        legend=:topleft)
plot!(Tz_new[1,:],levels,
        lw = 1,
        linecolor=:red,
        labels="Solution (Iteration) x" * string(itcnt))#,
#savefig(p,"HW1_Solution.svg")
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


# p = plot(Tz_new',levels,
#         xlabel="Temp",
#         ylabel="Depth",
#         yticks = 0:100:1000,
#         yflip=true,
#         fmt=png)
# ## ---------------------------------------
# ## Convergence test with different methods
# ## ---------------------------------------
# storedecay = Any[]
# storemaxit = Array{Float32}(undef,3)
# for M = 1:3
#     _,tot_it,errdecay = FD_calc_T(C,B_new,x_g,tol,M,max_iter,1.5,Tz_inv)
#     push!(storedecay,errdecay)
#     storemaxit[M] = tot_it
#     @printf("%i iterations for method %i\n",tot_it,M)
# end
#
# plot_conv = plot(1:maximum(storemaxit),storedecay[1],
#     title="Convergence by Iteration Method",
#     xlabel="Number of Iterations",
# #    xlims = [0,4e4],
#     xaxis=:log,
#     ylabel="Residual",
#     lw=2.5,
#     linecolor=:red,
#     label="Jacobi x" * string(storemaxit[1]))
# plot_conv = plot!(storedecay[2],
#     lw=2.5,
#     linecolor=:blue,
#     label="Gauss-Seidel x" * string(storemaxit[2]))
#
# plot_conv = plot!(storedecay[3],
#     lw=2.5,
#     linecolor=:black,
#     label="SOR (w=1.5) x" * string(storemaxit[3]))
# savefig(plot_conv,"HW1_Convergence.svg")
#
#
#
# ## ---------------------------------------
# ## Convergence test with SOR, varying ω
# ## ---------------------------------------
# ωvar = collect(0.1:0.1:1.9)
# storemaxit = Array{Float32}(undef,length(ωvar))
# for w = 1:length(ωvar)
#     wval = ωvar[w]
#     _,tot_it,_ = FD_calc_T(C,B_new,x_g,tol,3,1000000,wval,Tz_inv)
#     storemaxit[w] = tot_it
#     @printf("%i iterations for ω %f\n",tot_it,wval)
# end
#
# plot_w = plot(ωvar,storemaxit,
#     title="Iterations to Convergence vs Omega",
#     xlabel="Omega",
#     ylabel="Iterations to Convergence",
# #    yaxis=:log,
#     lw=2.5,
#     linecolor=:red)
# savefig(plot_w,"HW1_SOR.svg")
#
#
# ## ---------------------------------------
# ## Varying κ
# ## ---------------------------------------
# K_const  = ones(Float64,1,kmax)*10^-2
# K_linear = collect(κ_int:((κ_mld-κ_int)/(kmax-1)):κ_mld)
#
# Kprofs = [K_const,K_linear]
#
# T_byK = Any[]
# maxit_byK=Array{Float32}(undef,3)
# for K = 1:2
#     k_use = Kprofs[K]
#
#     #Calculate Coefficients (varing k)
#     C_vk,B_vk = FD_calc_coeff(kmax,z_f,z_c,k_use,S,BC_top,val_top,BC_bot,val_bot,z_t,z_b)
#
#     # Get Solution via inverting matrix
#     _,Tz_inv_vk = FD_inv_sol(C_vk,B_vk)
#
#     # Calculate T
#     Tz_vk,itcnt,_ = FD_calc_T(C_vk,B_vk,x_g,tol,3,1000000,1.5,Tz_inv_vk)
#
#     # Store K
#
#     push!(T_byK,Tz_vk[end,:])
#
#     # Store iterations
#     maxit_byK[K] = itcnt
# end
#
# # Assign third K
# maxit_byK[3] = itcnt
# push!(T_byK,Tz_new[end,:])
#
# plot_K = plot(T_byK[1],levels,
#     title="Temp Profiles (Varying K)",
#     xlabel="Temperature(°C)",
#     ylabel="Depth (m)",
#     yticks = 0:100:1000,
#     yflip=true,
#     lw=2.5,
#     linecolor=:red,
#     label="Constant x" * string(maxit_byK[1]),
#     legend
#     layout=(1,2))
# plot_K = plot!(T_byK[2],levels,
#     lw=2.5,
#     linecolor=:blue,
#     label="Linear x" * string(maxit_byK[2]),
#     layout=(1,2))
#
# plot_K = plot!(T_byK[3],levels,
#     lw=2.5,
#     linecolor=:black,
#     label="Step x" * string(maxit_byK[3]),
#     layout=(1,2))
#
# #savefig(plot_byK,"HW1_Kvary.svg")
#
#
# plot_K = plot(K_const',levels,
#     title="Temp Profiles (Varying K)",
#     xlabel="Eddy Diffusivity (m2/s)",
#     ylabel="Depth (m)",
#     yticks = 0:100:1000,
#     yflip=true,
#     lw=2.5,
#     linecolor=:red,
#     label="Constant 0.01",
#     legend=:bottomright,
#     layout=(1,2))
# plot_K = plot!(K_linear,levels,
#     lw=2.5,
#     linecolor=:blue,
#     label="Linear 0.01:1",
#     layout=(1,2))
#
# plot_K = plot!(κ[:,1:end-1]',levels,
#     lw=2.5,
#     linecolor=:black,
#     label="Step @ 300 m" * string(maxit_byK[3]),
#     layout=(1,2))
# savefig(plot_K,"HW1_K.svg")
