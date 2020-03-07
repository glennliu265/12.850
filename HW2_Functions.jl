using Plots
using Printf
using LinearAlgebra
#using PyPlot
#using Seaborn
include("AllFunctions12850.jl")
# -------------------------------------------------------
# User Edits
# -------------------------------------------------------
# Grid Parameters (Convention will be 1 = bottom)
levels    = collect(1000:-100:0)     #
#lowerlev = collect(1000:-100:300)
#upperlev = collect(390:-10:0)
#levels = vcat(lowerlev,upperlev)
kmax      = length(levels)         # maximum index of k
z_f       = ones(Int8,1,kmax)*1   # Height of each cell
z_c       = ones(Int8,1,kmax+1)*1  # Height of each cell # Distance between cell midpoints
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
κ_mld     = 10^-6 # Eddy Diffusivity in mixed-layer
κ_int     = 10^-6 #Eddy Diffusivity in the interior

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

x_init    = ones(Float64,kmax)*5 # Initial Temperature Profile

# Timestepping Options (current model is in seconds)
# Use mean MLD (200)
Δt     = 3600*24#1 * 3600 * 24 * 30     # One Month Resolution
ts_max = 30#3 * 12   * Δt      # Years of integration
#mld_cycle   = sin.((pi/6).*[1:12;]).*100 .+ 200
mld_cycle = ones(Int8,12)*200
θ = 0.5

# Vary Forcing
# Use mean FSNS (110) and max/min of around 20-200
#I_cycle   = sin.((pi/6).*[1:12;].-(pi/2)).*90 .+ 110
I_cycle = sin.((pi/6).*[1:12;].-(pi/2)).*90
# -------------------------------------------------------
# Make κ (seasonal)
# -------------------------------------------------------
κ_seas = zeros(Float64,12,kmax+1) # Preallocate
for imon = 1:12
        κ_seas[imon,:] = ocnmod.FD_make_κ(mld_cycle[imon],levels,kmax,κ_mld,κ_int)
end
# plot(κ_seas')

# -------------------------------------------------------
# Calculate Q_seasonal [month x profile]
# -------------------------------------------------------
Q_seas  = zeros(Float64,12,kmax) # Preallocate
val_top = zeros(Float64,12)
for imon = 1:12
        Q_seas[imon,:] = ocnmod.FD_calc_I(levels,z_f,z_att,I_cycle[imon],ocn_trns,rho,cp0,mld_cycle[imon])
        val_top        = I_cycle ./ (cp0.*rho.*mld_cycle)
end
# plot(Q_seas')

# -------------------------------------------------------
# Calculate seasonal matrices
# -------------------------------------------------------
C     = zeros(Float64,12,3,kmax)
S_new = zeros(Float64,12,kmax)
A_in  = zeros(Float64,12,kmax,kmax)
# Compute matrices for each mmonth
for imon = 1:12
        # Calculate Coefficients
        C[imon,:,:],S_new[imon,:],A_in[imon,:,:] = ocnmod.FD_calc_coeff(kmax,z_f,z_c,κ_seas[imon,:],Q_seas[imon,:],BC_top,val_top[imon],BC_bot,val_bot,z_t,z_b)
end

# -------------------------------------------------------
# Crank Nicolson Method...
# -------------------------------------------------------
#Preallocate
Tprof  = zeros(Float64,ts_max+1,kmax)
Tz_inv = zeros(Float64,ts_max+1,kmax)

Tprof[1,:] = x_init
Tz_inv[1,:] = x_init
for i = 1:ts_max
    loopstart = time()
    m = i%12

    # Get Next Month (for eqn L.H.S.)
    nm = m + 1

    # Set December month to 12
    if m == 0
        m = 12
    end

    # Get Corresponding Matrices
    C_in  = C[m,:,:]
    Sr_in = S_new[m,:,:]
    Sl_in = S_new[nm,:,:]
    B_in  = C[nm,:,:]

    # Prepare matrices
    A,b = ocnmod.CN_make_matrix(Δt,θ,Tprof[i,:],C_in,B_in,Sr_in,Sl_in,2);

    # Get Solution via inverting matrix
    ~,Tz_inv[i+1,:]       = ocnmod.FD_inv_sol(A,b)
    Tprof[i+1,:] = Tz_inv[i+1,:]
    elapsed = loopstart - time()
    if i%1==0
        @printf("Loop %i finished in %fs\n",i,elapsed)
    end
end



anim = @animate for i ∈ 1:ts_max
    p=plot(Tz_inv[i,:],levels,
            title="Temperature Profiles t="*string(i),
            xlabel="Temperature(°C)",
            ylabel="Depth (m)",
            yticks = 0:100:1000,
            yflip=true,
            xlims=(-10, 30),
            lw=2.5,
            linestyle=:dot,
            linecolor=:black,
            labels="Solution (Inversion)")
end
gif(anim,"TempProfanim_m2.gif",fps=2)

i=83456
p=plot(Tz_inv[i ,:],levels,
        title="Temperature Profiles t="*string(i),
        xlabel="Temperature(°C)",
        ylabel="Depth (m)",
        yticks = 0:100:1000,
        yflip=true,
        xlims=(-10, 30),
        lw=2.5,
        linestyle=:dot,
        linecolor=:black,
        labels="Solution (Inversion)")

    # # Iterate to find solution
    # Tprof[m,:],itcnt,resid = ocnmod.CN_solver(kmax,A,Sr_in,Sl_in,B_in,
    #                                             x_g,x_init,printint,
    #                                             Method,ω,θ)

                                                # Calculate T
  # Tprof[m,:],itcnt,resid = ocnmod.FD_calc_T(A,b,x_g,tol,method,max_iter,ω,Tz_inv,printint,A_in)
  #





end



end

#
# # Get Solution via inverting matrix
# C_new,Tz_inv       = ocnmod.FD_inv_sol(C,B_new)

# Calculate T
Tz_new,itcnt,resid = ocnmod.FD_calc_T(kmax,C,B_new,x_g,tol,method,max_iter,ω,printint,A_in)

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
