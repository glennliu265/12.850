using Plots
using Printf
using LinearAlgebra
#using Pyplot
#using PyPlot
#using Seaborn

include("AllFunctions12850.jl")

# -------------------------------------------------------
# User Edits
# -------------------------------------------------------
# Grid Parameters (Convention will be 1 = bottom)
levels    = collect(1000:-10:0)
mpts      = ocnmod.get_midpoints(levels)
#lowerlev = collect(1000:-100:300)
#upperlev = collect(390:-10:0)
#levels = vcat(lowerlev,upperlev)
δz        = (levels[2]-levels[1])*-1  # Cell Height (Constant for now)
kmax      = length(mpts)          # maximum index of k
z_f       = ones(Int8,1,kmax)*δz    # Height of each cell
z_c       = ones(Int8,1,kmax)*δz  # Height of each cell # Distance between cell midpoints
z_c0      = δz                      # Distance to bottom midpoint
#z_f =hcat(ones(Int8,1,length(lowerlev)-1)*100,ones(Int8,1,length(upperlev)+1)*10)
#z_c = hcat(ones(Int8,1,length(lowerlev))*100,ones(Int8,1,length(upperlev))*10)

# Source/Sink Options --------------------------
z_att     =  400   # Attenuation depth for source
ocn_trns  = 0.43   # Transmitted portion through ocn surface
S0        =  110   # Constant multiplying source term
cp0       = 3850   # J(kg*C)
rho       = 1025   # kg/m^3

# Eddy Diffusivity Options --------------------------
mld       =  300  # Mixed Layer Depth
κ_mld     = 10^-6 # Eddy Diffusivity in mixed-layer
κ_int     = 10^-6 #Eddy Diffusivity in the interior
κ0        = κ_int

# Iteration Parameters ------------------------------
tol       = 1e-4
ω         = 1.9
max_iter  = 10379
printint  = 1e6
method    = 3

# Setting Boundary Conditions --------------------------
# Currently Hard coded to work with 1-D F_diff, Temperature
# Where F_diff = -κ(ΔT/Δz)
BC_top    = 1    # Top BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)
BC_bot    = 1    # Bot BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)

# Value at top/bottom
val_top   = 10 #S0/(cp0*rho*mld) # For this case, I use a flux
val_bot   = 2 # In this case, I just prescribe a temperature at the boundary

# Distance to the temperature value
z_t       = δz/2
z_b       = δz/2

# Initial Conditions
x_init    = ones(Float64,kmax)*7.5 # Initial Temperature Profile

# Timestepping Options (current model is in seconds)
# Use mean MLD (200)
Δt        = 3600*24*30 #1 * 3600 * 24 * 30     # One Month Resolution
ts_max    = 36 #3 * 12   * Δt      # Years of integration
mld_cycle = sin.((pi/6).*[1:12;]).*100 .+ 200
mld_cycle = ones(Int8,12)*200
θ         = 0.5 # 0 for Forward Euler, 0.5 for CN, 1 for BW Euler



# Plotting Options
plotseas = 0 # Set to 1 to plot seasonal cycle of MLD/Qsw

# Vary Forcing
# Use mean FSNS (110) and max/min of around 20-200
#I_cycle   = sin.((pi/6).*[1:12;].-(pi/2)).*90 .+ 110
I_cycle = sin.((pi/6).*[1:12;].-(pi/2)).*S0
I_cycle = ones(Int8,12)*S0

# -------------------------------------------------------
# Make κ (seasonal)
# -------------------------------------------------------
κ_seas = zeros(Float64,12,kmax) # Preallocate
for imon = 1:12
        κ_seas[imon,:] = ocnmod.FD_make_κ(mld_cycle[imon],mpts,kmax,κ_mld,κ_int)
end

# plot(κ_seas')

# -------------------------------------------------------
# Calculate Q_seasonal [month x profile]
# -------------------------------------------------------
Q_seas  = zeros(Float64,12,kmax) # Preallocate
for imon = 1:12
        Q_seas[imon,:] = ocnmod.FD_calc_I(mpts,z_f,z_att,I_cycle[imon],ocn_trns,rho,cp0,mld_cycle[imon])
end
Q_seas  = zeros(Float64,12,kmax) # Preallocate

# Change surface flux to cycle seasonally
if BC_top == 2
        val_top        = I_cycle ./ (cp0.*rho.*mld_cycle)
end

# Dirichlet BCs
if BC_top == 1
    val_top = ones(Int8,12).*val_top
end

if BC_bot == 1
    val_bot = ones(Int8,12).*val_bot
end

# plot(Q_seas')
#plot([1:12;],val_top[:])

# -------------------------------------------------------
# Calculate seasonal matrices
# -------------------------------------------------------
C     = zeros(Float64,12,3,kmax)
S_new = zeros(Float64,12,kmax)
A_in  = zeros(Float64,12,kmax,kmax)

# Compute matrices for each mmonth
for imon = 1:12
        # Calculate Coefficients
        C[imon,:,:],S_new[imon,:],A_in[imon,:,:] = ocnmod.FD_calc_coeff(kmax,z_f,z_c,κ_seas[imon,:],Q_seas[imon,:],BC_top,val_top[imon],BC_bot,val_bot[imon],z_t,z_b,z_c0,κ0)
end

# -------------------------------------------------------
# Crank Nicolson Method...
# -------------------------------------------------------
#Preallocate
Tprof  = zeros(Float64,ts_max+1,kmax)
Tz_inv = zeros(Float64,ts_max+1,kmax)

#Iteration Information
itall  = zeros(ts_max)

#Time series for debugging
b_ts   = zeros(Float64,ts_max,kmax)
sr_ts =  zeros(Float64,ts_max,kmax)
sl_ts = zeros(Float64,ts_max,kmax)
t0_ts = zeros(Float64,ts_max,kmax)


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

    # if any(x -> x < 0, Sr_in)
    #         @printf("\tLoop %i has negative in Sr\n",i)
    #         #b = b .* -1
    # end
    # if any(x -> x < 0, Sl_in)
    #         @printf("\tLoop %i has negative in Sl\n",i)
    #         #b = b .* -1
    # end
    # if any(x -> x < 0, B_in)
    #         @printf("\tLoop %i has negative in B\n",i)
    #         #b = b .* -1
    # end
    # if any(x -> x < 0, C_in)
    #         @printf("\tLoop %i has negative in C\n",i)
    #         #b = b .* -1
    # end

    # Prepare matrices
    kprint = 1
    A,b,t0 = ocnmod.CN_make_matrix(Δt,θ,Tprof[i,:],C_in,B_in,Sr_in,Sl_in,1,kprint,z_f);

    if any(x -> x < 0, b)
            #@printf("\tLoop %i has negative in b\n",i)
            #b = b .* -1
    end



    # Get Solution via inverting matrix
    ~,Tz_inv[i+1,:]       = ocnmod.FD_inv_sol(A,b)

    # Get Solution via iteration
    itsol,itall[i],resid = ocnmod.FD_itrsolve(A,b,x_init,tol,method,ω,50)
    Tprof[i+1,:] = itsol[1,:]

    # Save timeseries for debugging
    b_ts[i,:] = b
    sr_ts[i,:] = Sr_in
    sl_ts[i,:] = Sl_in
    t0_ts[i,:] = t0

    # if any(x -> x < 0, Tz_inv[i+1,:])
    #         @printf("\tLoop %i has negative in Tz_inv\n",i)
    #         Tz_inv[i+1,:] = Tz_inv[i+1,:] .* -1
    #         Tprof[i+1,:] = Tprof[i+1,:] .* -1
    #
    # end

    elapsed = loopstart - time()
    if i%12==0
        @printf("\nLoop %i finished in %fs\n",i,elapsed)
        #@printf("\tTop Temp is %.2E\n",Tz_inv[i,end])
        #@printf("\tFt       is %.2E\n",b[end])
        #@printf("\n")
    end
end




# Anim Example
anim = @animate for i ∈ 1:ts_max
    l = @layout[a{0.3w} grid(2,1)]
    # Get Month
    m = i%12
    if m == 0
        m=12
    end

    # Ival = val_top[m]#I_cycle[m]
    # Mval = mld_cycle[m]

    p=plot(Tprof[i,:],mpts,
            title="Temperature Profile \nt=" * lpad(string(i),2,"0") * "; Mon=" * lpad(string(m),2,"0"),
            xlabel="Temperature(°C)",
            ylabel="Depth (m)",
            yticks = 0:100:1000,
            yflip=true,
            xlims=(0, 20),
            lw=2.5,
            linestyle=:dot,
            linecolor=:black,
            labels="Solution (Inversion)",
            layout = l,
            subplot = 1,
            legend=:none
            )
    p=plot!(Tprof[i,:],mpts,
            lw = 1,
            linecolor=:red,
            labels="Solution (Iteration) x" * string(itall[i]),
            layout = l,
            subplot=1
            )
    p=plot!([1:12;],val_top,
            subplot=2,
            label    = "Forcing",
            xlabel   = "Months",
            xticks   = 3:3:12,
            ylabel   = "degC/s",
            title    = "Forcing (Top)",
            legend   = :none,
            linecolor= :blue
            )
    p=plot!([m],seriestype=:vline,
            subplot  = 2,
            linecolor= :black,
            linewidth= 2.5
            )
    # p=plot!(m,I_cycle[m],
    #         subplot=2)
    p=plot!([1:12;],mld_cycle,
            subplot  = 3,
            label    = "MLD",
            xlabel   = "Months",
            xticks   = 3:3:12,
            ylabel   = "m",
            title    = "Mixed Layer Depth",
            legend   = :none,
            linecolor=:green
            )
    p=plot!([m],seriestype=:vline,
            subplot=3,
            linecolor=:black,
            linewidth=2.5
            )
    # p=plot!([1:12;],I_cycle,
    #         inset= bbox(0,0,1,1, :bottom, :right),


end
gif(anim,"TempProf_forcinginset.gif",fps=4)
#
# # ------------
# # Contour Plot
# # ------------
# cplot=contourf([0:36;],reverse(mpts),Tz_inv',
#         title   ="Temperature Contours",
#         ylabel  ="Depth (m)",
#         yticks  =([200:200:800;],["800","600","400","200"]),
#         xlabel  ="Months",
#         fillcolor=:thermal)
# #clabel(cplot)
# savefig(cplot,"HW2_Contours.svg")
#
# # ------------
# # Seasonal Plot
# # ------------
# if plotseas == 1
#     splot=plot([1:12;],I_cycle,
#             xlabel   = "Months",
#             xticks   = 3:3:12,
#             grid    = false,
#             ylabel   = "Incoming SWrad (W/m2)",
#             yguidefontsize =12,
#             linewidth = 2.5,
#             linecolor= :blue,
#             title    = "Seasonal Cycle",
#             label   = "Q_sw",
#             )
#     splot=plot!(twinx(),mld_cycle,
#             xticks   = 3:3:12,
#             grid    = false,
#             label    = "MLD",
#             ylabel   = "MLD Depth (m)",
#             linewidth = 2.5,
#             linecolor=:green,
#             legend= :inside
#             )
#     savefig(splot,"HW2_Seasonal_Cycle.svg")
# end
#






# i = 6
#     # Get Month
#     m = i%12
#     if m == 0
#         m=12
#     end
#     l = @layout [a b]#; c{0.2h}]
#     p1=plot(Tz_inv[i,:],mpts,
#             subplot=1,
#             title="Temperature Profile \nt=" * lpad(string(i),2,"0") * "; Mon=" * lpad(string(m),2,"0"),
#             xlabel="Temperature(°C)",
#             ylabel="Depth (m)",
#             yticks = 0:100:1000,
#             yflip=true,
#             xlims=(0, 20),
#             lw=2.5,
#             linestyle=:dot,
#             linecolor=:black,
#             labels="Solution (Inversion)")
#     p1=plot!(Tprof[i,:],mpts,
#             lw = 1,
#             linecolor=:red,
#             labels="Solution (Iteration) x" * string(itall[i]))
#

    # plot([1:12;],I_cycle,            #inset= bbox(.1,.2,0.25,0.25, :bottom, :right),
    # title = "Forcing (SWRad)",
    # xlabel="Month",
    # ylabel="W/m2")
    #
    # scatter!(12,350,
    # markershape=:hexagon,
    # markersize=:100)
    #plot(p1,p2,layout=l)


#
#
# end
#
#
#
# end
#
# #
# # # Get Solution via inverting matrix
# # C_new,Tz_inv       = ocnmod.FD_inv_sol(C,B_new)
#
# # Calculate T
# Tz_new,itcnt,resid = ocnmod.FD_calc_T(kmax,C,B_new,x_g,tol,method,max_iter,ω,printint,A_in)
#
# # Plotting Comparison between inverted solution and iterated
# #gr()
# plot(Tz_inv[end,:],levels,
#         title="Temperature Profiles",
#         xlabel="Temperature(°C)",
#         ylabel="Depth (m)",
#         yticks = 0:100:1000,
#         yflip=true,
#         fmt=png,
#         lw=2.5,
#         linestyle=:dot,
#         linecolor=:black,
#         labels="Solution (Inversion)",
#         legend=:topleft)
# plot!(Tz_new[1,:],levels,
#         lw = 1,
#         linecolor=:red,
#         labels="Solution (Iteration) x" * string(itcnt))#,
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
#
# # Test CN Make Matrix
# IC = x_init
# C  = C_in
# B = B_in
# fr = Sr_in
# fl = Sl_in
# meth = 1
#
#
# # Determine LHS and RHS multipliers
# # LHS - For timestep (n+1), multiply by θ
# l_mult = Δt*(θ)
# # RHS - For timestep (n)  , multiply by 1-θ
# r_mult =  Δt*(1-θ)
#
# # Meth1: Add Timestep corrections first
# if meth == 1
# C[2,:] = C[2,:] .+ (1/r_mult)
# B[2,:] = B[2,:] .- (1/l_mult)
# end
#
# # Multiply variables by time and theta factors
# B      = B  .* l_mult
# fl     = fl .* l_mult
# C      = C  .* r_mult
# fr     = fr .* r_mult
#
# # Meth2: Add single digit post-multiplication
# if meth == 2
# C[2,:] = C[2,:] .+ 1
# B[2,:] = B[2,:] .+ 1
# end
#
# # Now combine terms
# C_tri = ocnmod.makeTridiag(C)
# b     = C_tri * IC + fr + fl # Bring Sl from left side
# A     = B
#
# elapsed = time() - start
