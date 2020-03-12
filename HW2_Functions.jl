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
S0        = 2000   # Constant multiplying source term
cp0       = 3850   # J(kg*C)
rho       = 1025   # kg/m^3

# Eddy Diffusivity Options --------------------------
mld       =  300  # Mixed Layer Depth
κ_mld     = 10^-1  # Eddy Diffusivity in mixed-layer
κ_int     = 10^-4  # Eddy Diffusivity in the interior
κ0        = κ_int

# Iteration Parameters ------------------------------
tol       = 1e-4
ω         = 1.9
printint  = 1e6
method    = 3

# Setting Boundary Conditions --------------------------
# Currently Hard coded to work with 1-D F_diff, Temperature
# Where F_diff = κ(ΔT/Δz)
BC_top    = 2   # Top BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)
BC_bot    = 1    # Bot BC type (1 = Dirichlet, 2 = Neumann, 3 = Robin)

# Value at top/bottom
#val_top   = S0/(cp0*rho*z_att) # For this case, I use a flux
val_top   = S0/(cp0*rho*z_att)
val_bot   = 3.5# In this case, I just prescribe a temperature at the boundary

# Simple interpolation, keep this value here
z_t       = 2#δz/2
z_b       = 2#δz/2

# Initial Conditions
x_init    = ones(Float64,kmax)*5 # Initial Temperature Profile

# Timestepping Options (current model is in seconds)
# Use mean MLD (200)
Δt        = 3600*24*30 #1 * 3600 * 24 * 30     # One Month Resolution
ts_max    = 36 #3 * 12   * Δt      # Years of integration
mld_cycle = sin.((pi/6).*[1:12;]).*100 .+ 200
#mld_cycle = ones(Int8,12)*200
θ         = 0.5 # 0 for Forward Euler, 0.5 for CN, 1 for BW Euler



# Plotting Options
plotseas = 0 # Set to 1 to plot seasonal cycle of MLD/Qsw

# Vary Forcing
# Use mean FSNS (110) and max/min of around 20-200
#I_cycle   = sin.((pi/6).*[1:12;].-(pi/2)).*90 .+ 110
I_cycle = sin.((pi/6).*[1:12;].-(pi/2)).*S0
#I_cycle = ones(Int8,12)*S0

# -------------------------------------------------------
# Make κ (seasonal)
# -------------------------------------------------------
κ_seas = zeros(Float64,12,kmax) # Preallocate
for imon = 1:12
        κ_seas[imon,:] = ocnmod.FD_make_κ(mld_cycle[imon],mpts,kmax,κ_mld,κ_int)
end

# plot(κ_seas')

# -------------------------------------------------------
# Calculate Source Term, varying by season [month x profile]
# -------------------------------------------------------
Q_seas  = zeros(Float64,12,kmax) # Preallocate
for imon = 1:12
        Q_seas[imon,:] = ocnmod.FD_calc_I(mpts,z_f,z_att,I_cycle[imon],ocn_trns,rho,cp0)
end

# Change surface flux to cycle seasonally
if BC_top == 2
        val_top  = I_cycle ./ (cp0.*rho.*z_att)
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

tot_time = time()
for i = 1:ts_max
    loopstart = time()
    m = i%12

    # Get Next Month (for eqn L.H.S.)
    nm = m + 1

    # Set December month to 12
    if m == 0
        m = 12
    end

    #Get Corresponding Matrices
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
    A,b,t0 = ocnmod.CN_make_matrix(Δt,θ,Tprof[i,:],C_in,B_in,Sr_in,Sl_in,2,kprint,z_f);


    # Get Solution via inverting matrix
    ~,Tz_inv[i+1,:]       = ocnmod.FD_inv_sol(A,b)

    # Get Solution via iteration
    itsol,itall[i],resid = ocnmod.FD_itrsolve(A,b,x_init,tol,method,ω,50)
    Tprof[i+1,:] = itsol[1,:]

    # Save timeseries for debugging
    # b_ts[i,:] = b
    # sr_ts[i,:] = Sr_in
    # sl_ts[i,:] = Sl_in
    # t0_ts[i,:] = t0

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
elapsed = time() - tot_time
@printf("\nAll Loops finished in %fs\n",elapsed)


# ------------
# Anim Example with forcing panels
# ------------
anim = @animate for i ∈ 1:ts_max
    l = @layout[a{0.3w} grid(2,1)]
    # Get Month
    m = i%12
    if m == 0
        m=12
    end

    # Ival = val_top[m]#I_cycle[m]
    #Mval = convert(Int16,mld_cycle[m])

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
            legend=:none,
            linewidth=2.5
            )
    p=plot!(Tprof[i,:],mpts,
            lw = 1,
            linecolor=:red,
            labels="Solution (Iteration) x" * string(itall[i]),
            layout = l,
            subplot=1,
            linewidth=2.5
            )
    p=hline!([mld_cycle[m]],
            layout=l,
            subplot=1,
            linewidth=1.0,
            linecolor=:black,
            line=:dash
            )
    p=plot!([1:12;],val_top,
            subplot=2,
            xlabel   = "Months",
            xticks   = 3:3:12,
            ylabel   = "degC/s",
            title    = "Forcing",
            legend   = :none,
            linecolor= :blue,
            linewidth=2.5
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
            xlabel   = "Months",
            xticks   = 3:3:12,
            ylabel   = "meters",
            title    = "Mixed Layer Depth",
            legend   = :none,
            linecolor=:green,
            linewidth=2.5
            )

    p=plot!([m],seriestype=:vline,
            subplot=3,
            linecolor=:black,
            linewidth=2.5
            )

end
gif(anim,"TempProf_forcinginset.gif",fps=5)

# ------------------
# Forward Euler Test
# ------------------


#
#
# # ------------
# # Anim (just profile)
# # ------------
# anim = @animate for i ∈ 1:ts_max
#     # Get Month
#     m = i%12
#     if m == 0
#         m=12
#     end
#
#     # Ival = val_top[m]#I_cycle[m]
#     # Mval = mld_cycle[m]
#
#     p=plot(Tprof[i,:],mpts,
#             title="Temperature Profile \nt=" * lpad(string(i),2,"0") * "; Mon=" * lpad(string(m),2,"0"),
#             xlabel="Temperature(°C)",
#             ylabel="Depth (m)",
#             yticks = 0:100:1000,
#             yflip=true,
#             xlims=(-20, 20),
#             lw=2.5,
#             linestyle=:dot,
#             linecolor=:black,
#             labels="Solution (Inversion)",
#             )
#     p=plot!(Tprof[i,:],mpts,
#             lw = 1,
#             linecolor=:red,
#             labels="Solution (Iteration) x" * string(itall[i]),
#             )
# end
# gif(anim,"TempProf_profonly.gif",fps=4)
# #

# ------------
# Contour Plot
# ------------
# Repeat mld cycle for whole year

mldplot = repeat(mld_cycle.-1000,convert(Int8,ts_max/12),1)
mldplot = mldplot .*-1
cplot=contourf([0:36;],reverse(mpts),Tz_inv',
        title   ="Temperature Contours",
        ylabel  ="Depth (m)",
        yticks  =([200:200:800;],["800","600","400","200"]),
        xlabel  ="Months",
        clims   =(2, 8),
        fillcolor=:thermal)
cplot=plot!(mldplot,
            linecolor=:white,
            legend=:none,
            line=:dot,
            linewidth=2.5)

#clabel(cplot)
savefig(cplot,"HW2_Contours.svg")


# --------------------
# Diff plot, upper 500m
# --------------------
# cplotdiff=contourf([0:36;],reverse(mpts),Tprofdiff',
#         title   ="Temperature Anomaly (MLDvar-noMLD)",
#         ylabel  ="Depth (m)",
#         ylims=(500, 1000),
#         yticks  =([200:200:800;],["800","600","400","200"]),
#         xlabel  ="Months",
#         fillcolor=:balance)
# savefig(cplotdiff,"HW2_MLDnoMLDdiff.svg")

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

#Tprof_FE_dz = zeros(Float64,3,ts_max+1,kmax)



# ----
# ForEu Testing deltaT
# -----

# # Create array to fill
# Tprof_FE = zeros(Float64,3,ts_max+1,kmax)
# # Fill in array after running
# Tprof_Fe[1,:,:] = Tprof;
# save("Tprof_FE.jld","data",Tprof_FE)
#
# t = 37
# l = @layout[a b]
# pe4 = Tprof_FE[1,37,:];
# pe5 = Tprof_FE[2,37,:];
# pe6 = Tprof_FE[3,37,:];
# Feu = plot(pe4,mpts,
#         title="Temperature Profile (Forward Euler)",
#         xlabel="Temperature(°C)",
#         ylabel="Depth (m)",
#         yticks = 0:100:1000,
#         yflip=true,
#         xlims=(0, 20),
#         lw=2.5,
#         labels="t = 5e4",
#         layout=l,
#         subplot=1,
#         )
# Feu = plot!(pe5,mpts,
#         labels="t = 5e5",
#         lw=2.5,
#         subplot=1)
# Feu = plot!(pe6,mpts,
#         yticks = 0:100:1000,
#         xlabel="Temperature(°C)",
#         ylabel="Depth (m)",
#         yflip=true,
#         xlims=(-5e45, 5e45),
#         lw=2.5,
#         labels="t = 5e6",
#         layout=l,
#         subplot = 2,
#         legend=:left,
#         )
# savefig(Feu,"HW2_ForwardEuler_deltaT.svg")
#



# ----
# ForEu Testing deltaz
# -----
#
# l = @layout[a b]
# pe4 = Tprof100[37,:];
# pe5 = Tprof10[37,:];
# pe6 = Tprof31[37,:];
# Feu = plot(pe4,mp100,
#         title="Temperature Profile (Forward Euler)",
#         xlabel="Temperature(°C)",
#         ylabel="Depth (m)",
#         yticks = 0:100:1000,
#         yflip=true,
#         xlims=(0, 20),
#         lw=2.5,
#         labels="dz = 100m",
#         layout=l,
#         subplot=1,
#         )
# Fedz = plot!(pe5,mp10,
#         labels="dz = 10m",
#         lw=2.5,
#         subplot=1)
# Fedz = plot!(pe6,mp31,
#         yticks = 0:100:1000,
#         xlabel="Temperature(°C)",
#         ylabel="Depth (m)",
#         yflip=true,
#         xlims=(-20, 20),
#         lw=2.5,
#         labels="dz = 3.1m",
#         layout=l,
#         subplot = 2,
#         legend=:left,
#         )
# savefig(Fedz,"HW2_ForwardEuler_deltaZ.svg")
