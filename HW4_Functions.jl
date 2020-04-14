using Plots
using Printf
using LinearAlgebra
using PyPlot

include("AllFunctions12850.jl")

## Grid Set-up  -----------------------------------------------
# X and Y Grids
xgrid = [0:1:50;]
ygrid = [0:1:150;]

# Get Midpoints
mx      = ocnmod.get_midpoints(xgrid)
my      = ocnmod.get_midpoints(ygrid)

# Get Cell Width
δx        = (xgrid[2]-xgrid[1])
δy        = (ygrid[2]-ygrid[1])

# Get End value and indices
xmax        = length(mx)
ymax        = length(my)
Lx          = mx[xmax]
Ly          = my[ymax]

# Get Cell Spacing and Widths (Currently Fixed)
x_f       = ones(Int8,1,xmax)*δx   # Cell Face Length (x)
x_c       = ones(Int8,1,xmax)*δx   # Cell Midpoint Distance (x)

y_f       = ones(Int8,1,ymax)*δy   # Cell Face Length (x)
y_c       = ones(Int8,1,ymax)*δy

# Set distance to 0 cell midpoint (from cell 1)
x_c0      = δx
y_c0      = δy

## Set Diffusivity Parameter -------------------------------
# For vorticity
κx        = ones(Float64,1,xmax) .* 1
κy        = ones(Float64,1,ymax) .* 1
κx0       = κx[1]
κy0       = κy[1]

# For streamfunction, set diffusivity to 1
κxp        = ones(Float64,1,xmax)
κyp        = ones(Float64,1,ymax)
κxp0       = κxp[1]
κyp0       = κyp[1]

## Iteration Parameters ------------------------------
tol       = 1e-7
ω         = 1.6
printint  = 1e4
method    = 3
save_iter = 1000
max_iter  = 1e5

## Time parameters
dt        = 3600*24*30 #1 * 3600 * 24 * 30     # One Month Resolution
ts_max    = 12
θ         = 0.5

## Source Term Settings -------------------------------
ρ0 = 1025 #kg/m3
H  = 1000 #m
τ0 = 0.1 #.* ones(xmax,ymax)  #N/m2

# Set the wind directions
# Zonal Wind

# ------------------------------------------------------------
# Testing different approximation schemes for first derivative
# ------------------------------------------------------------
    # τx = [ -τ0*cos(pi*y/Ly) for y in my]
    # τxy = [ τ0*pi/Ly*sin(pi*y/Ly) for y in my]
    # a1,b1,c=ocnmod.ddx_1d(τx,y_f,1)
    # a2,b2,c=ocnmod.ddx_1d(τx,y_f,2)
    # a3,b3,c=ocnmod.ddx_1d(τx,y_f,3)
    # plot(b1,a1,label="Fw")
    # plot!(b2,a2,label="Bw")
    # plot!(b3,a3,label="Cen")
    # plot!(τxy,label="Ana")
    # It seems that the 3rd method works the best (but lose the ends)


 # Switch from westward to eastward
# contourf(mx,my,τx',
#         levels=10,
#         clabels=true,
#         title="Eastward Wind Stress (N/m2)",
#         xlabel="Zonal Distance (m)",
#         ylabel="Meridional Distance (m)")
#
# txplot = convert(Array{Float32},τx ./ findmax(τx)[1])

# Meridional Wind
#τy = [ -τ0*cos(2*pi*x/Lx)*(y/y) for x in mx, y in my] # Switch from westward to eastward
# τy = zeros(xmax,ymax)
# contourf(mx,my,τy',
#         levels=10,
#         clabels=true,
#         title="Northward Wind Stress (N/m2)",
#         xlabel="Zonal Distance (m)",
#         ylabel="Meridional Distance (m)")
# include("AllFunctions12850.jl")
#dτx_dy = ocnmod.ddx_1d(τx,y_f)


# Approximate the x and y derivatives of the wind directions
# Note: Just repeat the first elemtent for now
#dτx_dy = ocnmod.ddx_1d(vcat(τx[1],τx),y_f)

# Create Zonal Wind Stress Terms
xedge = [0:1:xmax+1;]
yedge = [0:1:ymax+1;]

# Case 1: Constal Zonal Wind Only, varying cosinusoidally
dτx = [ -τ0*cos(2*pi*y/Ly) for x in xedge, y in yedge]
dτy = [ 0 for x in xedge, y in yedge]

x_f2       = ones(Int8,1,length(xedge))*δx
y_f2       = ones(Int8,1,length(yedge))*δx
Txdy,x1,y1 = ocnmod.ddx_2d(dτx,x_f2,xedge,yedge,1)
Tydx,x1,y1 = ocnmod.ddx_2d(dτy,y_f2,xedge,yedge,2)

xi = 5
yi = 10
qscale = 1
pts,uv = ocnmod.quiverprep_2d(xedge,yedge,dτx,dτy,xi,yi,qscale)
p=contourf(xedge,yedge,dτx')
p=Plots.quiver!(pts,quiver=(uv),
    ylims=(0,150),
    xlims=(0,50),
    lc="black"
    )

pts,uv = ocnmod.quiverprep_2d(mx,my,Txdy,Tydx,xi,yi,qscale)
p=contourf(mx,my,Txdy')
p=Plots.quiver!(pts,quiver=(uv),
    ylims=(0,150),
    xlims=(0,50),
    lc="black"
    )



# -----------------------------------
# Directly Prescribe Wind Stress Curl
# -----------------------------------
casen = 1

# Case 01: Forcing Maxima that translates northwest
if casen == 1
    curlTx = [ τ0*pi/Ly*sin(1*(pi*y/Ly+m/6))*x/x for x in mx, y in my, m in 1:12]
    curlTy = [ τ0*pi/Lx*sin(1*(pi*x/Lx+m/6))*y/y for x in mx, y in my, m in 1:12]
end


# Case 02: Cosinusoidal zonal wind, no meriodional
if casen == 2
    curlTx = [ τ0*pi/Ly*sin(1*(pi*y/Ly))*sin(0.5*m/12*pi) for x in mx, y in my, m in 1:12]
    curlTy = [ 0 for x in mx, y in my, m in 1:12]
end


# -----------------------------------
# CALCULATE FORCING TERM
# -----------------------------------
    S = 1/(ρ0*H) * (curlTy - curlTx)


# -----------------------------------
#  Visualize Forcing Term, Quiver plot on wind curl
# -----------------------------------
    casen = 2
    xi = 2
    yi = 10
    qscale = 1e3

    aniquiv2 = @animate for t ∈ 1:12
        u = curlTx[:,:,t]
        v = curlTy[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(mx,my,u,v,xi,yi,qscale)
        Splot = S[:,:,t] ./ findmax(abs.(S[:,:,t]))[1]
        a = heatmap(mx,my,Splot',
#            clims=(-1,1),
            title="Forcing Term Case "*string(casen)*"; t= "*string(t),
            seriescolor=:balance,
            xlabel="x (meters)",
            ylabel="y (meters)",
            )
        a = Plots.quiver!(pts,quiver=(uv),
            ylims=(0,150),
            xlims=(0,50),
            lc="black"
            )
    end
    gif(aniquiv2,"HW4_Forcing_Case"*string(casen)*".gif",fps=2)

## Boundary Conditions
# 1 = "Dirichlet", 2 = "Neumann", 3="Periodic"
# West
# 1 = Dirichlet, 2 = Neumann, 3 = Periodic

WBC = 1
wb_val = [y/y*0 for y in my]#[y/y for y in mx]

# Eastern
EBC = 1
eb_val = [y/y*0 for y in my]#[y/y for y in my]

# North
NBC = 3
nb_val = [x/x*0 for x in mx] #[sin(3*(x/Lx*pi)) for x in mx]

# South
SBC = 3
sb_val = [x/x*0 for x in mx]#[sin(3*(x/Lx*pi)) for x in mx]


## Run the script
# Set periodic boundary options DO NOT EDIT
#___________________________________________
        chk_per = zeros(Int,4)
        wper =0
        eper =0
        sper =0
        nper =0
        if WBC == 3
            wper = 1
        end
        if EBC == 3
            eper = 1
        end
        if NBC == 3
            nper = 1
        end
        if SBC == 3
            sper = 1
        end
        chk_per[1] = nper
        chk_per[2] = sper
        chk_per[3] = eper
        chk_per[4] = wper


# Preallocate
vort_t = zeros(Float64,xmax,ymax,ts_max)
ψ_t = zeros(Float64,xmax,ymax,ts_max)
ut = zeros(Float64,xmax-2,ymax-2,ts_max)
vt = zeros(Float64,xmax-2,ymax-2,ts_max)


for t = 1:ts_max

        if t > 1
                x_init = vort_t[:,:,t-1]
        else
                x_init = zeros(Float64,xmax,ymax)
        end
        loopstart = time()
        m = t%12

        # Get Next Month (for eqn L.H.S.)
        nm = m + 1

        # Set December month to 12
        if m == 0
            m = 12
        end

        # Get Corresponding forcing term
        S0 = S[:,:,m]
        S1 = S[:,:,nm]

        # Compute Coefficients and modifications to BCs
        Cx0,Bx0 = ocnmod.FD_calc_coeff_2D(xmax,x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
        Cy0,By0 = ocnmod.FD_calc_coeff_2D(ymax,y_f,y_c,κy,NBC,nb_val,SBC,sb_val,y_c0,κy0)
        Cx1,Bx1 = ocnmod.FD_calc_coeff_2D(xmax,x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
        Cy1,By1 = ocnmod.FD_calc_coeff_2D(ymax,y_f,y_c,κy,NBC,nb_val,SBC,sb_val,y_c0,κy0)

        # Modify Boundary Conditions
        S0[1,:]    -= Bx0[1,:] # South BC
        S0[xmax,:] -= Bx0[2,:] # North BC
        S0[:,1]    -= By0[1,:] # West BC
        S0[:,ymax] -= By0[2,:] # East BC
        S1[1,:]    -= Bx1[1,:] # South BC
        S1[xmax,:] -= Bx1[2,:] # North BC
        S1[:,1]    -= By1[1,:] # West BC
        S1[:,ymax] -= By1[2,:] # East BC

        # Set up Crank Nicolson
        ug = zeros(Float64,xmax,ymax)
        Ax,Ay,b = ocnmod.CN_make_matrix_2d(dt,θ,ug,Cx0,Cy0,Cx1,Cy1,S0,S1,1e5,chk_per)

        #u_out,itcnt,r,u_scrap,err_scrap,errmap = ocnmod.FD_itrsolve_2D(Cx,Cy,S,ug,tol,ω,method,wper,eper,sper,nper,max_iter,save_iter)
        vort_t[:,:,t],itcnt2,r2 = ocnmod.cgk_2d(Ax,Ay,b,ug,chk_per,tol,max_iter)

        # Solving for streamfunction (first pass)
        # S_psi = vort_t[:,:,t]
        # Cxp,Bxp = ocnmod.FD_calc_coeff_2D(xmax,x_f,x_c,κxp,EBC,eb_val,WBC,wb_val,x_c0,κxp0)
        # Cyp,Byp = ocnmod.FD_calc_coeff_2D(ymax,y_f,y_c,κyp,NBC,nb_val,SBC,sb_val,y_c0,κyp0)
        # S_psi[1,:]    -= Bxp[1,:] # South BC
        # S_psi[xmax,:] -= Bxp[2,:] # North BC
        # S_psi[:,1]    -= Byp[1,:] # West BC
        # S_psi[:,ymax] -= Byp[2,:] # East BC
        # ψ_t[:,:,t],~,~ = ocnmod.cgk_2d(Cxp,Cyp,S_psi,ug,chk_per,tol,max_iter)

        # Solve for streamfunction, using vorticity as the input (functionized)
        S_psi = vort_t[:,:,t]
        ψ_t[:,:,t],~,~ =ocnmod.invPV_2d(S_psi,
                          x_f,y_f,x_c,y_c,x_c0,y_c0,
                          κxp,κyp,κxp0,κyp0,
                          NBC,nb_val,SBC,sb_val,EBC,eb_val,WBC,wb_val,
                          ug,tol,max_iter)
        # Solve for uv
        ut[:,:,t],vt[:,:,t],~,~=ocnmod.psi2uv(ψ_t[:,:,t],x_f,y_f,mx,my,1)


end

# ---------------------------------------------------
# Animate Vorticity (Left) and Streamfunction (Right)
# ---------------------------------------------------
# Prepare for animating vorticity and streamfunction

xi = 2
yi = 10
if casen == 1
    qscale = 1e3
elseif casen == 2
    qscale = 1e6
end
x_u = mx[2:end-1]
y_v = my[2:end-1]
anim3 = @animate for t ∈ 1:ts_max
        l = @layout[a b]

        # Make quivers
        u = ut[:,:,t]
        v = vt[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(x_u,y_v,u,v,xi,yi,qscale)

        c = contourf(mx,my,vort_t[:,:,t]'/findmax(abs.(vort_t))[1],
                clabels=true,
                levels=[-1:0.1:1;],
                clims=(-1,1),
                title="Vorticity at t="*string(t),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )

        h = contourf(mx,my,ψ_t[:,:,t]'/findmax(abs.(ψ_t))[1],
                clabels=true,
                levels=[-1:0.1:1;],
                clims=(-1,1),
                title="Psi at t="*string(t),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
        h = Plots.quiver!(pts,quiver=(uv),
            ylims=(0,150),
            xlims=(0,50),
            lc="black"
            )

        plot(c,h,layout = l)
end
gif(anim3,"HW4_Vort_Psi_Case"*string(casen)*".gif",fps=5)


# -----------------------
# Psi only Plots
# -----------------------
xi = 2
yi = 10
if casen == 1
    qscale = 1e3
elseif casen == 2
    qscale = 1e6
end
x_u = mx[2:end-1]
y_v = my[2:end-1]
anim3 = @animate for t ∈ 1:ts_max
        # Make quivers
        u = ut[:,:,t]
        v = vt[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(x_u,y_v,u,v,xi,yi,qscale)

        h = contourf(mx,my,ψ_t[:,:,t]'/findmax(abs.(ψ_t))[1],
                clabels=true,
                levels=[-1:0.1:1;],
                clims=(-1,1),
                title="Psi at t="*string(t),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
        h = Plots.quiver!(pts,quiver=(uv),
            ylims=(0,150),
            xlims=(0,50),
            lc="black"
            )
end
gif(anim3,"HW4_Psionly_Case"*string(casen)*".gif",fps=5)
#
# anim4 = @animate for t ∈ 1:ts_max
#         contourf(mx,my,ψ_t[:,:,t]'./findmax(ψ_t)[1],
#                 clabels=true,
# #                clims=(-1e-3,1e-3),
#                 title="Streamfunction at t="*string(t),
#                 xlabel="x (meters)",
#                 ylabel="y (meters)",
#                 levels=10,
#                 fillcolor=:inferno
#                 )
#         end
# gif(anim4,"HW4_psifirsttest.gif",fps=10)

xi = 2
yi = 10
qscale = 1e3
aniquiv2 = @animate for t ∈ 1:12
    u = curlTx[:,:,t]
    v = curlTy[:,:,t]
    pts,uv = ocnmod.quiverprep_2d(mx,my,u,v,xi,yi,qscale)
    a = heatmap(mx,my,Splot[:,:,t]',
        clims=(-1,1),
        title="Forcing Term Case "*string(casen)*"; t= "*string(t),
        seriescolor=:balance,
        xlabel="x (meters)",
        ylabel="y (meters)",
        )
    a = Plots.quiver!(pts,quiver=(uv),
        ylims=(0,150),
        xlims=(0,50),
        lc="black"
        )
end
gif(aniquiv2,"HW4_Forcing_Case"*string(casen)*".gif",fps=2)

## F4 - Error Map, EW Periodic, 1e6 iteration
fig5= contourf(mx,my,u_out',
#        clims = (-.1,.1),
        clabels=true,
        #levels=20,
        title="SORStaggered Pt Vortices (No Flow E/W Periodic N/S), t="*string(t),
        xlabel="x (meters)",
        ylabel="y (meters)",
        seriescolor=:balance
        #seriescolor=:inferno
        )


diffme = u_out-uo2;
fig5= contourf(mx,my,diffme',
#        clims = (-.1,.1),
        clabels=true,
        levels=20,
        title="SORStaggered Pt Vortices (No Flow E/W Periodic N/S), It#"*string(itcnt),
        xlabel="x (meters)",
        ylabel="y (meters)",
        seriescolor=:balance
        #seriescolor=:inferno
        )
savefig(fig5,"HW3_F9_StaggeredVortices.svg")




## Heatmap (Vorticity) + Streamfunction Overlay
f4 = plot(
        contourf(mx,my,ζ',
                title="Vorticity",
                fillcolor=:blues,
                clabels=true,
#                xlabel="x (meters)",
#                ylabel="y (meters)",
                levels=10,
#                clims=(-2,2),
                ),
        contourf(mx,my,u_out',
                clabels=true,
                seriescolor=:black,
                title="Psi (No Flow N/S/E/W ), It#"*string(itcnt),
#                xlabel="x (meters)",
#                ylabel="y (meters)",
                levels=10,
                fillcolor=:inferno,
                )
        )
savefig(f4,"HW3_QuadSin2.svg")


# Fig 4 (ζ=10, No flow E/W)
fig4= contourf(mx,my,u_out',
        clabels=true,
        title="Psi (No Flow N/S Periodic E/W; Vorticity=10), It#"*string(itcnt),
        xlabel="x (meters)",
        ylabel="y (meters)",
        levels=10,
        fillcolor=:inferno
        )
savefig(fig3,"HW3_NoFlowNS_PerEW_ConstVort.svg")

# Animate Convergence
anim = @animate for i ∈ 1:save_iter
        c = contourf(mx,my,u_scrap[i,:,:]',
                title="Streamfunction: Iteration "*string(i),
                clabels=true,
                levels=20,
                fillcolor=:balance
                )
end
gif(anim,"HW3_Convergence.gif",fps=10)

# Animate Convergence
anim2 = @animate for i ∈ 1:save_iter
        c = heatmap(mx,my,abs.(err_scrap[i,:,:]'),
                title="Error: Iteration "*string(i),
                clims=(0, 1)
                )
end
gif(anim2,"HW3_Err.gif",fps=10)



anim3 = @animate for i ∈ 1:save_iter
        l = @layout[a b]
        c = contourf(mx,my,u_scrap[i,:,:]',
                title="Streamfunction: Iteration "*string(i),
                clabels=true,
                levels=20,
                fillcolor=:balance,
                )
        h = heatmap(mx,my,abs.(err_scrap[i,:,:]'),
                title="Error: Iteration "*string(i),
                )
        plot(c,h,layout = l)

end
gif(anim3,"HW3_Together_sinusoidal_vert.gif",fps=10)

# # Plot Outside (Since BCs are wonky)
# contourf(mx[2:end-1],my,u_out[2:end-1,:])
# contourf(mx[2:end-1],my,u_out[2:end-1,:],
#         xlims=(1, 5),
#         ylims=(1, 5),
#         )

# For Small Grids
# Currently it seems that z must by (y,x)
# heatmap(mx,my,u_out)
# heatmap(mx[2:end-1],my,u_out[:,2:end-1])
#
#
# xgr= mx[2:end-1]
# ygr= my
# ugr = u_out[2:end-1,:]
# heatmap(ygr,xgr,ugr)


## Fig 1 (ζ=10, No flow)
# fig1 = contourf(mx,my,u_out',
#         clabels=true,
#         title="Psi (No Flow N/S/E/W; Vorticity=10), It#"*string(itcnt),
#         xlabel="x (meters)",
#         ylabel="y (meters)",
#         levels=10,
#         fillcolor=:inferno
#         )
# savefig(fig1,"HW3_NoFlow_ConstVort.svg")


## Fig 2 (ζ=10, No flow N/S Periodic E/W)
# fig2= contourf(mx,my,u_out',
#         clabels=true,
#         title="Psi (No Flow N/S Periodic E/W; Vorticity=10), It#"*string(itcnt),
#         xlabel="x (meters)",
#         ylabel="y (meters)",
#         levels=10,
#         fillcolor=:inferno
#         )
# savefig(fig2,"HW3_NoFlowNS_PerEW_ConstVort.svg")

## Periodic Convergence Test
# pct = plot(r2,
#     label = "N/S Periodic",
#     title = "Iterations to Convergence for Different BCs",
#     xlabel="Iterations to Convergence",
#     ylabel="Residual"
#     )
# pct = plot!(r3,
#     label = "E/W Periodic"
#     )
# savefig(pct, "PeriodicConvergence.svg")

## Convergence Test
# par = [0.5:.1:1.9;]
# itall = zeros(Float64,length(par))
# for i = 1:length(par)
#         ω = par[i]
#         method = 3
#         ~,itall[i],~,~, = ocnmod.FD_itrsolve_2D(Cx,Cy,S,ug,tol,ω,method,wper,eper,sper,nper,1e5,save_iter)
# end
#
# ct = plot(par,itall,
#     title="SOR Convergence Test",
#     xlabel="Omega",
#     xlim = (0.5,2),
#     ylabel="Iterations to Convergence",
#     labels="Iterations",
#     legend=:topleft,
#     linewidth=2.5
#     )
# savefig(ct,"HW3_SORConvergence.svg")

## Fig4: Sinusoidal Vorticity, No Flow N/S/E/W
# f4 = plot(
#         contourf(mx,my,ζ',
#                 title="Vorticity",
#                 fillcolor=:blues,
#                 clabels=true,
# #                xlabel="x (meters)",
# #                ylabel="y (meters)",
#                 levels=10,
#                 clims=(0.25,2),
#                 ),
#         contourf(mx,my,u_out',
#                 clabels=true,
#                 seriescolor=:black,
#                 title="Psi (No Flow N/S/E/W ), It#"*string(itcnt),
# #                xlabel="x (meters)",
# #                ylabel="y (meters)",
#                 levels=10,
#                 fillcolor=:inferno,
#                 )
#         )
# savefig(f4,"HW3_sinVort_NoFlowNSEW.svg")

##

##F2 - Plot Convergence Error
# l = @layout[a; b]
#
# p1 = plot(r_1,
#         xaxis = :log,
# #        title="Convergence Speed for Iterative Methods",
#         title="No Flow N-S-E-W",
#         xlabel="Iterations to Convergence",
#         ylabel="Residual",
#         lw = 2.5,
#         label="Jacobi",
#         )
# p1 = plot!(r_2,
#         lw = 2.5,
#         label="Gauss-Seidel",
#         )
# p1 = plot!(r_3,
#         lw = 2.5,
#         label="SOR, w = 1.6",
#         )
#
# p2 = plot(pb_1,
#         xaxis = :log,
# #        title="Convergence Speed for Iterative Methods",
#         title="Periodic E-W, No Flow N-S",
#         xlabel="Iterations to Convergence",
#         ylabel="Residual",
#         lw = 2.5,
#         label="Jacobi",
#         )
# p2 = plot!(pb_2,
#         lw = 2.5,
#         label="Gauss-Seidel",
#         )
# p2 = plot!(pb_3,
#         lw = 2.5,
#         label="SOR, w = 1.6",
#         )
# pce = plot(p1,p2,layout = l)
# savefig(pce,"HW3_ItMethConv2.svg")

##F3 - Non convergence plots, last 100 iterations
# l = @layout[a; b; c]
# p1 = plot(r_3[end-100:end],
#         title="SOR, No Flow EW, 3961 Iterations",
#         lw = 2.5,
# #        label="SOR, No Flow EW",
#         legend=:none
#         )
# p2 = plot(pb_3[end-100:end],
#         lw = 2.5,
#         title="SOR, Periodic EW, 1e5 Iterations",
#         ylabel="Residual",
#         legend=:none
#         )
# p3 = plot(rmil[end-100:end],
#         lw = 2.5,
#         title="SOR, Periodic EW, 1e6 Iterations",
#         legend=:none,
#         xlabel="Last N Iterations",
#         )
# pcnum = plot(p1,p2,p3,layout = l)
# savefig(pcnum,"HW_Last_Iterations.svg")


## Non Convergence, Periodic Direction

p1 = plot(r_EW,
        label="EW-Periodic R="*string(r_EW[end]),
        xlabel="Iterations",
        ylabel="Residual",
        title="SOR Convergence Rate for Different BCs",
        xaxis=:log,
        lw=2.5
        )
p1 = plot!(r_NS,
        label="NS-Periodic R="*string(r_NS[end]),
        lw=2.5)
p1 = plot!(r_nf,
        label="All-No Flow R="*string(r_nf[end]),
        ls=:dot,
        lw=2.5)
p1 = plot!(r_ap,
        label="All-Periodic R="*string(r_ap[end]),
        ls=:dash,
        lw=2.5)
savefig(p1,"HW3_F4_DiffBCsConvergence.svg")

## F4 - Error Map, EW Periodic, 1e6 iteration


fig4= heatmap(mx,my,errmap',
        title="Residual (No Flow E/W Periodic N/S), It#"*string(itcnt),
        xlabel="x (meters)",
        ylabel="y (meters)",
        fillcolor=:inferno,
        cb=:bottom,
        clims=[-1e-13 1e-13]
        )
savefig(fig4,"HW3_F4_ErrMapSOR.svg")


#
