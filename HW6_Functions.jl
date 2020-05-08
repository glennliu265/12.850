using Plots
using Printf
using LinearAlgebra
#using PyPlot
using JLD

include("AllFunctions12850.jl")


## Choice Parameters to toggle
K = 1e2 # Horizontal Eddy Diffusivity (for Vorticity)



## Grid Set-up  -----------------------------------------------
# X and Y Grids
#xgrid = [0:1e5:5e6;]
#ygrid = [0:1e5:1e7;]
xgrid = [0:1e5:5e6;]
ygrid = [0:1e5:1e7;]

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

## Set Diffusivity Parameter
# For vorticity
κx        = Float64[K for x in mx, y in my]
κy        = Float64[K for x in mx, y in my]
κx0       = κx[1]
κy0       = κy[1]

# For streamfunction, set diffusivity to 1
κxp        = ones(Float64,1,xmax,ymax)
κyp        = ones(Float64,1,xmax,ymax)
κxp0       = κxp[1]
κyp0       = κyp[1]


## Iteration Parameters ------------------------------
tol       = 1e-7
ω         = 1.6
printint  = 1e3
method    = 3
save_iter = 1000
max_iter  = 1e5

## Time parameters
dt        = 3600*24*30     # Timestep for model integration
ts_max    = 36          # Number of timesteps to take
tottime   = ts_max*dt   # Total time elapsed in seconds
fdt       = 3600*24*30  # Time Resolution of forcing (seconds)
ftmax    = convert(Int32,ceil(tottime/fdt))
ftall     = 1:ftmax #
θ         = 0.5

## Source Term Settings -------------------------------
ρ0 = 1025 #kg/m3
H  = 1000 #m
τ0 = 0.1  #.* ones(xmax,ymax)  #N/m2


## Case Setting
casen = 2

# Create Zonal Wind Stress Terms
xedge = [0:δx:(Lx+2δx);]
yedge = [0:δy:(Ly+2δy);]

# Case 1: Constant Zonal Wind Only, increasing in time
# Seasonal fluctuation of magnitude with time
if casen ==1
    dτx = [ -τ0*cos(2*pi*y/Ly)* sin(pi/2*(m)/(ftmax)) for x in xedge, y in yedge, m in ftall]
    dτy = [ 0 for x in xedge, y in yedge, m in ftall]
end

# Case 3: Vortex Translation
if casen == 3
    dτx = [ -τ0*sin(1*pi*y/Ly-m/6) for x in xedge, y in yedge, m in ftall]
    dτy = [ -τ0*sin(1*pi*x/Lx+m/6) for x in xedge, y in yedge, m in ftall]
end

# Case 2: Meridional Component (v sinusoidal with x) which weakens with time
if casen == 2
    dτx = [ -τ0*cos(2*pi*y/Ly) for x in xedge, y in yedge, m in ftall]
    dτy = [ -τ0*cos(2*pi*x/Lx) * sin(pi*y/Ly) * sin(2*pi/2*m/ftmax) for x in xedge, y in yedge, m in ftall]
    #maskout = findall(x->(x<50 || x > 100),my)
    #dτy[:,maskout,:] .= 0
end

if casen == 4
    dτx = [-uval*sin(pi*x/Lx)*cos(2*pi*y/Ly)*sin(pi*m/ftmax) for x in xedge, y in yedge, m in ftall]
    dτy = [vval*sin(pi*x/(2*Lx))*sin(pi*m/ftmax) for x in xedge, y in yedge, m in ftall]
end
x_f2       = ones(Int8,1,length(xedge))*δx
y_f2       = ones(Int8,1,length(yedge))*δx
curlTx     = zeros(Float64,xmax,ymax,ftmax)
curlTy     = zeros(Float64,xmax,ymax,ftmax)
# Take x along y
for m = 1:ftmax
    curlTx[:,:,m],x1,y1 = ocnmod.ddx_2d(dτx[:,:,m],y_f2,xedge,yedge,2)
    curlTy[:,:,m],x1,y1 = ocnmod.ddx_2d(dτy[:,:,m],x_f2,xedge,yedge,1)
end


# -----------------------------------
# CALCULATE FORCING TERM
# -----------------------------------
    S = 1/(ρ0*H) * (curlTy - curlTx)


# -----------------------------------
#  Visualize Forcing Term, Quiver plot on wind curl
# -----------------------------------

    if casen == 2
        qscale = 5e6
    else
        qscale = 1
    end
    S0    = findmax(abs.(S))[1]
    Splot = S ./ findmax(abs.(S))[1]
    aniquiv2 = @animate for t ∈ 1:ftmax
        u = dτx[:,:,t]
        v = dτy[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(mx,my,u,v,qscale)

        a = heatmap(mx,my,Splot[:,:,t]',
            clims=(-1,1),
            title="Forcing Term Case "*string(casen)*"; t= "*string(t)*"; S0="*@sprintf("%.2e",S0),
            seriescolor=:balance,
            xlabel="x (meters)",
            ylabel="y (meters)",
            )
        a = Plots.quiver!(pts,quiver=(uv),
            ylims=(0,Ly),
            xlims=(0,Lx),
            lc="black"
            )
    end
    gif(aniquiv2,"HW4_Forcing_Case"*string(casen)*".gif",fps=5)

## Boundary Conditions
# 1 = "Dirichlet", 2 = "Neumann", 3="Periodic"

WBC = 1
wb_val = [y/y*0 for y in my]#[y/y for y in mx]

# Eastern
EBC = 1
eb_val = [y/y*0 for y in my]#[y/y for y in my]

# North
NBC = 1
nb_val = [x/x*0 for x in mx] #[sin(3*(x/Lx*pi)) for x in mx]

# South
SBC = 1
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
if NBC == 1
    nper = 1
end
if SBC == 1
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

allstart=time()
for t = 1:ts_max

        if t > 1
                x_init = vort_t[:,:,t-1]
        else
                x_init = zeros(Float64,xmax,ymax)
        end
        m = t
        if t == ts_max
            nm = 1
        else
            nm = 1
        end

        # Get Corresponding forcing term
        S0 = S[:,:,m]
        S1 = S[:,:,nm]

        # Compute Diffusivity Coefficients
        Cx0,Bx0 = ocnmod.CD_diffu_calc_coeff(x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
        Cy0,By0 = ocnmod.CD_diffu_calc_coeff(y_f',y_c',κy',NBC,nb_val,SBC,sb_val,y_c0,κy0)

        Cx1,Bx1 = ocnmod.CD_diffu_calc_coeff(x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
        Cy1,By1 = ocnmod.CD_diffu_calc_coeff(y_f',y_c',κy',NBC,nb_val,SBC,sb_val,y_c0,κy0)

        # Permute dimensions
        Cy0p = permutedims(Cy0,[1,3,2])
        Cy1p = permutedims(Cy1,[1,3,2])
        By0p = permutedims(By0,[2,1])
        By1p = permutedims(By1,[2,1])


        # Modify Forcing Term
        S0 .-= (Bx0 + By0p)
        S1 .-= (Bx1 + By1p)

        # Set up Crank Nicolson
        ug = zeros(Float64,xmax,ymax)
        Ax,Ay,b = ocnmod.CN_make_matrix_2d(dt,θ,ug,Cx0,Cy0p,Cx1,Cy1p,S0,S1,1e5,chk_per,2)

        #u_out,itcnt,rSOR,u_scrap,err_scrap,errmap = ocnmod.FD_itrsolve_2D(Ax,Ay,b,ug,tol,ω,method,wper,eper,sper,nper,max_iter,save_iter)
        #vort_t[:,:,t],itcnt2,rCGK = ocnmod.cgk_2d(Ax,Ay,b,ug,chk_per,tol,max_iter)
        vort_t[:,:,t],itcnt,r = ocnmod.cgk_2d(Ax,Ay,b,ug,chk_per,tol,max_iter,2)

        # Solve for streamfunction, using vorticity as the input (functionized)
        S_psi = vort_t[:,:,t]
        ψ_t[:,:,t],~,~ =ocnmod.invPV_2d(S_psi,
                          x_f,y_f,x_c,y_c,x_c0,y_c0,
                          NBC,nb_val,SBC,sb_val,EBC,eb_val,WBC,wb_val,
                          chk_per,
                          ug,tol,max_iter)
        # Solve for uv
        ut[:,:,t],vt[:,:,t],~,~=ocnmod.psi2uv(ψ_t[:,:,t],x_f,y_f,mx,my,1)


end
elapsed = time() - allstart
@printf("\nCompleted calculations in %s",elapsed)

# ---------------------------------------------------
# Animate Vorticity (Left) and Streamfunction (Right)
# ---------------------------------------------------
# Prepare for animating vorticity and streamfunction

if casen == 1
    qscale = 1e3
elseif casen == 2
    qscale = 5e3
elseif casen == 3
    qscale = 1e3
end
x_u = mx[2:end-1]
y_v = my[2:end-1]

# Determine maximum values
ζ0 = findmax(abs.(vort_t))[1]
plotvort = vort_t/ζ0

# Determine maximum values
ψ0     = findmax(abs.(ψ_t))[1]
plotsf = ψ_t/ψ0

anim3 = @animate for t ∈ 1:ts_max
        l = @layout[a b]

        # Prepare wind stress quicvers
        wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)

        # Make quivers
        u = ut[:,:,t]
        v = vt[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(x_u,y_v,u,v,qscale)

        #c = Plots.contourf(mx,my,vort_t[:,:,t]'/findmax(abs.(vort_t))[1],
        c = Plots.contourf(mx,my,plotvort[:,:,t]',
                clabels=true,
                levels=[-1:0.1:1;],
                clims=(-1,1),
                title="Vorticity; t="*string(t)*"; z0="*@sprintf("%.2e",ζ0),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
        c=Plots.quiver!(wspts,quiver=(wsuv),
            ylims=(0,Ly),
            xlims=(0,Lx),
            lc="black",
            )

        #h = Plots.contourf(mx,my,ψ_t[:,:,t]'/findmax(abs.(ψ_t))[1],
        h = Plots.contourf(mx,my,plotsf[:,:,t]',
                clabels=true,
                levels=[-1:0.1:1;],
                clims=(-1,1),
                title="Psi; t="*string(t)*"; psi0="*@sprintf("%.2e",ψ0),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
        h = Plots.quiver!(pts,quiver=(uv),
            ylims=(0,Ly),
            xlims=(0,Lx),
            lc="black"
            )

        Plots.plot(c,h,layout = l)
end
gif(anim3,"HW4_Vort_Psi_Case"*string(casen)*"e0_2dtest.gif",fps=5)


## save ata
#JLD.save("/Users/gyl/Case2_10m.jld","sf10",ψ_t,"u10",ut,"v10",vt,"vort10",vort_t)






### Figure Library
# # ---------------------------
# # Fig 1. Timescale vs Diffusivity Plots
# # ---------------------------
    #
    # kp = [10.0^x for x in (-2):1:5]
    # tp8 = 1e8 ./ kp
    # tp14 = 1e14 ./ kp
    # #tp4 = 1e4 ./ kp
#     tkplot = Plots.plot(kp,tp8,
#                 title = "Timescale vs Diffusivity; Domain = 1E8 m2",
#                 label = "Time Scale",
#                 xaxis = :log,
#                 yaxis = :log,
#                 xtick = (kp),
#                 ylabel= "Timescale (s)",
#                 xlabel="Diffusivity m2/s",
#                 lw    = 2.5,
#                 lc    = :black
#                 )
#     tkplot = hline!([3.145e8],label="Decadal",ls=:dot,lw=2.5,lc=:purple)
#     tkplot = hline!([3.145e7],label="Annual",ls=:dot,lw=2.5,lc=:red)
#     tkplot = hline!([2.592e6],label="Monthly",ls=:dot,lw=2.5,lc=:blue)
#     tkplot = hline!([8.840e4],label="Daily",lw=2.5,ls=:dot,lc=:green)
#     Plots.savefig(tkplot,"HW4_TimeScalevsDiffusivity.png")

# # ---------------------------------------------------------------------
# #  F2. Plot Wind Stress Curl Case 1
# # ---------------------------------------------------------------------
#     l = @layout[a b]
#     # Plot Wind Stress
#     qscale = 5e3
#
#     t = 1
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(abs.(curlTy-curlTx))[1]
#     #curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     p1=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     p1=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     t = 36
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     p2=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     p2=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     # Plot Wind Stress Curl
#     pws = Plots.plot(p1,p2,layout=l)
#
#     Plots.savefig(pws,"HW4_WindStressPlotCase"*string(casen)*".svg")


# ---------------------------------------------------------------------
#  F3. Plot Wind Stress Curl (Case 2)
# ---------------------------------------------------------------------
#     l = @layout[a b]
#     # Plot Wind Stress
#     qscale = 5e3
#
#     t = 1
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(abs.(curlTy-curlTx))[1]
#     #curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     p1=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     p1=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     t = 18
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     p2=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     p2=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     # Plot Wind Stress Curl
#     pws = Plots.plot(p1,p2,layout=l)
#
#     Plots.savefig(pws,"HW4_WindStressPlot.svg")


# # ---------------------------------------------------------------------
# #  F4. Plot Wind Stress Curl Case 3
# # ---------------------------------------------------------------------
#     l = @layout[a b ; c d]
#     # Plot Wind Stress
#     qscale = 5e3
#
#     t = 1
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(abs.(curlTy-curlTx))[1]
#     #curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     pa=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     pa=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     t = 9
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(abs.(curlTy-curlTx))[1]
#     #curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     pb=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     pb=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     t = 28
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(abs.(curlTy-curlTx))[1]
#     #curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     pc=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     pc=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     t = 36
#     plotcurl = curlTy[:,:,t]-curlTx[:,:,t]
#     curl0 = findmax(plotcurl)[1]
#     plotcurl = plotcurl / curl0
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
#     pd=Plots.heatmap(mx,my,plotcurl',seriescolor=:balance)
#     pd=Plots.quiver!(wspts,quiver=(wsuv),
#         title="Tau @ t="*string(t)*"; curl0="*@sprintf("%.2e",curl0),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         clims=(-1,1),
#         lc="black",
#         ylabel="y (meters)",
#         xlabel="x (meters)"
#         )
#
#     # Plot Wind Stress Curl
#     pws = Plots.plot(pa,pb,pc,pd,layout=l)
#
#     Plots.savefig(pws,"HW4_WindStressPlotCase"*string(casen)*".svg")


# # ---------------------------
# # Fig 3. Comparison of Iterative Methods
# # ---------------------------
    #
    # p1 = Plots.plot(rSOR,
    #         xaxis = :log,
    #         yaxis = :log,
    # #        title="Convergence Speed for Iterative Methods",
    #         title="Comparison of Iterative Methods",
    #         xlabel="Iterations to Convergence",
    #         ylabel="Residual",
    #         lw = 2.5,
    #         label="SOR",
    #         )
    # p1 = Plots.plot!(rCGK,
    #         lw = 2.5,
    #         label="CGK",
    #         )
    # Plots.savefig(p1,"HW4_Convergence_SpeedTest2_10m.svg")

#
# #
# # Figure: Case 2 Difference
# #
#
#
#     qscale = 5e3
#     # Determine maximum values
#     ζ0 = findmax(abs.(vort_t))[1]
#     plotvort = vort_t/ζ0
#
#     # Determine maximum values
#     ψ0     = findmax(abs.(ψ_t))[1]
#     plotsf = ψ_t/ψ0
#
#     # Take Differences for plotting
#     pvdiff = plotvort[:,:,18] - plotvort[:,:,1]
#     sfdiff = plotsf[:,:,18]   - plotsf[:,:,1]
#     dtx = dτx[:,:,18] - dτx[:,:,1]
#     dty = dτy[:,:,18] - dτy[:,:,1]
#     u = ut[:,:,18] - ut[:,:,1]
#     v = vt[:,:,18] - vt[:,:,1]
#
#     # Prepare wind stress quicvers
#     wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dtx,dty,qscale)
#     pts,uv = ocnmod.quiverprep_2d(x_u,y_v,u,v,qscale)
#
#      l = @layout[a b]
#     c = Plots.contourf(mx,my,pvdiff',
#             clabels=true,
#             levels=[-.5:0.1:.5;],
#             clims=(-0.5,0.5),
#             title="Vorticity Difference (t = 1 to 18); z0="*@sprintf("%.2e",ζ0),
#             xlabel="x (meters)",
#             ylabel="y (meters)",
#             fillcolor=:balance
#             )
#     c=Plots.quiver!(wspts,quiver=(wsuv),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         lc="black",
#         )
#
#     h = Plots.contourf(mx,my,sfdiff',
#             clabels=true,
#             levels=[-.5:0.1:.5;],
#             clims=(-0.5,0.5),
#             title="Psi difference (t = 1 to 18); psi0="*@sprintf("%.2e",ψ0),
#             xlabel="x (meters)",
#             ylabel="y (meters)",
#             fillcolor=:balance
#             )
#     h = Plots.quiver!(pts,quiver=(uv),
#         ylims=(0,Ly),
#         xlims=(0,Lx),
#         lc="black"
#         )
#
#     dp2 =Plots.plot(c,h,layout = l)
# Plots.savefig(dp2,"HW4_Vort_Psi_Difference_Case"*string(casen)*".svg")





#
# Figure: Case 2 Testing K values
#
    # Run once with k = 10^2 (vortt2, psit2, etc)
    # Run again with k = 10^0 (vort_t,psi_t,etc)

    qscale = 5e3

    vortd = vortt_2 - vort_t;
    psid  = psit_2 - ψ_t;
    utd   = ut2 - ut;
    vtd   = vt2 - vt;


    # Determine maximum values
    ζ0 = findmax(abs.(vortd))[1]
    plotvort = vortd/ζ0

    # Determine maximum values
    ψ0     = findmax(abs.(psid))[1]
    plotsf = psid/ψ0

    t = 18
    # Prepare wind stress quicvers
    wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dtx,dty,qscale)
    pts,uv = ocnmod.quiverprep_2d(x_u,y_v,utd[:,:,t],vtd[:,:,t],qscale)

     l = @layout[a b]
    c = Plots.contourf(mx,my,plotvort[:,:,t]',
            clabels=true,
            levels=[-1:0.1:1;],
            clims=(-1,1),
            title="Vorticity(2 - 0); z0="*@sprintf("%.2e",ζ0),
            xlabel="x (meters)",
            ylabel="y (meters)",
            fillcolor=:balance
            )
    c=Plots.quiver!(wspts,quiver=(wsuv),
        ylims=(0,Ly),
        xlims=(0,Lx),
        lc="black",
        )

    h = Plots.contourf(mx,my,plotsf[:,:,t]',
            clabels=true,
            levels=[-1:0.1:1;],
            clims=(-1,1),
            title="Psi (2 - 0); psi0="*@sprintf("%.2e",ψ0),
            xlabel="x (meters)",
            ylabel="y (meters)",
            fillcolor=:balance
            )
    h = Plots.quiver!(pts,quiver=(uv),
        ylims=(0,Ly),
        xlims=(0,Lx),
        lc="black"
        )

    dp2 =Plots.plot(c,h,layout = l)
Plots.savefig(dp2,"HW4_Vort_Psi_Difference_Case"*string(casen)*".svg")
