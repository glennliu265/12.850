using Plots
using Printf
using LinearAlgebra

include("AllFunctions12850.jl")


## Grid Set-up  -----------------------------------------------
# X and Y Grids
xgrid = [0:1:50;]
ygrid = [0:1:150;]
# xgrid = [0:.1:2;]
# ygrid = [0:.1:1;]

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
# Constant-value
κx        = ones(Float64,1,xmax)
κy        = ones(Float64,1,ymax)

κx0       = κx[1]
κy0       = κy[1]

## Iteration Parameters ------------------------------
tol       = 1e-6
ω         = 1.6
printint  = 1e4
method    = 3
save_iter = 1000
max_iter  = 1e5

## Source Term Settings
ζ = zeros(Float64,xmax,ymax)


#ζ = [sin(x/Lx*pi)^2 + sin(y/Ly*pi)^2 for x in mx, y in my] #.*0
#[sin(x/Lx*pi)^2 + cos(y/Ly*pi)^2 for x in mx, y in my]
# ζ[1,75]  = 10
# ζ[40,150] = -10
ζ[2,7] = -10
ζ[3,7] = 10
#ζ[25,75] = 2

#ones(Float64,xmax,ymax).*10#
contourf(mx,my,ζ')
 #heatmap(mx,my,ζ')
## Boundary Conditions
# 1 = "Dirichlet", 2 = "Neumann", 3="Periodic"
# West
# 1 = Dirichlet, 2 = Neumann, 3 = Periodic

WBC = 3
wb_val = [y/y*0 for y in my]#[y/y for y in mx]

# Eastern
EBC = 3
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

# Compute Coefficients and modifications to BCs
Cx,Bx = ocnmod.FD_calc_coeff_2D(xmax,x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
Cy,By = ocnmod.FD_calc_coeff_2D(ymax,y_f,y_c,κy,NBC,nb_val,SBC,sb_val,y_c0,κy0)

# Modify Boundary Conditions
ζ[1,:]    -= Bx[1,:] # South BC
ζ[xmax,:] -= Bx[2,:] # North BC
ζ[:,1]    -= By[1,:] # West BC
ζ[:,ymax] -= By[2,:] # East BC

# # # note that this assumes xmax = ymax (equal sized spacing on x and y)
# A = zeros(Float64,5,xmax)
# A[1,:] = Cy[1,:]            # [  i , j-1]
# A[2,:] = Cx[1,:]            # [i-1 ,   j]
# A[3,:] = Cx[2,:] .+ Cy[2,:] # [  i ,   j]
# A[4,:] = Cx[3,:]            # [i+1 ,   j]
# A[5,:] = Cy[3,:]            # [  i , j+1]

S=ζ
ug = ones(Float64,xmax,ymax)

u_out,itcnt,r,u_scrap,err_scrap,errmap = ocnmod.FD_itrsolve_2D(Cx,Cy,S,ug,tol,ω,method,wper,eper,sper,nper,max_iter,save_iter)

## F4 - Error Map, EW Periodic, 1e6 iteration


fig5= contourf(mx,my,u150',
#        clims = (-.1,.1),
        clabels=true,
        levels=20,
        title="Staggered Pt Vortices (No Flow E/W Periodic N/S), It#"*string(itcnt),
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
