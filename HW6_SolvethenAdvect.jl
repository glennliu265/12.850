using Plots
using Printf
using LinearAlgebra
#using PyPlot
using JLD

include("AllFunctions12850.jl")

## Choice Parameters to toggle
K = 1e2 # Horizontal Eddy Diffusivity (for Vorticity)
saveiter = 1;


## Grid Set-up  -----------------------------------------------
# X and Y Grids
#xgrid = [0:1e5:5e6;]
#ygrid = [0:1e5:1e7;]
xgrid = [0:2.5e5:5e6;]
ygrid = [0:2.5e5:1e7;]

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
dt        = 3600*24     # Timestep for model integration
ts_max    = 365        # Number of timesteps to take
θ         = 0.5
tmaxval   = dt*ts_max   # Maximum timestep value in seconds

## Source Term Settings -------------------------------
ρ0 = 1025 #kg/m3
H  = 1000 #m
τ0 = 0.1  #.* ones(xmax,ymax)  #N/m2


## Case Setting
casen = 1

# Create Zonal Wind Stress Terms
xedge = [0:δx:(Lx+2δx);]
yedge = [0:δy:(Ly+2δy);]

# Case 1: Constant Zonal Wind Only, increasing in time
# Seasonal fluctuation of magnitude with time
if casen ==1
    dτx = [ -τ0*cos(2*pi*y/Ly) for x in xedge, y in yedge, t in 1:ts_max]
    dτy = [ 0 for x in xedge, y in yedge, t in 1:ts_max]
end

# # Case 2: Meridional Component (v sinusoidal with x) which weakens with time
# if casen == 2
#     dτx = [ -τ0*cos(2*pi*y/Ly) for x in xedge, y in yedge, m in ftall]
#     dτy = [ -τ0*cos(2*pi*x/Lx) * sin(pi*y/Ly) * sin(2*pi/2*m/ftmax) for x in xedge, y in yedge, m in ftall]
#     #maskout = findall(x->(x<50 || x > 100),my)
#     #dτy[:,maskout,:] .= 0
# end

# if casen == 4
#     dτx = [-uval*sin(pi*x/Lx)*cos(2*pi*y/Ly)*sin(pi*m/ftmax) for x in xedge, y in yedge, m in ftall]
#     dτy = [vval*sin(pi*x/(2*Lx))*sin(pi*m/ftmax) for x in xedge, y in yedge, m in ftall]
# end
x_f2       = ones(Int8,1,length(xedge))*δx
y_f2       = ones(Int8,1,length(yedge))*δx
curlTx     = zeros(Float64,xmax,ymax,ts_max)
curlTy     = zeros(Float64,xmax,ymax,ts_max)
# Take x along y
for m = 1:ts_max
    curlTx[:,:,m],x1,y1 = ocnmod.ddx_2d(dτx[:,:,m],y_f2,xedge,yedge,2)
    curlTy[:,:,m],x1,y1 = ocnmod.ddx_2d(dτy[:,:,m],x_f2,xedge,yedge,1)
end


# -----------------------------------
# CALCULATE FORCING TERM
# -----------------------------------
    S = 1/(ρ0*H) * (curlTy - curlTx)

## Boundary Conditions
# 1 = "Dirichlet", 2 = "Neumann", 3="Periodic"

WBC = 1
wb_val = [y/y*0 for y in my, t in 1:ts_max]#[y/y for y in mx]

# Eastern
EBC = 1
eb_val = [y/y*0 for y in my, t in 1:ts_max]#[y/y for y in my]

# North
NBC = 1
nb_val = [x/x*0 for x in mx, t in 1:ts_max] #[sin(3*(x/Lx*pi)) for x in mx]

# South
SBC = 1
sb_val = [x/x*0 for x in mx, t in 1:ts_max]#[sin(3*(x/Lx*pi)) for x in mx]


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

tsdim = Int16[]
push!(tsdim,ceil(ts_max/saveiter))


# Save variables
zall = zeros(Float64,xmax,ymax,tsdim[1])
ψ_t = zeros(Float64,xmax,ymax,tsdim[1])
uall = zeros(Float64,xmax,ymax,tsdim[1])
vall = zeros(Float64,xmax,ymax,tsdim[1])

# Dummy Variables
x_guess = zeros(Float64,xmax,ymax)
dummyC = zeros(3,xmax,ymax)
dummyS = zeros(xmax,ymax)
allstart=time()

# Previous Iteration
zprev    = zeros(Float64,xmax,ymax) # Vorticity IC
uprev    = zeros(Float64,xmax,ymax)
vprev    = zeros(Float64,xmax,ymax)

for t = 1:ts_max
        global uprev, vprev, zprev
        u0 = uprev
        v0 = vprev
        z0 = zprev
        # # Take vorticity and velocities from last timestep
        # if t > 1
        #         x_init = z0
        #         u0     = u0
        #         v0     = v0
        # # On first timestep, begin from state of rest
        # else
        #         x_init = zeros(Float64,xmax,ymax)
        #         u0 = zeros(Float64,xmax,ymax)
        #         v0 = zeros(Float64,xmax,ymax)
        # end

        m = t
        if t == ts_max
            nm = 1
        else
            nm = t+1
        end

        ## ------------------------
        # Solve for Vorticity (HW4)
        ## ------------------------

        # Compute Diffusivity Coefficients (x and y (and permute y))
        Cx0,Bx0 = ocnmod.CD_diffu_calc_coeff(x_f ,x_c ,κx ,EBC,eb_val[:,t],WBC,wb_val[:,t],x_c0,κx0)
        Cy0,By0 = ocnmod.CD_diffu_calc_coeff(y_f',y_c',κy',NBC,nb_val[:,t],SBC,sb_val[:,t],y_c0,κy0)
        Cx1,Bx1 = ocnmod.CD_diffu_calc_coeff(x_f ,x_c ,κx ,EBC,eb_val[:,nm],WBC,wb_val[:,nm],x_c0,κx0)
        Cy1,By1 = ocnmod.CD_diffu_calc_coeff(y_f',y_c',κy',NBC,nb_val[:,nm],SBC,sb_val[:,nm],y_c0,κy0)


        # Permute dimensions (p for permute)
        Cy0p   = permutedims(Cy0,[1,3,2])
        By0p   = permutedims(By0,[2,1])
        Cy1p   = permutedims(Cy1,[1,3,2])
        By1p   = permutedims(By1,[2,1])

        # Modify Forcing Term
        S0 = S[:,:,t]
        S1 = S[:,:,nm]
        S0 .-= (Bx0 + By0p)
        S1 .-= (Bx0 + By0p)

        # Set up Crank Nicolson
        Ax,Ay,b = ocnmod.CN_make_matrix_2d(dt,θ,z0,Cx0,Cy0p,Cx1,Cy1p,S0,S1,1e5,chk_per,2)

        # Solve for Vorticity (use last timestep as guess)
        z0,itcnt,r = ocnmod.cgk_2d(Ax,Ay,b,z0,chk_per,tol,max_iter,2)

        # Invert for streamfunction
        ψ_t,~,~ =ocnmod.invPV_2d(z0,
                          x_f,y_f,x_c,y_c,x_c0,y_c0,
                          NBC,nb_val,SBC,sb_val,EBC,eb_val,WBC,wb_val,
                          chk_per,
                          z0,tol,max_iter,1)
        # Solve for UV
        u0,v0,~,~=ocnmod.psi2uv(ψ_t,x_f,y_f,mx,my,
                                    0,nb_val[:,t],sb_val[:,t],eb_val[:,t],wb_val[:,t])

        # Save velocities at specified iterations


        # --------------------
        # Advect PV
        # --------------------
        Ad_x0,Ad_Bx0 = ocnmod.UW_calc_coeff_2D(x_f,xmax,u0 ,EBC,eb_val[:,t],WBC,wb_val[:,t])
        Ad_y0,Ad_By0 = ocnmod.UW_calc_coeff_2D(y_f,ymax,v0',NBC,nb_val[:,t],SBC,sb_val[:,t])


        # Permute dimensions (p for permute)
        Ad_y0p   = permutedims(Ad_y0,[1,3,2])
        Ad_By0p  = permutedims(Ad_By0,[2,1])

        # Forcing Term is not there, so it will just be Ad_By0
        Ad_S0    = Ad_By0p .*-1


        # FEM advection
        ~,~,z0 = ocnmod.CN_make_matrix_2d(dt,0,z0,Ad_x0,Ad_y0p,dummyC,dummyC,Ad_S0,dummyS,1e5,chk_per,2)

        if t%saveiter == 0
            idx = Int16[]
            push!(idx,ceil(t/saveiter))
            uall[:,:,idx[1]] = u0
            vall[:,:,idx[1]] = v0
            zall[:,:,idx[1]] = z0
            elapsed = time() - allstart
            @printf("\nCompleted iteration %i in %s",t,elapsed)
        end

        zprev = z0
        vprev = v0
        uprev = u0

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
qscale = 5e6
x_u = mx[2:end-1]
y_v = my[2:end-1]

# Determine maximum values
ζ0 = findmax(abs.(zall))[1]
plotvort = vort_t/ζ0

# Determine maximum values
ψ0     = findmax(abs.(zall))[1]
plotsf = ψ_t/ψ0

anim3 = @animate for t ∈ 1:10
        l = @layout[a b]
        #
        # wspts,wsuv = ocnmod.quiverprep_2d(xedge,yedge,dτx[:,:,t],dτy[:,:,t],qscale)
        #
        # Make quivers
        u = uall[:,:,t]
        v = vall[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(mx,my,u,v,qscale)

        c = Plots.contourf(mx,my,zall[:,:,t]',
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
        # h = Plots.contourf(mx,my,plotsf[:,:,t]',
        #         clabels=true,
        #         levels=[-1:0.1:1;],
        #         clims=(-1,1),
        #         title="Psi; t="*string(t)*"; psi0="*@sprintf("%.2e",ψ0),
        #         xlabel="x (meters)",
        #         ylabel="y (meters)",
        #         fillcolor=:balance
        #         )
        h = Plots.quiver(pts,quiver=(uv),
            ylims=(0,Ly),
            xlims=(0,Lx),
            lc="black"
            )

        Plots.plot(c,h,layout = l)
end
gif(anim3,"HW6_Case"*string(casen)*"e0_2dtest.gif",fps=5)




var = zall
anim3 = @animate for t ∈ 1:5:365

        c = Plots.contourf(mx,my,var[:,:,t]',
                clabels=true,
                title="Vorticity; t="*string(t),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
end
gif(anim3,"HW6_Case"*string(casen)*"e0_2dtest.gif",fps=20)




qscale = 5e5
anim3 = @animate for t ∈ 1:5:365
        u = uall[:,:,t]
        v = vall[:,:,t]
        pts,uv = ocnmod.quiverprep_2d(mx,my,u,v,qscale)
        c = Plots.contourf(mx,my,zall[:,:,t]',
                clabels=true,
                title="Vorticity; t="*string(t),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
        h = Plots.quiver!(pts,quiver=(uv),
            title="Vorticity; t="*string(t),
            ylims=(0,Ly),
            xlims=(0,Lx),
            lc="black"
            )
end
gif(anim3,"HW6_Case"*string(casen)*"e0_quivtest.gif",fps=20)
