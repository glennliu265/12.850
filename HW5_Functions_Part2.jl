using Plots
using LinearAlgebra
using Printf
include("AllFunctions12850.jl")


## Grid Set-up  -----------------------------------------------
# X and Y Grids
#xgrid = [0:1e5:5e6;]
#ygrid = [0:1e5:1e7;]
xgrid = [0:100:5000;]
ygrid = [0:100:10000;]

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

## Velocity Field and IC -------------------------------


u = [-0.1*cos(pi*y/Ly) for x in xgrid, y in ygrid, t in t:ts_max]
v = [0 for x in xgrid, y in ygrid]

T0 = [5*sin(pi*y/Ly) + 5*sin(pi*x/Lx) for x in mx, y in my, t in t:ts_max]
#T0[T0 .< 2] .= 0

wcpts,wcuv = ocnmod.quiverprep_2d(xgrid,ygrid,u[:,:,t],v[:,:,t],5e3)
pwc = heatmap(mx,my,T0',
        seriescolor=:dense,
        clims=(0,12))
pwc=Plots.quiver!(wcpts,quiver=(wcuv),
        linecolor=:black,
        title="Velocity Field",
        ylims=(0,1e4),
        xlims=(0,5e3),
        xlabel="Zonal (m)",
        ylabel="Meridional (m)")


## Forcing Term -------------------------------

S = [0 for x in mx, y in my, t in 1:ts_max]

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
temps = zeros(Float64,xmax,ymax,ts_max)


allstart=time()
for t = 1:ts_max

        if t == ts_max
                tp1 = 1
        else
                tp1 = t+1
        end

        # Get Corresponding forcing term
        S0 = S[:,:,t]
        S1 = S[:,:,tp1]

        # Compute Coefficients

        # x- Coefficients
        Cx0,Bx0   = ocnmod.UW_calc_coeff_2D(x_f,xmax,u[:,:,t],EBC,eb_val,WBC,wb_val)
        Cx1,Bx1   = ocnmod.UW_calc_coeff_2D(x_f,xmax,u[:,:,tp1],EBC,eb_val,WBC,wb_val)


        Cy0,By0   = ocnmod.UW_calc_coeff_2D(y_f,ymax,u[:,:,t]',NBC,nb_val,SBC,sb_val)
        Cy1,By1   = ocnmod.UW_calc_coeff_2D(y_f,ymax,u[:,:,tp1]',NBC,nb_val,SBC,sb_val)

        # Modify Boundary Conditions
        S0 .+= (Bx0 .+ By0')
        S1 .+= (Bx1 .+ By1')


        # Set up Crank Nicolson
        Ax,Ay,b = ocnmod.CN_make_matrix_2d(dt,θ,T0,Cx0,Cy0,Cx1,Cy1,S0,S1,1e5,chk_per)

        #u_out,itcnt,rSOR,u_scrap,err_scrap,errmap = ocnmod.FD_itrsolve_2D(Ax,Ay,b,ug,tol,ω,method,wper,eper,sper,nper,max_iter,save_iter)
        #vort_t[:,:,t],itcnt2,rCGK = ocnmod.cgk_2d(Ax,Ay,b,ug,chk_per,tol,max_iter)
        vort_t[:,:,t],~,~ = ocnmod.cgk_2d(Ax,Ay,b,ug,chk_per,tol,max_iter)

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
elapsed = time() - allstart
@printf("\nCompleted calculations in %s",elapsed)
