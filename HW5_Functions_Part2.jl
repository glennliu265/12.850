using Plots
using LinearAlgebra
using Printf
include("AllFunctions12850.jl")


## Grid Set-up  -----------------------------------------------
# X and Y Grids
#xgrid = [0:1e5:5e6;]
#ygrid = [0:1e5:1e7;]
xgrid = [0:1:50;]
ygrid = [0:1:100;]

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
max_iter  = 1e3

## Time parameters
dt        = 1#3600*24*30     # Timestep for model integration
ts_max    = 20         # Number of timesteps to take
θ         = 0

## Velocity Field and IC -------------------------------


#u = [1e3*cos(pi*y/Ly) for x in xgrid, y in ygrid, t in t:ts_max]
u = [-1 for x in xgrid, y in ygrid, t in 1:ts_max]

v = [1 for x in xgrid, y in ygrid, t in 1:ts_max]

T0 = [5*sin(pi*y/Ly) + 5*sin(pi*x/Lx) for x in mx, y in my]
#T0[T0 .< 2] .= 0

wcpts,wcuv = ocnmod.quiverprep_2d(xgrid,ygrid,u[:,:,t],v[:,:,t],2)
pwc = heatmap(mx,my,T0',
        seriescolor=:dense,
        clims=(0,12))
pwc=Plots.quiver!(wcpts,quiver=(wcuv),
        linecolor=:black,
        title="Velocity Field",
        ylims=(0,ymax),
        xlims=(0,xmax),
        xlabel="Zonal (m)",
        ylabel="Meridional (m)")


## Forcing Term -------------------------------

S = [0 for x in mx, y in my, t in 1:ts_max]

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
temps = zeros(Float64,xmax,ymax,ts_max)


allstart=time()
for t = 1:ts_max

        if t == ts_max
                tp1 = 1
        else
                tp1 = t+1
        end

        if t == 1
                IC = T0
        else
                IC = temps[:,:,t-1]
        end

        # Get Corresponding forcing term
        S0 = S[:,:,t]
        S1 = S[:,:,tp1]

        # Compute Coefficients
        # x- Coefficients
        Cx0,Bx0   = ocnmod.UW_calc_coeff_2D(x_f,xmax,u[:,:,t],EBC,eb_val,WBC,wb_val)
        Cx1,Bx1   = ocnmod.UW_calc_coeff_2D(x_f,xmax,u[:,:,tp1],EBC,eb_val,WBC,wb_val)

        Cy0,By0   = ocnmod.UW_calc_coeff_2D(y_f,ymax,v[:,:,t]',NBC,nb_val,SBC,sb_val)
        Cy1,By1   = ocnmod.UW_calc_coeff_2D(y_f,ymax,v[:,:,tp1]',NBC,nb_val,SBC,sb_val)

        # Permute dimensions
        Cy0 = permutedims(Cy0,[1,3,2])
        Cy1 = permutedims(Cy1,[1,3,2])
        By0 = By0'
        By1 = By1'

        # Modify Boundary Conditions
        S0 .+= (Bx0 .+ By0)
        S1 .+= (Bx1 .+ By1)


        # Set up Crank Nicolson
        Ax,Ay,b = ocnmod.CN_make_matrix_2d(dt,θ,IC,Cx0,Cy0,Cx1,Cy1,S0,S1,1e5,chk_per,2)

        temps[:,:,t],itcnt,r = ocnmod.cgk_2d(Ax,Ay,b,IC,chk_per,tol,max_iter,2)
        #temps[:,:,t],itcnt,rSOR,~,~,~ = ocnmod.FD_itrsolve_2D(Ax,Ay,b,IC,tol,ω,3,wper,eper,sper,nper,max_iter,save_iter,2)

end
elapsed = time() - allstart
@printf("\nCompleted calculations in %s",elapsed)



anim3 = @animate for t ∈ 1:ts_max
        c = heatmap(mx,my,temps[:,:,t]',
                clabels=true,
                #levels=[-1:0.1:1;],
                clims=(0,10),
                title="Temp; t="*string(t),
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:balance
                )
end
gif(anim3,"HW5_TempAdv_test.gif",fps=5)
