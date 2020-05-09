using Plots
using LinearAlgebra
using Printf
include("AllFunctions12850.jl")

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


## Iteration Parameters ------------------------------
tol       = 1e-7
ω         = 1.6
printint  = 1e3
method    = 3
save_iter = 1000
max_iter  = 1e3

## Time parameters
dt        = 3600*24*14     # Timestep for model integration
ts_max    = 100        # Number of timesteps to take
θ         = 0


## Courant to evaluate stability
uval = 0.04
vval = 0.04
courant = uval*dt/δx +  vval*dt/δy


## Velocity Field and IC -------------------------------
case     = 3


#u = [1e3*cos(pi*y/Ly) for x in xgrid, y in ygrid, t in t:ts_max]
#u = [uval*sin(2*pi*t/ts_max)*cos(2*pi*y/Ly) for x in xgrid, y in ygrid, t in 1:ts_max]
u = [-uval*cos(2*pi*y/Ly)*sin(2*pi*t/ts_max)  for x in xgrid, y in ygrid, t in 1:ts_max]
v = [vval*cos(2*pi*x/Lx)*sin(2*pi*t/ts_max) for x in xgrid, y in ygrid, t in 1:ts_max]
#v = [vval*sin(2*pi*t/ts_max)*cos(2*pi*x/Lx)*sin(2*pi*t/ts_max) for x in xgrid, y in ygrid, t in 1:ts_max]

#u =  [-uval*sin(pi*x/Lx)*cos(2*pi*y/Ly)*sin(2*pi*t/ts_max) for x in xgrid, y in ygrid, t in 1:ts_max]
#v =  [vval*sin(2*pi*x/(Lx))*cos(pi*y/Ly)*sin(2*pi*t/ts_max) for x in xgrid, y in ygrid, t in 1:ts_max]

# Sinusoidal
#if case == 2
        #T0 = [0 for x in mx, y in my]
        #T0 = [10*sin(pi*y/Ly) + 5*sin(pi*x/Lx) for x in mx, y in my]
# Static Box
if case == 1 || case == 2
        xmask = xmax/2-5 .<  mx .< xmax/2+5
        ymask = ymax/2-5 .<  my .< ymax/2+5
        mask = [x+y for x in xmask, y in ymask]
        mask = mask .>= 2
        T0 = ones(Float64,xmax,ymax)*10 .* mask
        #T0[T0 .< 2] .= 0
elseif case == 4
        T0 = [0 for x in mx, y in my]
elseif case == 3
        T0 = [10-10*(y/Ly) for x in mx, y in my]
end


t = ts_max
wcpts,wcuv = ocnmod.quiverprep_2d(xgrid,ygrid,u[:,:,t],v[:,:,t],0.5)
pwc = heatmap(mx,my,T0',
        seriescolor=:dense,
        clims=(0,12))
pwc=Plots.quiver!(wcpts,quiver=(wcuv),
        linecolor=:black,
        title="Velocity Field",
        ylims=(0,Ly),
        xlims=(0,Lx),
        xlabel="Zonal (m)",
        ylabel="Meridional (m)")


## Forcing Term -------------------------------
# First case: no forcing.
S = Float64[0 for x in mx, y in my, t in 1:ts_max]

# Second, sinusoidal point
xmask = xmax/2-5 .<  mx .< xmax/2+5
ymask = ymax/2-5 .<  my .< ymax/2+5
mask = [x+y for x in xmask, y in ymask]
mask = mask .>= 2
#S = Float64[5 for x in mx, y in my, t in 1:ts_max]
#S = Float64[1 + sin(pi*y/Ly)*sin(pi*x/Lx)*sin(2*pi*t/ts_max) for x in mx, y in my, t in 1:ts_max]
#S .*= mask #[S[:,:,t] .* mask for t in 1:ts_max]
## Boundary Conditions
# 1 = "Dirichlet", 2 = "Neumann", 3="Periodic"
# West
# 1 = Dirichlet, 2 = Neumann, 3 = Periodic

WBC = 1
if case == 3
        wb_val = [10-y/Ly*10 for y in my, t in 1:ts_max]
elseif case == 2
        wb_val = [10*sin(2*pi*t/ts_max) for y in my, t in 1:ts_max]
else
        wb_val = [y/y*0 for y in my, t in 1:ts_max]#[y/y for y in mx]
#wb_val = [5*sin(2*pi*t/ts_max) for y in my, t in 1:ts_max]
end
# Eastern
EBC = 1
if case == 3
        eb_val = [10-y/Ly*10 for y in my, t in 1:ts_max]
else
        eb_val = [y/y*0 for y in my, t in 1:ts_max]
end

# North
NBC = 1
if case == 3
        nb_val = [0 for x in mx, t in 1:ts_max]
else

        nb_val = [x/x*0 for x in mx, t in 1:ts_max] #[sin(3*(x/Lx*pi)) for x in mx]
        #nb_val = [5*sin(2*pi*t/ts_max) for x in mx, t in 1:ts_max]
end
# South
SBC = 1
if case == 3
        sb_val = [10 for x in mx, t in 1:ts_max]
else

        #sb_val = [x/x*0 for x in mx, t in 1:ts_max] #[sin(3*(x/Lx*pi)) for x in mx]
        sb_val = [x/x*0 for x in mx, t in 1:ts_max]#[sin(3*(x/Lx*pi)) for x in mx]
end


#courant =
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
        Cx0,Bx0   = ocnmod.UW_calc_coeff_2D(x_f,xmax,u[:,:,t],EBC,eb_val[:,t],WBC,wb_val[:,t])
        Cx1,Bx1   = ocnmod.UW_calc_coeff_2D(x_f,xmax,u[:,:,tp1],EBC,eb_val[:,tp1],WBC,wb_val[:,tp1])

        Cy0,By0   = ocnmod.UW_calc_coeff_2D(y_f,ymax,v[:,:,t]',NBC,nb_val[:,t],SBC,sb_val[:,t])
        Cy1,By1   = ocnmod.UW_calc_coeff_2D(y_f,ymax,v[:,:,tp1]',NBC,nb_val[:,tp1],SBC,sb_val[:,tp1])

        # Permute dimensions
        Cy0p = permutedims(Cy0,[1,3,2])
        Cy1p = permutedims(Cy1,[1,3,2])
        By0p = permutedims(By0,[2,1])
        By1p = permutedims(By1,[2,1])

        # Modify Boundary Conditions
        S0 .+= (Bx0 .+ By0p)
        S1 .+= (Bx1 .+ By1p)


        # Set up Crank Nicolson
        Ax,Ay,b = ocnmod.CN_make_matrix_2d(dt,θ,IC,Cx0,Cy0p,Cx1,Cy1p,S0,S1,1e5,chk_per,2)

        #temps[:,:,t],itcnt,r = ocnmod.cgk_2d(Ax,Ay,b,IC,chk_per,tol,max_iter,2)
        temps[:,:,t],itcnt,rSOR,~,~,~ = ocnmod.FD_itrsolve_2D(Ax,Ay,b,IC,tol,ω,1,wper,eper,sper,nper,max_iter,save_iter,2)

end
elapsed = time() - allstart
@printf("\nCompleted calculations in %s",elapsed)



anim3 = @animate for t ∈ 1:ts_max
        wcpts,wcuv = ocnmod.quiverprep_2d(xgrid,ygrid,u[:,:,t],v[:,:,t],0.5)
        c = contourf(mx,my,temps[:,:,t]',
                clabels=true,
                #levels=[-1:0.1:1;],
                clims=(0,10),
                title="Temperature (degC); t="*string(t)*"s",
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:dense,
                dpi=200
                )
        pwc=Plots.quiver!(wcpts,quiver=(wcuv),
                linecolor=:black,
                #title="Velocity Field",
                ylims=(0,Ly),
                xlims=(0,Lx))
end
gif(anim3,"HW5_Temperature Advection_case"*string(case)*".gif",fps=10)



# Static Plot

l = @layout[a b; c d]

        t = 1
        wcpts,wcuv = ocnmod.quiverprep_2d(xgrid,ygrid,u[:,:,t],v[:,:,t],0.5)
        a = contourf(mx,my,temps[:,:,t]',
                clabels=true,
                #levels=[-1:0.1:1;],
                clims=(0,10),
                title="t="*string(t)*"s",
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:dense,
                dpi=200
                )
        pwc=Plots.quiver!(wcpts,quiver=(wcuv),
                linecolor=:black,
                #title="Velocity Field",
                ylims=(0,ymax),
                xlims=(0,xmax))

        t = 1
        wcpts,wcuv = ocnmod.quiverprep_2d(xgrid,ygrid,u[:,:,t],v[:,:,t],0.5)
        a = contourf(mx,my,temps[:,:,t]',
                clabels=true,
                #levels=[-1:0.1:1;],
                clims=(0,10),
                title="t="*string(t)*"s",
                xlabel="x (meters)",
                ylabel="y (meters)",
                fillcolor=:dense,
                dpi=200
                )
        pwc=Plots.quiver!(wcpts,quiver=(wcuv),
                linecolor=:black,
                #title="Velocity Field",
                ylims=(0,Ly),
                xlims=(0,Ly))
