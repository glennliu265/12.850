using Plots
using Printf
using LinearAlgebra

include("AllFunctions12850.jl")


## Grid Set-up  -----------------------------------------------
# X and Y Grids
xgrid = [0:1:100;]
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
Lx = mx[xmax]
Ly = my[ymax]


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
tol       = 1e-10
ω         = 1.9
printint  = 1e6
method    = 3

## Source Term Settings
ζ = [sin(x/Lx*pi)^2 + sin(y/Ly*pi)^2 for x in mx, y in my] #.*0
contourf(mx,my,ζ)

## Boundary Conditions
# 1 = "Dirichlet", 2 = "Neumann", 3="Periodic"
# West
WBC = 3
wb_val = [y/y for y in mx]

# East
EBC = 3
eb_val = [y/y for y in my]

# North
NBC    = 3
nb_val = [x/x for x in mx]#[sin(3*(x/Lx*pi)) for x in mx]

# South
SBC    = 3
sb_val = [x/x for x in mx]#[sin(3*(x/Lx*pi)) for x in mx]


## Run the script

# Set periodic boundary options
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
ζ[1,:]    += By[1,:] # South BC
ζ[xmax,:] += By[2,:] # North BC
ζ[:,1]    += Bx[1,:] # West BC
ζ[:,ymax] += Bx[2,:] # East BC

# # note that this assumes xmax = ymax (equal sized spacing on x and y)
A = zeros(Float64,5,xmax)
A[1,:] = Cy[1,:]            # [  i , j-1]
A[2,:] = Cx[1,:]            # [i-1 ,   j]
A[3,:] = Cx[2,:] .+ Cy[2,:] # [  i ,   j]
A[4,:] = Cx[3,:]            # [i+1 ,   j]
A[5,:] = Cy[3,:]            # [  i , j+1]


S=ζ
ug = ones(Float64,xmax,ymax)

u_out,itcnt,r = ocnmod.FD_itrsolve_2D(Cx,Cy,S,ug,tol,ω,method,wper,eper,sper,nper,10^4)
contourf(mx,my,u_out)

# Plot Outside (Since BCs are wonky)
contourf(mx[2:end-1],my,u_out[2:end-1,:])
contourf(mx[2:end-1],my,u_out[2:end-1,:],
        xlims=(1, 5),
        ylims=(1, 5),
        )



# For Small Grids
# Currently it seems that z must by (y,x)
heatmap(mx,my,u_out)
heatmap(mx[2:end-1],my,u_out[:,2:end-1])

heatmap(ygr,xgr,ugr)
xgr= mx[2:end-1]
ygr= my
ugr = u_out[2:end-1,:]
