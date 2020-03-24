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
ymax       = length(my)
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
κx        = ones(Float16,1,xmax)
κy        = ones(Float16,1,ymax)

κx0       = κx[1]
κy0       = κy[1]

## Iteration Parameters ------------------------------
tol       = 1e-4
ω         = 1.9
printint  = 1e6
method    = 3

## Source Term Settings
ζ = [sin(x/Lx*pi)^2 + sin(y/Ly*pi)^2 for x in mx, y in my]
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
NBC    = 1
nb_val = [sin(3*(x/Lx*pi)) for x in mx]

# South
SBC = 1
sb_val = [sin(3*(x/Lx*pi)) for x in mx]


## Run the script


# Compute Coefficients and modifications to BCs
Cx,Bx = ocnmod.FD_calc_coeff_2D(xmax,x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
Cy,By = ocnmod.FD_calc_coeff_2D(ymax,y_f,y_c,κy,NBC,nb_val,SBC,sb_val,y_c0,κy0)

# Modify Boundary Conditions
ζ[1,:]    += By[1,:] # South BC
ζ[xmax,:] += By[2,:] # North BC
ζ[:,1]    += Bx[1,:] # West BC
ζ[:,ymax] += Bx[2,:] # East BC

# note that this assumes xmax = ymax (equal sized spacing on x and y)
A = zeros(Float64,5,xmax)
A[1,:] = Cy[1,:]            # [  i , j-1]
A[2,:] = Cx[1,:]            # [i-1 ,   j]
A[3,:] = Cx[2,:] .+ Cy[2,:] # [  i ,   j]
A[4,:] = Cx[3,:]            # [i+1 ,   j]
A[5,:] = Cy[3,:]            # [  i , j+1]

for j = 1:ymax

    # Get Coefficients (y)
    B1  = Cy[1,j]
    B3y = Cy[2,j]
    B5  = Cy[3,j]

    for i = 1:xmax

        # Get Coefficients (x)
        B2 = Cx[1,i]
        B3 = B3y + Cx[2,i]
        B4 = Cx[2,i]

        f = S()





for y = 1:ymax

    for x = 1:xmax



    end

end

# Modify Forcing Term and combine other terms



#Cmake-> run for x, then for y,putting same sourcing term in both
# combine the middle terms

# X conditionals
if x == 0

elseif x == Lx

else

end

# Y conditionals
if y == 0

elseif y == Ly

else

end
