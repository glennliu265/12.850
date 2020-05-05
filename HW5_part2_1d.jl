# Testing

using Plots
using LinearAlgebra
using Printf
include("AllFunctions12850.jl")
include("ocn5.jl")

# Set up grid
xgrid = [0:1:100;]

# Get Midpoints
mx      = ocnmod.get_midpoints(xgrid)

# Get Cell Width
δx        = (xgrid[2]-xgrid[1])


# Get End value and indices
xmax        = length(mx)
Lx          = mx[xmax]

# Get Cell Spacing and Widths (Currently Fixed)
x_f       = ones(Int8,1,xmax)*δx   # Cell Face Length (x)

# Velocity Field
u = [1 for x in xgrid, t in 1:ts_max]


# Set BCs
WBC = 3
wb_val = [0]#[y/y for y in mx]

# Eastern
EBC = 3
eb_val = [0]#[y/y for y in my]

chk_per = [1,1]

## Forcing

S = [0 for x in mx, t in 1:ts_max]



## Time parameters
dt        = 1    # Timestep for model integration
ts_max    = 25         # Number of timesteps to take
θ         = 0


## IC

T0 = 5*ones(Float64,xmax)
mask    = 20 .< mx .< 50
T0 = T0 .* mask
plot(mx,T0)

## Script start?

# Preallocate
temps = zeros(Float64,xmax,ts_max)
tinv = zeros(Float64,xmax,ts_max)


for t = 1:ts_max

        if t == ts_max
                tp1 = 1
        else
                tp1 = t+1
        end

        if t == 1
                IC = T0
        else
                IC = tinv[:,t-1]
        end

        # Get Corresponding forcing term
        S0 = S[:,t]
        S1 = S[:,tp1]

        # Get Coefficients (UW1) and modifications to forcing term
        C0,B0 = ocn5.UW_calc_coeff_1D(x_f,xmax,u[:,1],EBC,eb_val,WBC,wb_val)
        C1,B1 = ocn5.UW_calc_coeff_1D(x_f,xmax,u[:,tp1],EBC,eb_val,WBC,wb_val)

        # Modify Forcing term
        S0 += B0
        S1 += B1

        # Set up CN Matrices
        A,b = ocn5.CN_make_matrix_1D(dt,θ,IC,C0,C1,S0,S1,x_f,chk_per)

        # Solve matrix by inversion
                # Tridiag method
        du = A[3,1:end-1] # Upper diagonal
        dl = A[1,2:end]   # Lower diagonal
        d  = A[2,:]       # Diagonal

        A_tri = Tridiagonal(dl,d,du)

        # Save solution by inversion
        tinv[:,t] = b' * inv(A_tri)


end

animtest = @animate for t ∈ 1:ts_max
    #if t%2 == 0
        plot(mx,tinv[:,t],
            title ="t="*string(t),
            ylims=(-6,6)
            )
        #end

end
gif(animtest,"HW5_1D_uw1_type.gif",fps=5)
