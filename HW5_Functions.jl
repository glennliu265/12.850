using Plots
using LinearAlgebra

include("AllFunctions12850.jl")

x0      = 0
xn      = 10
dx      = 1

# Time
dt    = 0.1
tsmax = 100

# Grid Setup (u as edges, t as midpoints)
ugrid  = [x0:dx:xn;]
tgrid  = ocnmod.get_midpoints(ugrid)
xmax   = length(tgrid)

# Prescribe velocity field
u      = ones(size(tgrid))*1

# Prescribe tracer field
c0      = ones(size(tgrid)) * 5
mask    = 2 .< tgrid .< 5
c0 = c0 .* mask
plot(tgrid,c0)

# Evaluate stability
cour = u[1] * dt / dx


# FEM and UW1
c = zeros(Float64,xmax,tsmax)
for t = 1:tsmax
    dt = 1

    if t == 1
        cin = c0
    else
        cin = c[:,t-1]
    end

    ct = zeros(Float64,size(tgrid))
    for i = 1:xmax
        ui = u[i]
        ci = cin

        # Assume Periodic
        if i == 1
            ip1 = i+1
            im1 = xmax
        elseif i == xmax
            ip1 = 1
            im1 = i-1
        else
            ip1 = i+1
            im1 = i-1
        end

        # Calculate u+ and u-
        up = (ui + abs(ui))/2
        um = (ui - abs(ui))/2

        ct[i] = cin[i] - dt/dx * (up*(ci[i] - ci[im1]) + um*(ci[ip1]-ci[i]))

    end

    c[:,t] = ct
end


animtest = @animate for t ∈ 1:tsmax
    plot(tgrid,c[:,t],
        title ="t="*string(t)
        )
        #ylims=(0,5))

end
gif(animtest,"HW5_test.gif",fps=5)
