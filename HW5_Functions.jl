using Plots
using LinearAlgebra

include("AllFunctions12850.jl")

#immpt
dx      = 1
dt      = .1
u_value = .1

# Grid
x0      = 0
xn      = 100


# Time
tsmax = 100

# Question 1 type
q1t = 1


# Grid Setup (u as edges, t as midpoints)
ugrid  = [x0:dx:xn;]
tgrid  = ocnmod.get_midpoints(ugrid)
xmax   = length(tgrid)

# Prescribe velocity field
u      = ones(size(tgrid))*u_value

# Prescribe tracer field
if q1t == 1
    c0      = ones(size(tgrid)) * 5
    mask    = 2 .< tgrid .< 5
    c0 = c0 .* mask
elseif q1t == 2
    c0 = sin.(pi.* tgrid ./ tgrid[xmax])
end

plot(tgrid,c0)


## FEM - UW1
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
        title ="t="*string(t),
        ylims=(0,6))

end
gif(animtest,"HW5_test_uw1_type"*string(q1t)*".gif",fps=5)



## FEM - UW2
# Evaluate stability
#cour = u[1] * dt / dx

# FEM and UW2
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
            ip2 = i+2
            im1 = xmax
            im2 = xmax-1

        elseif i == 2

            ip1 = i+1
            ip2 = i+2
            im1 = i-1
            im2 = xmax

        elseif i == xmax

            ip1 = 1
            ip2 = 2
            im1 = i-1
            im2 = i-2

        elseif i == xmax-1

            ip1 = i+1
            ip2 = 1
            im1 = i-1
            im2 = i-2

        else

            ip1 = i+1
            ip2 = i+2
            im1 = i-1
            im2 = i-2

        end

        # Calculate u+ and u-
        up = (ui + abs(ui))/2
        um = (ui - abs(ui))/2



        ct[i] = ci[i] - dt/(dx*2) * (up* (3*ci[i] - 4*ci[im1] + 1*ci[im2])
                + um*(-1*ci[ip2]+4*ci[ip1]-3*ci[i]) )

    end

    c[:,t] = ct
end


animtest = @animate for t ∈ 1:tsmax
    plot(tgrid,c[:,t],
        title ="t="*string(t),
        ylims=(0,6)
        )

end
gif(animtest,"HW5_test_uw1_type"*string(q1t)*".gif",fps=5)





# Amplification Factor Modulus
# FEM-UW1
# Courant Number
nu = u[1]*dt/dx

# Create k vector
kdx = [0:(pi/12):pi;]


a2 = (1 .- nu .+ nu .* cos.(kdx)).^2
b2 = (nu .* sin.(kdx)).^2
lambda = sqrt.(a2.+b2)

plot(kdx,lambda,proj=:polar,lims=(0,pi))
