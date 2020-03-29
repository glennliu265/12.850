module ocnmod
    using Plots
    using Printf
    using LinearAlgebra

    """
    -----------------------------------------------------------
    get_midpoints
    -----------------------------------------------------------
        Get midpoints for an array with the height of each cell face

        Inputs
            1) levels   = contains height at N levels
        Outputs
            1) midpoints = height of each midpoint (N-1 points)

        """
    function get_midpoints(levels)
        numpts = length(levels)-1
        midpoints = zeros(Float64,numpts)
        for k = 1:numpts
            midpoints[k] = levels[k] + (levels[k+1]-levels[k])/2
        end
       return midpoints
    end


        """
        -----------------------------------------------------------
        FD_calc_coeff_2D
        -----------------------------------------------------------
            Compute both C (Matrix of Coefficients) and modification to
            corresponding forcing term

            NOTE: Assumes Forcing term is on the same side of the equation
            as the diffusion term.

            Inputs:
                ||~ Box Geometry and Indices ~||
                kmax   = maximum iteration size
                z_f    = Height of each grid cell
                z_c    = Midpoint distance between cells
                z_c0   = Midpoint distance to zeroth cell

                ||~ Initial Profiles         ~||
                κ      = Diffusivity Profile
                κ0     = Diffusivity for zeroth cell

                ||~ Boundary Conditions      ~||
                BC_top  = 1 for Dirichlet, 2 for Neumann
                val_top = Input value for top (matches size of other dim)
                BC_bot  = 1 for Dirichlet, 2 for Neumann
                val_bot = Input value for bot (matches size of other dim)

            Outputs:
                1) C - Matrix of Coefficients [1:3,kmax]
                    1 = C(k,k-1); 2 = C(k,k); 3 = C(k,k+1)
                2) B - Column vector of modification to Source/Sink Terms, with
                    first and last entries modified for BCs along the other dimension
            """
        function FD_calc_coeff_2D(kmax,z_f,z_c, # Cell geometry/indices
            κ,                             # κ and Source/Sink Terms
            BC_top,val_top,                  # Bottom BCs
            BC_bot,val_bot,                  # Top BCs
            z_c0,κ0)                         # bottom diffusivity and midpoint dist

            #@printf("Applying the following BCs:")
            #@printf("\n\tBottom - Type %i of value %.2E",BC_bot,val_bot)
            #@printf("\n\tTop    - Type %i of value %.2E",BC_top,val_top)

            # Preallocate C: length of dimension
            # B: length of boundary
            C = zeros(Float64,3,kmax)
            B = zeros(Float64,2,length(val_top)) # 1st dim, 1=bot,2=top,


            # Options for bottom boundary ------------------------------------
            k = 1

            # Dirichlet Bottom
            if BC_bot == 1
                # Compute Source Term with BC (Prescribed Temp)
                B[1,:]   = val_bot .* (2 * κ0) / (z_f[k] * z_c0)

                # Ck0 remains the same. Multiplier goes on bottom term
                C[2,1] = (( 2 * κ0   / z_c0    ) +
                          (     κ[k] / z_c[k]  )) * -1/z_f[k] # Prescribed Temp
            # Newmann Bottom
            elseif BC_bot == 2
                # Compute Source Term with BC (Prescribed Flux)
                B[1,:]   = val_bot ./ z_f[k]

                # Ck2 remains the same, Ck1 dep. on BCs Need negative here!
                C[2,1] = - κ[k] / (z_f[k] * z_c[k])


            #Periodic Bottom
            elseif BC_top == 3
                C[1,k] = κ[kmax]     / (z_f[k] * z_c[kmax]  )
                C[3,k] = κ[k]       / (z_f[k] * z_c[k])
                C[2,k] = (C[1,k] + C[3,k]) * -1
            end

            # Set C(k,k-1) to 0, and set C(k,k+1) to the usual
            C[1,1] = 0
            C[3,1] = κ[k]   / (z_f[k] * z_c[k])

            # Options for top boundary ----------------------------------------
            k = kmax

            # Dirichlet Top
            if BC_top == 1

                # Compute Change to Source Term with BC (Prescribed Temp)
                B[2,:]   = val_top .* ( 2 * κ[k]) / (z_f[k] * z_c[k])

                # Calculate C(k,k)
                C[2,kmax] = (( 2 * κ[k]   / z_c[k]  ) +
                             (     κ[k-1] / z_c[k-1])) * -1/z_f[k] # Prescribed Temp

            # Neumann Top
            elseif BC_top == 2
                # Compute Source Term with BC (Prescribed Flux)
                B[2,:]   = val_top ./ z_f[k]

                # Calculate C(k,k) (need negative here!)
                C[2,kmax] = -κ[k-1] / (z_f[k] * z_c[k-1])# This depends on the type of BC

            #Periodic Top
            elseif BC_top == 3
                C[1,k] = κ[k-1]     / (z_f[k] * z_c[k-1]  )
                C[3,k] = κ[1]       / (z_f[k] * z_c[1])
                C[2,k] = (C[1,k] + C[3,k]) * -1
            end

            # Set C(k,k+1) to 0, and set C(k,k-1) to the usual
            C[1,kmax] = κ[k-1]     / (z_f[k] * z_c[k-1])
            C[3,kmax] = 0

            # Options for interior ----------------------------------------------------
            for k = 2:kmax-1
                # Compute Coefficients
                C[1,k] = κ[k-1]     / (z_f[k] * z_c[k-1]  )
                C[3,k] = κ[k]       / (z_f[k] * z_c[k])
                C[2,k] = (C[1,k] + C[3,k]) * -1 # Note this might only work in certain cases
            end

            return C,B
        end

        """
        FD_itrsolve_2D
         Cx     = x coefficients (3 x n)
         Cy     = y coefficients (3 x m)
         S      = Modified Forcing Term (n x m)
         ug     = Initial guess at quantity u (n x m)
         tol    = Error Tolerance
         ω      = Weight for SOR
         method = [1,Jacobi] [2,Gauss-Siedel] [3,SOR]
         wper   = Periodic Western Boundary (1 = periodic)
         eper   = Periodic Eastern Boundary
         sper   = Periodic Southern Boundary
         nper   = Periodic Northern Boundary

         Out:
         u_out = Final approximation of quantity[n x m]
         itcnt = # of iterations to convergence
         r     = Array of residuals per iteration

        """
        function FD_itrsolve_2D(Cx,Cy,S,ug,tol,ω,method,wper,eper,sper,nper,maxiter,saveiter)
            xmax = size(Cx,2)
            ymax = size(Cy,2)

            # Preallocate
            u = zeros(Float64,2,xmax,ymax)
            r = Float64[]

            # Scrap (Delete Later: Save first 10 iterations)
            u_scrap = zeros(Float64,saveiter,xmax,ymax)
            err_scrap = zeros(Float64,saveiter,xmax,ymax)

            # Assign ug to the first entry
            u[1,:,:] = ug # 1 will store the previous guess
            u[2,:,:] = ug # 2 will strore the updated guess

            itcnt = 0
            start = time()
            while itcnt == 0 || r[itcnt] > tol
                # Loop for each column, by row
                for j = 1:ymax

                    # Get Coefficients (y)
                    B1  = Cy[1,j]
                    B3y = Cy[2,j]
                    B5  = Cy[3,j]

                    for i = 1:xmax

                        # Get Coefficients (x)
                        B2 = Cx[1,i]
                        B3 = B3y + Cx[2,i]
                        B4 = Cx[3,i]

                        # Retrieve value from Source Term
                        f = S[i,j]

                        # Set indexing for iteration method (applies for only u1,u2)
                        if method == 1
                            midx = 1 # Take previous term (Jacobi)
                        else
                            midx = 2 # Take newly updated term (Gauss Siedel and SOR)
                        end

                        # First, assume periodic boudaries
                        # Make i indices periodic
                        i2 = i-1
                        i4 = i+1
                        if i == 1
                            i2 = xmax
                        end
                        if i == xmax
                            i4 = 1
                        end
                        # Make j indices periodic
                        j1 = j-1
                        j5 = j+1
                        if j == ymax
                            j5 = 1
                        end

                        if j == 1
                            j1 = ymax
                        end

                        # Interior Points, Periodic
                        u1 = u[midx,i,j1] # Take j-1 depending on method
                        u2 = u[midx,i2,j] # Take i-1 depending on method
                        u4 = u[1   ,i4,j]
                        u5 = u[1   ,i,j5]

                        # Modifications for nonperiodic cases
                        if wper != 1 && i == 1
                            u2 = 0 # All i-1 terms = 0
                        end

                        if eper != 1 && i == xmax
                            u4 = 0 # All i+1 terms = 0
                        end

                        if sper != 1 && j == 1
                            u1 = 0 # All j-1 terms = 0
                        end

                        if nper != 1 && j == ymax
                            u5 = 0 # All j+1 terms = 0
                        end

                        # Now calculate the central point
                        u3 = (f-B1*u1-B2*u2-B4*u4-B5*u5)/B3

                        if method == 3
                            u[2,i,j] = ω * u3 + (1-ω) * u[1,i,j]
                        else
                            u[2,i,j] = u3
                        end

                    end
                # End loop for the point, i,j
                end

                ##  Repeat process to compute residual-------------
                err = zeros(Float64,xmax,ymax)

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

                        # Retrieve value from Source Term
                        f = S[i,j]

                        # First, assume periodic boudaries
                        # Make i indices periodic
                        i2 = i-1
                        i4 = i+1
                        if i == 1
                            i2 = xmax
                        end
                        if i == xmax
                            i4 = 1
                        end
                        # Make j indices periodic
                        j1 = j-1
                        j5 = j+1
                        if j == ymax
                            j5 = 1
                        end

                        if j == 1
                            j1 = ymax
                        end
                        # Interior Points, Periodic
                        u1 = u[2,i,j1]
                        u2 = u[2,i2,j]
                        u3 = u[2,i,j]
                        u4 = u[2,i4,j]
                        u5 = u[2,i,j5]

                        # Modifications for nonperiodic cases
                        if wper != 1 && i == 1
                            u2 = 0 # All i-1 terms = 0
                        end

                        if eper != 1 && i == xmax
                            u4 = 0 # All i+1 terms = 0
                        end

                        if sper != 1 && j == 1
                            u1 = 0 # All j-1 terms = 0
                        end

                        if nper != 1 && j == ymax
                            u5 = 0 # All j+1 terms = 0
                        end

                        err[i,j] = (B1*u1+B2*u2+B3*u3+B4*u4+B5*u5) - f
                    end
                end

                # Assign this iteration to the last
                u[1,:,:] = u[2,:,:]
                #u[2,:,:] .*= 0

                push!(r,norm(err))
                itcnt += 1
                if itcnt > maxiter
                    break
                end

                if itcnt%10^5 == 0
                    elapsed = time()-start
                    @printf("\nOn Iteration %i in %s s",itcnt,elapsed)
                end

                # Scrap: Save first 10 iterations and error
                if itcnt <=saveiter
                    u_scrap[itcnt,:,:]=u[1,:,:]
                    err_scrap[itcnt,:,:]=err
                end


            end
            u_out = u[1,:,:]

            elapsed = time()-start
            @printf("\nFinished %.2e iterations in %s",itcnt,elapsed)
            return u_out, itcnt, r, u_scrap, err_scrap

        end



    # Module End
    end
