module ocnmod
    using Plots
    using Printf
    using LinearAlgebra
    """
    -----------------------------------------------------------
    FD_make_κ
    -----------------------------------------------------------
        Create vertical profile of eddy diffusivity (κ) based on
        2 layers (mixed layer with higher κ and interior with lower
        κ). κ0 is included in the interior layers (extra entry for
        κ_mld)

        Eventually: Implement "transition" layer rather than
        discontinuous diffusivity values

        Inputs
            1) mld     = Mixed Layer Depth
            2) levels  = Depth levels
            3) kmax    = Number of cells
            4) κ_mld   = Eddy Diffusivity in mixed-layer
            5) κ_int   = Eddy Diffusivity in the interior
        Outputs
            1) κ       = Vertical diffusivity profile

        """
    function FD_make_κ(mld,levels,kmax,κ_mld,κ_int)
        # Find # of levels deeper than mld
        #mldout    = argmin(levels.-mld)
        mldout    = length(findall(x->x>mld,levels))
        # Create κ profile [κ0 as part of interior, κ_int]
        κ         = hcat(ones(Float64,1,mldout + 1)    * κ_int,
                         ones(Float64,1,kmax - mldout) * κ_mld)
       return κ
    end

    """
    -----------------------------------------------------------
    FD_calc_I
    -----------------------------------------------------------
        Calculate the source/sink term for incoming shortwave radiation
        Inputs
            1) levels = depth levels
            2) z_f    = height of each cell
            3) z_att  = attentuation depth
            4) S0     = Initial incoming solar radiation
            5) cp0  (Default set to 3850)
            6) rho  (Default set to 1025)
        Outputs
            1) S      = Source/Sink Term
        """
    function FD_calc_I(levels,z_f,z_att,S0,ocn_trns,rho,cp0,mld)
        #Q = (ocn_trns * S0 * exp.(-1 * levels ./ z_att))
        #S = Q ./ (z_f .* cp0 .* rho)'
        S = (ocn_trns * S0 * exp.(-1 * levels ./ z_att)) ./ (mld .* cp0 .* rho)'
        return S
    end

    """
    -----------------------------------------------------------
    FD_calc_coeff
    -----------------------------------------------------------
        Compute both C (Matrix of Coefficients) and corresponding
        forcing terms B
        Inputs:
            ||~ Box Geometry and Indices ~||
            kmax   = maximum iteration size
            z_f    = Height of each grid cell
            z_c    = Midpoint distance between cells

            ||~ Initial Profiles         ~||
            κ      = Diffusivity Profile
            S  = Initial forcing/sink term (initial)

            ||~ Boundary Conditions      ~||
            BC_top  = 1 for Dirichlet, 2 for Neumann
            val_top = Input value for top
            BC_bot  = 1 for Dirichlet, 2 for Neumann
            val_bot = Input value for bot
            z_t     = Height to top cell
            z_b     = Height to bottom cell
        Outputs:
            1) C - Matrix of Coefficients [1:3,kmax]
                1 = C(k,k-1); 2 = C(k,k); 3 = C(k,k+1)
            2) B - Column vector of Source/Sink Terms, with
                first and last entries modified for BCs
        """
    function FD_calc_coeff(kmax,z_f,z_c, # Cell geometry/indices
        κ,S,                             # κ and Source/Sink Terms
        BC_top,val_top,                  # Bottom BCs
        BC_bot,val_bot,                  # Top BCs
        z_t,z_b)                         # Top and Bottom heights

        # Preallocate
        C = zeros(Float64,3,kmax)
        B = zeros(Float64,kmax)
        A = zeros(Float64,kmax,kmax)

        # Options for bottom boundary ------------------------------------
        k = 1

        # Dirichlet Bottom
        if BC_bot == 1
            # Compute Source Term with BC (Prescribed Temp)
            B[1]   = S[k] - val_bot *
                     ( (z_b * κ[k-1]) / (z_f[k] * z_c[k-1]) )

            # Ck0 remains the same. Multiplier goes on bottom term
            C[2,1] = (( 2 * κ[k-1]   / z_c[k-1]) +
                      (     κ[k]     / z_c[k]  )) * -1/z_f[k] # Prescribed Temp
        # Newmann Bottom
        elseif BC_bot == 2
            # Compute Source Term with BC (Prescribed Flux)
            global B[1]   = S[k] - val_bot / z_f[k]

            # Ck2 remains the same, Ck1 dep. on BCs
            C[2,1] = κ[k] / (z_f[k] * z_c[k])
        end

        # Set C(k,k-1) to 0, and set C(k,k+1) to the usual
        C[1,1] = 0
        C[3,1] = κ[k]   / (z_f[k] * z_c[k])

        # Options for top boundary ----------------------------------------
        k = kmax

        # Dirichlet Top
        if BC_top == 1

            # Compute Source Term with BC (Prescribed Temp)
            B[kmax]   = S[k] - val_top * ( (z_t * κ[k]) / (z_f[k] * z_c[k]) )

            # Calculate C(k,k)
            C[2,kmax] = (( 2 * κ[k]   / z_c[k]  ) +
                         (     κ[k-1] / z_c[k-1])) * -1/z_f[k] # Prescribed Temp

        # Neumann Top
        elseif BC_top == 2
            # Compute Source Term with BC (Prescribed Flux)
            global B[kmax]   = S[k] - val_top / z_f[k]

            # Calculate C(k,k)
            C[2,kmax] = κ[k-1] / (z_f[k] * z_c[k-1])# This depends on the type of BC

        end

        # Set C(k,k+1) to 0, and set C(k,k-1) to the usual
        C[1,kmax] = κ[k-1]     / (z_f[k] * z_c[k-1])
        C[3,kmax] = 0

        # Options for interior ----------------------------------------------------
        for k = 2:kmax-1
            B[k] = S[k]
            # Compute Coefficients
            C[1,k] = κ[k-1]     / (z_f[k] * z_c[k-1]  )
            C[3,k] = κ[k]   / (z_f[k] * z_c[k])
            C[2,k] = (C[1,k] + C[3,k]) * -1 # Note this might only work in certain cases
        end

        # Output A Matrix (optional)
        du = C[3,1:end-1]
        dl = C[1,2:end]
        d  = C[2,:]
        A  = Tridiagonal(dl,d,du)

        return C,B,A
    end

    """
    -----------------------------------------------------------
    FD_calc_T
    -----------------------------------------------------------
        For a steady state, 1-D Diffusion Problem
        Solves Ax = B using one of the following iterative methods
        (1) Jacobi (2) Gauss-Siedel (3) Successive Overrelaxation,
        where A is the sparse matrix of coefficients given by C,
        x is the temperature profile and B is the source term.

        Begins with an initial guess, x and iterates until the residual
        (sum of(Ax-B)^2) is below the tolerance level

        Note that input 8 is temporary (need to figure out correct
        way to compute residuals)

        Input
            1) C       - coefficients
            2) B       - source term
            3) x_g     - array of initial guesses
            4) tol     - tolerance level for residual
            5) Method  - (1) Jacobi, (2) Gauss-Siedel, (3) SOR
            6) max_iter - maximum number of iterations
            7) ω       - weights for successive overrelaxation
            8) x_inv   - Inverted solution for error calculation
        """
    function FD_calc_T(kmax,C,B,x_g,tol,Method,max_iter,ω,printint,A_in)
        allstart = time()

        # Preallocate Arrays [1 x # of cells]
        Tz  = zeros(Float64,2,length(x_g)) # Array for current profile (itr=m)

        # Input Initial Guess(2nd row of Tz)
        Tz[2,:] = x_g

        r   = Float64[]#zeros(Float64,max_iter+1)

        #Tz_all = zeros(Float64,max_iter,length(x_g))

        # Set up counters and indexing
        itcnt = 0     # Count of Iterations
        m     = 0     # Iteration Storage Number
        err   = tol+1 # initial value, pre error calculation (need to set up for the diff between guess and actual)

        while  err > tol #m <= max_iter #|| err > tol
            start = time()
            #@printf("Starting iter %i...",m)
            r_m  = zeros(Float64,1,kmax+1)
            for k = 1:kmax

                # First Cell
                if k == 1
                    T0 = 0
                    T2 = Tz[2,k+1] # Take T2 from last iteration (or init guess)

                # Last Cell
                elseif k == kmax
                    # For Jacobi Method, use previous values (Row 2 of Tz)
                    if Method == 1
                        T0 = Tz[2,k-1]

                    # For GS and SOR, use new value (Row 1 of Tz)
                    else
                        T0 = Tz[1,k-1]
                    end

                    T2 = 0

                # Interior Cells
                else
                    # For Jacobi Method, use previous values (row 2 Tz)
                    if Method == 1
                        T0 = Tz[2,k-1]
                    # For GS and SOR, use new value (row 1 Tz)
                    else
                        T0 = Tz[1,k-1]
                    end

                    # Use old T2 value for all cases
                    T2 = Tz[2,k+1]
                end

                # Get Coefficients and Source Term
                # C[Coeff#,level]
                C0 = C[1,k]
                C1 = C[2,k]
                C2 = C[3,k]
                B1 = B[k]

                # Compute Temperature
                T1 = 1/C1 * (B1 - C2*T2 - C0*T0)


                # need to compute this outside (with Tm+1 rather than with components)
                r_m[k] = C0*T0+C1*T1+C2*T2 - B1

                # Apply weights for sucessive overrelaxation
                if Method == 3 && itcnt > 0
                    global Tz[1,k] = ω * T1 + (1-ω) * Tz[2,k]
                else
                    global Tz[1,k] = T1
                end

                #Residuals
                #r_m[k] = (C0*T0 + C1*Tz[m,k] + C2*T2 - B1).^2
                if mod(m,printint)==0
                    #@printf("\tIT: %i k=%i with err %e\n",m,k,r_m[k])
                end
            end

            itcnt += 1
            m += 1

            # Cora's Method
            #r   = A_in*Tz[1,:]-B
            #err = norm(r)
            #err = sqrt(sum(r_m.^2)/length(r))

            # Tridiag method
            du = C[3,1:end-1] # Upper diagonal
            dl = C[1,2:end]   # Lower diagonal
            d  = C[2,:]       # Diagonal

            C_tri = Tridiagonal(dl,d,du)

            # Calculate residual (can do this more efficiently iteratively)
            err = norm(C_tri*Tz[1,:] - B)

            # Compute residual but differencing from last iteration
            #err = maximum(broadcast(sqrt,abs2.((Tz[1,:]-Tz[2,:])./Tz[1,:])))
            #err = sqrt(sum(norm(Tz[1,:]-Tz[2,:])))#./Tz[1,:])))
            #err = norm(Tz[1,:]-Tz[2,:])/length(Tz)
            #r[m] = sum(abs.((Tz[1,:]-Tz[2,:])))#./Tz)))
            #err = maximum(abs.((Tz[1,:]-Tz[2,:])))#./Tz)))

            #r[m] = sum(r_m,dims=2)[1]
            #x_itr = Tz[1,:]
            #x_err = sqrt.((x_itr' .- x_inv).^2)
            #r[m] =sum(x_err,dims=2)[1]
            #err = sum(x_err,dims=2)[1]

            # Push calculated error to residual variable
            push!(r,err)

            # Save new Tz profile to Tz0
            Tz[2,:] = Tz[1,:] # Note that index 2 is the prev

            # Save Tz profile to corresponding index on storage
            #Tz_all[m,:] = Tz[1,:]

            #@printf("Now on iteration %i with residual %.5E \n",m,err)

            elapsed = time() - allstart
            if mod(m,printint)==0
                @printf("Now on iteration %i with resid %.2E in %fs\n",m,err,elapsed)
                #@printf("\tTz ")
                #print(Tz)
                #@printf("\n\tTz0 ")
                #print(Tz0)

                #@printf("Now on iteration %i with err %f\n",m,err)
            end

        end
        # Restrict to iteration
        #r  =  r[1:itcnt]
        #Tz = Tz[1:itcnt,:]
        elapsed = time()-allstart
        @printf("Completed in %i iterations with resid %.2E. Only took %fs\n",itcnt,err,elapsed)
        return Tz, itcnt,r#, Tz_all
    end

    """
    -----------------------------------------------------------
    FD_inv_sol
    -----------------------------------------------------------
        Find the solution of Ax=B by inverting A, which contains
        the coefficients of a tridiagonal matrix.
        The first dim of a indicates...
            1 = subdiagonal (lower)
            2 = diagonal
            3 = superdiagonal (upper)

        Inputs
            1) A       = Coefficients (coeff#,k)
            2) B       = R.H.S. (1xk)

        Outputs
            1) A_tri   = Tridiagonal Matrix


        """
    function FD_inv_sol(A,B)
        # Create Tridiagonal Matrix
        du = A[3,1:end-1] # Upper diagonal
        dl = A[1,2:end]   # Lower diagonal
        d  = A[2,:]       # Diagonal
        A_tri = Tridiagonal(dl,d,du)
        x = B' * inv(A_tri)
       return A_tri,x
    end
    """
    -----------------------------------------------------------
    makeTridiag
    -----------------------------------------------------------
        Convert C[3,:] to a sparse tridiagonal matrix where:
            C[1,:] = Lower Diagonal
            C[2,:] = Diagonal
            C[3,:] = Upper Diagonal

        Inputs:
            1) C[3,:] = 3 x n array of values

        Output
            1) C_tri  = n x n Tridiagonal matrix

        """
    function makeTridiag(C)

        # Tridiag method
        du = C[3,1:end-1] # Upper diagonal
        dl = C[1,2:end]   # Lower diagonal
        d  = C[2,:]       # Diagonal

        # Make Tridiagonal
        C_tri = Tridiagonal(dl,d,du)
        return C_tri
    end

    """
    -----------------------------------------------------------
    CN_make_matrix
    -----------------------------------------------------------
        Given the finite difference form of Crank-Nicolson:

            u[n+1] - Δt/θ*(B*u[n+1] + f[n+1]) = u[n] + Δt/(θ-1)*(C*u[n] + f[n])

        Rewrites the Equation to the form:

            B'*u[n+1] = b

        Where the prime ' indicates the term has been "premultiplied" by the
        corresponding weight. (ex. B' = B*(-Δt/θ))

            b = C'*u[n] + f'[n] + f'[n+1]

        This ultimately prepares the matrices to enter into a iterative solver

        Inputs
            1) Δt = Timestep size
            2) θ  = Weight applied to n+1 step [θ=0.5, CN; =1,BWEuler; =0,FWEuler]
            3) IC = Initial Conditions
            4) C  = Matrix of Coefficients on r.h.s. (step n)
            5) B  = Matrix of Coefficients on l.h.s. (step n+1)
            6) fr = Forcing/Source Term on r.h.s. (step n  )
            7) fl = Forcing/Source Term on l.h.s. (step n+1)
        Output
            1) b  = premultiplied "Forcing Term" that combines Inputs(3,4,6,7)
            2) A  = premultiplied B matrix
        """
    function CN_make_matrix(Δt,θ,IC,C,B,fr,fl,meth)

        start = time()

        # Determine LHS and RHS multipliers
        # LHS - For timestep (n+1), multiply by θ
        l_mult =  Δt*(θ)
        # RHS - For timestep (n)  , multiply by 1-θ
        r_mult =  Δt*(1-θ)

        # Meth1: Add Timestep corrections first
        if meth == 1
            C[2,:] = C[2,:] .+ (1/r_mult)
            B[2,:] = B[2,:] .- (1/l_mult)
        end

        # Multiply variables by time and theta factors
        B      = B  .* l_mult
        fl     = fl .* l_mult
        C      = C  .* r_mult
        fr     = fr .* r_mult

        # Meth2: Add single digit post-multiplication
        # To diagonal (Ck,Bk)
        if meth == 2
            C[2,:] = C[2,:] .+ 1
            B[2,:] = B[2,:] .+ 1
        end

        # Now combine terms
        C_tri = makeTridiag(C)
        b     = C_tri * IC + fr + fl # Bring Sl from left side
        A     = -B

        elapsed = time() - start
        #@printf("Completed matrix prep in %fs\n",elapsed)
        return A,b
    end

        #
        #
        #
        # """
        # -----------------------------------------------------------
        # CN_solver
        # -----------------------------------------------------------
        #     For a steady state, 1-D Diffusion Problem
        #     Solves Ax = B using one of the following iterative methods
        #     (1) Jacobi (2) Gauss-Siedel (3) Successive Overrelaxation,
        #     where A is the sparse matrix of coefficients given by C,
        #     x is the temperature profile and B is the source term.
        #
        #     Begins with an initial guess, x and iterates until the residual
        #     (sum of(Ax-B)^2) is below the tolerance level
        #
        #     Note that input 8 is temporary (need to figure out correct
        #     way to compute residuals)
        #
        #     Input
        #         1) C       - coefficients
        #         2) B       - source term
        #         3) x_g     - array of initial guesses
        #         4) tol     - tolerance level for residual
        #         5) Method  - (1) Jacobi, (2) Gauss-Siedel, (3) SOR
        #         6) max_iter - maximum number of iterations
        #         7) ω       - weights for successive overrelaxation
        #         8) x_inv   - Inverted solution for error calculation
        #     """
        # function CN_solver(kmax,C,Sr,Sl,B,x_g,x_init,printint,Method,ω,θ)
        #     allstart = time()
        #
        #
        #
        #     # Preallocate Arrays [1 x # of cells]
        #     Tz  = zeros(Float64,2,kmax) # Array for current profile (itr=m)
        #
        #     # Input Initial Guess(2nd row of Tz)
        #     Tz[2,:] = x_g
        #
        #     r   = Float64[]#zeros(Float64,max_iter+1)
        #
        #     #Tz_all = zeros(Float64,max_iter,length(x_g))
        #
        #     # Set up counters and indexing
        #     itcnt = 0     # Count of Iterations
        #     m     = 0     # Iteration Storage Number
        #     err   = tol+1 # initial value, pre error calculation (need to set up for the diff between guess and actual)
        #
        #     while  err > tol #m <= max_iter #|| err > tol
        #         start = time()
        #         for k = 1:kmax
        #
        #             # First Cell
        #             if k == 1
        #                 # Set bottom flux to 0 (included with Forcing Term)
        #                 Tr0 = 0
        #                 Tl0 = 0
        #
        #                 # For top flux, rhs = IC, lhs = init. guess
        #                 Tr2 = x_init[k+1]
        #                 Tl2 = Tz[2,k+1]
        #
        #             # Last Cell
        #             elseif k == kmax
        #
        #                 # For Jacobi Method, use previous values (Row 2 of Tz)
        #                 if Method == 1
        #                     Tr0 =
        #                     Tl0 = Tz[2,k-1]
        #
        #                 # For GS and SOR, use new value (Row 1 of Tz)
        #                 else
        #                     T0 = Tz[1,k-1]
        #                 end
        #
        #                 T2 = 0
        #
        #             # Interior Cells
        #             else
        #                 # For Jacobi Method, use previous values (row 2 Tz)
        #                 if Method == 1
        #                     T0 = Tz[2,k-1]
        #                 # For GS and SOR, use new value (row 1 Tz)
        #                 else
        #                     T0 = Tz[1,k-1]
        #                 end
        #
        #                 # Use old T2 value for all cases
        #                 T2 = Tz[2,k+1]
        #             end
        #
        #             # RHS Variables
        #             Tr0 = x_init[k-1] # Need to adjust for the first and last iteration
        #             Tr1 = x_init[k]
        #             Tr2 = x_init[k+1]
        #             Srr = Sr[k] #* ((1-θ)/Δt)
        #             C0  = C[1,k]
        #             C1  = C[2,k] + 1#((1-θ)/Δt) # Premultiplied C with Δt/(1-θ)
        #             C2  = C[3,k]
        #
        #             # LHS Variables
        #             if itcnt == 0
        #                 Tl0 = x_g[k-1]
        #                 Tl1 = x_g[k]
        #                 Tl2 = x_g[k+1]
        #             else
        #                 Tl0 = Tz[2,k-1]
        #                 #Tl1 = Tz[2,k]
        #                 Tl2 = Tz[2,k+1]
        #             end
        #             Sll= Sl[k] #* (θ/Δt)
        #             B0 = B[1,k]
        #             B1 = B[2,k]  -1 #- (θ/Δt) # Premultiplied B with Δt/(θ)
        #             B2 = B[3,k]
        #
        #             # Compute the temperature at the grid cell
        #             rhs = C0*Tr0 + C1*Tr1 + C2*Tr2 + Srr
        #             Tl1 = (1/B1) * (rhs - B0*Tl0 - B2*Tl2 - Sll) # CHECK SLL term
        #
        #             # Apply weights for sucessive overrelaxation
        #             if Method == 3 && itcnt > 0
        #                 global Tz[1,k] = ω * Tl1 + (1-ω) * Tz[2,k]
        #             else
        #                 global Tz[1,k] = Tl1
        #             end
        #
        #         end
        #
        #         itcnt += 1
        #         m += 1
        #
        #         # Cora's Method
        #         #r   = A_in*Tz[1,:]-B
        #         #err = norm(r)
        #         #err = sqrt(sum(r_m.^2)/length(r))
        #
        #         # Tridiag method
        #         du = C[3,1:end-1] # Upper diagonal
        #         dl = C[1,2:end]   # Lower diagonal
        #         d  = C[2,:]       # Diagonal
        #
        #         C_tri = Tridiagonal(dl,d,du)
        #
        #         # Calculate residual
        #         err = norm(C_tri*Tz[1,:] - B)
        #
        #         # Push calculated error to residual variable
        #         push!(r,err)
        #
        #         # Save new Tz profile to Tz0
        #         Tz[2,:] = Tz[1,:] # Note that index 2 is the prev
        #
        #         # Save Tz profile to corresponding index on storage
        #         #Tz_all[m,:] = Tz[1,:]
        #
        #         #@printf("Now on iteration %i with residual %.5E \n",m,err)
        #
        #         elapsed = time() - allstart
        #         if mod(m,printint)==0
        #             @printf("Now on iteration %i with resid %.2E in %fs\n",m,err,elapsed)
        #             #@printf("\tTz ")
        #             #print(Tz)
        #             #@printf("\n\tTz0 ")
        #             #print(Tz0)
        #
        #             #@printf("Now on iteration %i with err %f\n",m,err)
        #         end
        #
        #     end
        #
        #
        #
        #     # Restrict to iteration
        #     #r  =  r[1:itcnt]
        #     #Tz = Tz[1:itcnt,:]
        #     elapsed = time()-allstart
        #     @printf("Completed in %i iterations with resid %.2E. Only took %fs\n",itcnt,err,elapsed)
        #     return Tz, itcnt,r#, Tz_all
        # end
        #



end
