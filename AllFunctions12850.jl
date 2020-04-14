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
        κ         = hcat(ones(Float64,1,mldout)        * κ_int,
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
    function FD_calc_I(levels,z_f,z_att,S0,ocn_trns,rho,cp0)
        #Q = (ocn_trns * S0 * exp.(-1 * levels ./ z_att))
        #S = Q ./ (z_f .* cp0 .* rho)'
        S = (ocn_trns * S0 * exp.(-1 * levels ./ z_att)) ./ (z_att .* cp0 .* rho)'
        return S
    end

    """
    -----------------------------------------------------------
    FD_calc_coeff
    -----------------------------------------------------------
        Compute both C (Matrix of Coefficients) and corresponding
        forcing terms B

        NOTE: Assumes Forcing term is on the same side of the equation
        as the diffusion term.

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
        z_t,z_b,                         # Dist to top/bot face from midpt
        z_c0,κ0)                         # bottom diffusivity and midpoint dist

        #@printf("Applying the following BCs:")
        #@printf("\n\tBottom - Type %i of value %.2E",BC_bot,val_bot)
        #@printf("\n\tTop    - Type %i of value %.2E",BC_top,val_top)
        # Preallocate
        C = zeros(Float64,3,kmax)
        B = zeros(Float64,kmax)
        A = zeros(Float64,kmax,kmax)

        # Options for bottom boundary ------------------------------------
        k = 1

        # Dirichlet Bottom
        if BC_bot == 1
            # Compute Source Term with BC (Prescribed Temp)
            B[1]   = S[k] + val_bot *
                     ( (z_b * κ0) / (z_f[k] * z_c0) )

            # Ck0 remains the same. Multiplier goes on bottom term
            C[2,1] = (( 2 * κ0   / z_c0    ) +
                      (     κ[k] / z_c[k]  )) * -1/z_f[k] # Prescribed Temp
        # Newmann Bottom
        elseif BC_bot == 2
            # Compute Source Term with BC (Prescribed Flux)
            B[1]   = S[k] + val_bot / z_f[k]

            # Ck2 remains the same, Ck1 dep. on BCs Need negative here!
            C[2,1] = - κ[k] / (z_f[k] * z_c[k])
        end

        # Set C(k,k-1) to 0, and set C(k,k+1) to the usual
        C[1,1] = 0
        C[3,1] = κ[k]   / (z_f[k] * z_c[k])

        # Options for top boundary ----------------------------------------
        k = kmax

        # Dirichlet Top
        if BC_top == 1

            # Compute Source Term with BC (Prescribed Temp)
            B[kmax]   = S[k] + val_top * ( (z_t * κ[k]) / (z_f[k] * z_c[k]) )

            # Calculate C(k,k)
            C[2,kmax] = (( 2 * κ[k]   / z_c[k]  ) +
                         (     κ[k-1] / z_c[k-1])) * -1/z_f[k] # Prescribed Temp

        # Neumann Top
        elseif BC_top == 2
            # Compute Source Term with BC (Prescribed Flux)
            B[kmax]   = S[k] + val_top / z_f[k]

            # Calculate C(k,k) (need negative here!)
            C[2,kmax]        = -κ[k-1] / (z_f[k] * z_c[k-1])# This depends on the type of BC

        end

        # Set C(k,k+1) to 0, and set C(k,k-1) to the usual
        C[1,kmax] = κ[k-1]     / (z_f[k] * z_c[k-1])
        C[3,kmax] = 0

        # Options for interior ----------------------------------------------------
        for k = 2:kmax-1
            B[k] = S[k]
            # Compute Coefficients
            C[1,k] = κ[k-1]     / (z_f[k] * z_c[k-1]  )
            C[3,k] = κ[k]       / (z_f[k] * z_c[k])
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
            #r_m  = zeros(Float64,1,kmax+1)
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
                #r_m[k] = C0*T0+C1*T1+C2*T2 - B1

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
        @printf("Completed in %i iterations with resid %.2E (%fs)\n",itcnt,err,elapsed)
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
    print_exeq
        Function to print a sample equation
    """
    function print_exeq(k,kmax,b,C,IC,fr,fl,A)
        if k == 1
            IC0 = 0
            IC2 = IC[k+1]
        elseif k == kmax
            IC0 = IC[k-1]
            IC2  = 0
        else
            IC0 = IC[k-1]
            IC2 = IC[k+1]
        end
        @printf("\n For this loop, we have for level %i of %i",k,kmax)
        @printf("\n\tlhs1: A Matrix:[%.2E, %.2E, %.2E]",A[1,k],A[2,k],A[3,k])
        @printf("\n\trhs2: C Matrix:[%.2E, %.2E, %.2E]",C[1,k],C[2,k],C[3,k])
        @printf("\n\trhs3: IC      :[%.2E, %.2E, %.2E]",IC0,IC[k],IC2)
        @printf("\n\trhs4: fr      :[%.2E]",fr[k])
        @printf("\n\trhs5: fl      :[%.2E]",fl[k])
        @printf("\n\trhs5: b       :[%.2E]",b[k])
        print("\nWhere we have: (1)*u = (2)*(3)+(4)+(5) and...")
        print("\n\tb =  (2)*(3)+(4)+(5)...")
        print("\n")
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
    function CN_make_matrix(Δt,θ,IC,C,B,fr,fl,meth,kprint,z_f)

        start = time()

        # Determine LHS and RHS multipliers
        # LHS - For timestep (n+1), multiply by θ
        l_mult =  -Δt*(θ)
        # RHS - For timestep (n)  , multiply by 1-θ
        r_mult =  Δt*(1-θ)

        # Meth1: Add Timestep corrections first
        # Note, negative signs for Bk, Ck already included
        if meth == 1
            C[2,:] = C[2,:] .+ abs((1/r_mult))
            B[2,:] = B[2,:] .- abs((1/l_mult))
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
        t0    = C_tri * IC
        b     = C_tri * IC + fr - fl # Bring Sl from left side
        A     = B

        k = kprint
        #print_exeq(k,length(fr),b,C,IC,fr,fl,A)

        elapsed = time() - start
        #@printf("Completed matrix prep in %fs\n",elapsed)
        return A,b,t0
    end

    """
    -----------------------------------------------------------
    # print_itrsolv
        Function to print a sample equation
    -----------------------------------------------------------
        Prints sample equation for a given k level and every
        itcnt iteration. Designed to run within the loop for
        FD_itrsolve.
    """
    function print_itrsolv(itcnt,A0,A1,A2,b1,x0,x1,x2,k,kmax)
        @printf("\nOn itr %i level %i of %i we have...",itcnt,k,kmax)
        @printf("\n\t%.2E = 1/%.2E * ( %.2E -...",x1,A1,b1)
        @printf("\n\t\t... %.2E * %.2E - %.2E * %.2E ) ",A0,x0,A2,x2)
        print("\n")
    end

    """
    -----------------------------------------------------------
    FD_itrsolve(A,b,x_init,tol,Method,ω)
    -----------------------------------------------------------
        For a steady state, 1-D Diffusion Problem
        Solves Ax = b using one of the following iterative methods
        (1) Jacobi (2) Gauss-Siedel (3) Successive Overrelaxation,
        where A is the sparse matrix of coefficients given by A,
        x is the temperature profile and B is the source term.

        Begins with an initial guess, x and iterates until the residual
        (sum of(Ax-B)^2) is below the tolerance level

        Input
            1) A        - coefficients [3 x n], 1=k-1; 2=k; 3=k+1
            2) b        - source term
            3) x_g      - array of initial guesses
            4) tol      - tolerance level for residual
            5) Method   - (1) Jacobi, (2) Gauss-Siedel, (3) SOR
            7) ω        - weights for successive overrelaxation
            8) printint - interval (iterations) to print text
        """
    function FD_itrsolve(A,b,x_init,tol,Method,ω,printint)
        allstart = time()

        # Preallocate Arrays [1 x # of cells]
        ncells = length(x_init)
        x        = zeros(Float64,2,ncells) # Array for current profile (itr=m)

        # Input Initial Guess(2nd row of Tz)
        x[2,:] = x_init

        # Preallocate for residuals
        r   = Float64[]#zeros(Float64,max_iter+1)

        # Set up counters and indexing
        itcnt = 0     # Count of Iterations
        m     = 0     # Iteration Storage Number
        err   = tol+1 # initial value, pre error calculation (need to set up for the diff between guess and actual)

        while  err > tol #m <= max_iter #|| err > tol
            start = time()
            #@printf("Starting iter %i...",m)
            #r_m  = zeros(Float64,1,kmax+1)
            for k = 1:ncells

                # First Cell
                if k == 1
                    x0 = 0
                    x2 = x[2,k+1] # Take T2 from last iteration (or init guess)

                # Last Cell
            elseif k == ncells
                    # For Jacobi Method, use previous values (Row 2 of Tz)
                    if Method == 1
                        x0 = x[2,k-1]

                    # For GS and SOR, use new value (Row 1 of Tz)
                    else
                        x0 = x[1,k-1]
                    end
                    x2 = 0
                # Interior Cells
                else
                    # For Jacobi Method, use previous values (row 2 Tz)
                    if Method == 1
                        x0 = x[2,k-1]
                    # For GS and SOR, use new value (row 1 Tz)
                    else
                        x0 = x[1,k-1]
                    end
                    # Use old T2 value for all cases
                    x2 = x[2,k+1]
                end

                # Get Coefficients and Source Term
                A0 = A[1,k]
                A1 = A[2,k]
                A2 = A[3,k]
                b1 = b[k]

                # Compute Temperature
                x1 = 1/A1 * (b1 - A2*x2 - A0*x0)

                # Apply weights for sucessive overrelaxation
                if Method == 3 && itcnt > 0
                    global x[1,k] = ω * x1 + (1-ω) * x[2,k]
                else
                    global x[1,k] = x1
                end

                if itcnt%printint == 0 && k in [1,ncells,ceil(ncells/2)]
                    #print_itrsolv(itcnt,A0,A1,A2,b1,x0,x1,x2,k,ncells)
                end
            end

            itcnt += 1
            m += 1

            # Make sparse tridiagonal matrix for residual calculation
            A_tri = makeTridiag(A)

            # Calculate residual (can do this more efficiently iteratively)
            err = norm(A_tri*x[1,:] - b)

            # Push calculated error to residual variable
            push!(r,err)

            # Save new Tz profile to Tz0
            x[2,:] = x[1,:] # Note that index 2 is the prev

            elapsed = time() - allstart
            if mod(m,printint)==0
                #@printf("Now on iteration %i with resid %.2E in %fs\n",m,err,elapsed)
            end

        end

        elapsed = time()-allstart
        @printf("Completed in %i iterations with resid %.2E. (%fs)\n",itcnt,err,elapsed)
        return x, itcnt,r
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

        # Set the k+1 cell (same for all cases)
        C[3,1] = κ[k]   / (z_f[k] * z_c[k])

        # Dirichlet Bottom
        if BC_bot == 1
            # Compute Source Term with BC (Prescribed Temp)
            B[1,:]   = val_bot .* (2 * κ0) / (z_f[k] * z_c0)

            # Ck0 remains the same. Multiplier goes on bottom term
            C[1,1] = 0
            C[2,1] = (( 2 * κ0   / z_c0    ) +
                      (     κ[k] / z_c[k]  )) * -1/z_f[k] # Prescribed Temp

        # Newmann Bottom
        elseif BC_bot == 2
            # Compute Source Term with BC (Prescribed Flux)
            B[1,:]   = val_bot ./ z_f[k]

            # Ck2 remains the same, Ck1 dep. on BCs Need negative here!
            C[1,1] = 0
            C[2,1] = - κ[k] / (z_f[k] * z_c[k])

            #Periodic Bottom
        elseif BC_bot == 3
            C[1,1] = κ[kmax] / (z_f[k] * z_c[kmax])
            C[2,1] = (C[1,k] + C[3,k]) * -1
        end

        # Set C(k,k-1) to 0, and set C(k,k+1) to the usual


        # Options for top boundary ----------------------------------------
        k = kmax

        # Set k-1 cells (same for all conditions)
        C[1,kmax] = κ[k-1]     / (z_f[k] * z_c[k-1])

        # Dirichlet Top
        if BC_top == 1

            # Compute Change to Source Term with BC (Prescribed Temp)
            B[2,:]   = val_top .* ( 2 * κ[k]) / (z_f[k] * z_c[k])

            # Calculate C(k,k)
            C[2,kmax] = (( 2 * κ[k]   / z_c[k]  ) +
                         (     κ[k-1] / z_c[k-1])) * -1/z_f[k] # Prescribed Temp
            C[3,kmax] = 0

        # Neumann Top
        elseif BC_top == 2

            # Compute Source Term with BC (Prescribed Flux)
            B[2,:]   = val_top ./ z_f[k]

            # Calculate C(k,k) (need negative here!)
            C[2,kmax] = -κ[k-1] / (z_f[k] * z_c[k-1])# This depends on the type of BC
            C[3,kmax] = 0

        #Periodic Top
        elseif BC_top == 3


            C[3,kmax] = κ[1]       / (z_f[k] * z_c[1])
            C[2,kmax] = (C[1,k] + C[3,k]) * -1

        end

        # Set C(k,k+1) to 0, and set C(k,k-1) to the usual


        # Options for interior ----------------------------------------------------
        for k = 2:kmax-1
            # Compute Coefficients
            C[1,k] = κ[k-1]     / (z_f[k] * z_c[k-1])
            C[3,k] = κ[k]       / (z_f[k] * z_c[k]  )
            C[2,k] = (C[1,k] + C[3,k]) * -1 # Note this might only work in certain cases
        end

        return C,B
    end

    """
    -----------------------------------------------------------
    FD_itrsolve_2D
    -----------------------------------------------------------
        # Iterative Solver for 2D finite-difference problems


        # Inputs:
             Cx       = x coefficients (3 x n)
             Cy       = y coefficients (3 x m)
             S        = Modified Forcing Term (n x m)
             ug       = Initial guess at quantity u (n x m)
             tol      = Error Tolerance
             ω        = Weight for SOR
             method   = [1,Jacobi] [2,Gauss-Siedel] [3,SOR]
             wper     = Periodic Western Boundary (1 = periodic, 0 = not)
             eper     = Periodic Eastern Boundary
             nper     = Periodic Northern Boundary
             maxiter  = Maximum amount of iterations permitted
             saveiter = Amount of iterations to save

         Out:
             u_out   = Final approximation of quantity[n x m]
             itcnt   = # of iterations to convergence
             r       = Array of residuals per iteration
             err_map = Map of the error for the final timestep
    """
    function FD_itrsolve_2D(Cx,Cy,S,ug,tol,ω,method,wper,eper,sper,nper,maxiter,saveiter)
        xmax = size(Cx,2)
        ymax = size(Cy,2)

        # Assign per to individual matrix
        chk_per = zeros(Int,4)
        chk_per[1] = nper
        chk_per[2] = sper
        chk_per[3] = eper
        chk_per[4] = wper

        # Preallocate
        u = zeros(Float64,2,xmax,ymax)
        r = Float64[]

        # Scrap (Delete Later: Save first 10 iterations)
        u_scrap   = zeros(Float64,saveiter,xmax,ymax)
        err_scrap = zeros(Float64,saveiter,xmax,ymax)
        err_map = zeros(Float64,xmax, ymax)

        # Assign ug to the first entry
        u[1,:,:] = ug # 1 will store the previous guess
        u[2,:,:] = ug # 2 will store the updated guess

        itcnt = 0
        start = time()
        while itcnt == 0 || r[itcnt] > tol
            # Loop for each column, by row
            for j = 1:ymax

                # Get Coefficients (y)
                B1  = Cy[1,j]
                #B3y = Cy[2,j]
                B5  = Cy[3,j]

                for i = 1:xmax

                    # Get Coefficients (x)
                    B2 = Cx[1,i]
                    B3 = Cy[2,j] + Cx[2,i]
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
                    if i == 1
                        i2 = xmax
                    else
                        i2 = i-1
                    end
                    if i == xmax
                        i4 = 1
                    else
                        i4 = i+1
                    end

                    # Make j indices periodic
                    if j == ymax
                        j5 = 1
                    else
                        j5 = j+1
                    end
                    if j == 1
                        j1 = ymax
                    else
                        j1 = j-1
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

                    # Now calculate the central point [i,j]
                    u3 = (f - B1*u1 - B2*u2 - B4*u4 - B5*u5) * (1/B3)

                    # Apply weights for SOR
                    if method == 3
                        u[2,i,j] = ω * u3 + (1-ω) * u[1,i,j]
                    else
                        u[2,i,j] = u3
                    end

                end
            # End loop for the point, i,j
            end

            # Multiply A (mat. of coeffs) by x (guess)
            Sg = Ax_2D(Cx,Cy,u[2,:,:],chk_per)

            # Calculate residual
            err = S - Sg
            #err = calc_res_2d(Cx,Cy,S,u[2,:,:],chk_per)

            # Assign this iteration to the last
            u[1,:,:] = u[2,:,:]

            push!(r,norm(err))

            itcnt += 1

            # Save error map and break on last iteration
            if itcnt > maxiter
                err_map = err;
                break
            end

            # Currently set to print on every 10,000th iteration
            if itcnt%10^5 == 0
                elapsed = time()-start
                @printf("\nOn Iteration %i in %s s",itcnt,elapsed)
            end

            # Scrap: Save first "saveiter" iterations and error for animation
            if itcnt <=saveiter
                u_scrap[itcnt,:,:]=u[1,:,:]
                err_scrap[itcnt,:,:]=err
            end

        end
        u_out = u[1,:,:]

        elapsed = time()-start
        @printf("\nFinished %.2e iterations in %s",itcnt,elapsed)
        return u_out, itcnt, r, u_scrap, err_scrap, err_map

    end

    """
    -----------------------------------------------------------
    2D_calc_residual(Cx,Cy,S,x,chk_per)
    -----------------------------------------------------------
        Calculate the residual (r = b - Ax), where:
        x       = some initial guess [i x j]
        b       = Source term (S [i x j]_
        A       = Matrix of coefficients in x and y directions (Cx [5 x i], Cy [5 x j])
        chk_per = boolean (N-per,S-per,E-per,W-per)

        Inputs:
            1) Cx      = Coefficients in x-direction [5 x i]
            2) Cy      = Coefficients in y-direction [5 x j]
            3) S       = Source term                 [i x j]
            4) x       = Guess                       [i x j]
            5) chk_per = Boolean for periodic BC     [4 (N,S,E,W)], 1 = period, 0 = not

        Output:
            1) res     = residual r = b - Ax [i x j]
    """
    function calc_res_2d(Cx,Cy,S,x,chk_per)
        xmax = size(Cx,2)
        ymax = size(Cy,2)

        nper = chk_per[1]
        sper = chk_per[2]
        eper = chk_per[3]
        wper = chk_per[4]

        res = zeros(Float64,xmax,ymax)

        for j = 1:ymax

            # Get Coefficients (y)
            B1  = Cy[1,j]
            B5  = Cy[3,j]

            for i = 1:xmax
                # Get Coefficients (x)
                B2 = Cx[1,i]
                B3 = Cy[2,j] + Cx[2,i]
                B4 = Cx[3,i]

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
                if j == 1
                    j1 = ymax
                end
                if j == ymax
                    j5 = 1
                end

                # Interior Points, Periodic
                u1 = x[i,j1]
                u2 = x[i2,j]
                u3 = x[i,j]
                u4 = x[i4,j]
                u5 = x[i,j5]

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

                res[i,j] = f - (B1*u1+B2*u2+B3*u3+B4*u4+B5*u5)
            end
        end
        return res
    end

    """
    -----------------------------------------------------------
    Ax_2D(Cx,Cy,x,chk_per)
    -----------------------------------------------------------
        # Function to compute Ax = b iteratively where
            x       = some initial guess [i x j]
            A       = Matrix of coefficients in x and y directions (Cx [5 x i], Cy [5 x j])
            chk_per = boolean (N-per,S-per,E-per,W-per)

        Inputs:
            1) Cx      = Coefficients in x-direction [5 x i]
            2) Cy      = Coefficients in y-direction [5 x j]
            4) x       = Guess                       [i x j]
            5) chk_per = Boolean for periodic BC     [4 (N,S,E,W)], 1 = period, 0 = not
    """
    function Ax_2D(Cx,Cy,x,chk_per)
        xmax = size(Cx,2)
        ymax = size(Cy,2)

        nper = chk_per[1]
        sper = chk_per[2]
        eper = chk_per[3]
        wper = chk_per[4]

        b = zeros(Float64,xmax,ymax)

        for j = 1:ymax

            # Get Coefficients (y)
            B1  = Cy[1,j]
            B5  = Cy[3,j]

            for i = 1:xmax
                # Get Coefficients (x)
                B2 = Cx[1,i]
                B3 = Cy[2,j] + Cx[2,i]
                B4 = Cx[3,i]

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
                if j == 1
                    j1 = ymax
                end
                if j == ymax
                    j5 = 1
                end

                # Interior Points, Periodic
                u1 = x[i,j1]
                u2 = x[i2,j]
                u3 = x[i,j]
                u4 = x[i4,j]
                u5 = x[i,j5]

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

                b[i,j] = B1*u1 + B2*u2 + B3*u3 + B4*u4 + B5*u5
            end
        end
        return b
    end

    """
    -----------------------------------------------------------
    cgk_2d(Cx,Cy,S,xg,chk_per,tol,maxiter)
    -----------------------------------------------------------

        Conjugate-Gradient Krylov Method, applied for solving the 2D
        problem Ax = b, where A contains the coefficients for x and y,
        x is the vector of the target quantities and b is the source term S

        Inputs:
            1) Cx      - Coefficients in the x-direction
            2) Cy      - Coefficients in the y-direction
            3) S       - Source term
            4) xg      - Guess at the quantity x
            5) chk_per - Boolean for periodic BCs
            6) tol     - Tolerance for residual
            7) maxiter - maximum number of iterations

        Outputs
            1) x_out   - final guess for x
            2) itcnt   - count of iterations
            3) res     - residual for each iteration
    """
    function cgk_2d(Cx,Cy,S,xg,chk_per,tol,maxiter)
        start = time()
        xmax = size(Cx,2)
        ymax = size(Cy,2)
        # Preallocate array of residuals
        res   = Float64[]           # Array to store stepsize
        r0    = zeros(xmax,ymax)
        d0    = zeros(xmax,ymax)    # Array to store directions [1=prev, 2=old]
        x_out = zeros(xmax,ymax)    # Array of the quantity x [1=prev, 2=old]
        itcnt = 0 # Iteration count
        ridx  = 1 # Count for indexing residuals
        push!(res,tol+1)

        # Compute first residual and direction
        Ax  = Ax_2D(Cx,Cy,xg,chk_per)
        r0  = S - Ax
        d0  = r0
        x0  = xg
        if norm(r0) < tol
            x_out = x0
            push!(res,norm(r0))
            @printf("\nGuess is already below tolerance with residual %.2e",norm(r0))
            return x_out,itcnt,res
        end

        while itcnt < maxiter

            # # ------------------------------------------
            # # 1. Compute Step Size α = d^2 / (d * A * d)
            # # ------------------------------------------
            Ad0 = Ax_2D(Cx,Cy,d0,chk_per)
            α   = sum(r0.^2) / sum(d0 .* Ad0)

            # # ------------------------------------------
            # # 2. Get new x and new residual
            # # ------------------------------------------

            x   = x0 + α*d0

            Ax1 = Ax_2D(Cx,Cy,x,chk_per)
            r1  = S - Ax1
            push!(res,norm(r1))
            if norm(r1) < tol
                x_out = x
                break
            end

            # # ------------------------------------------
            # # 3. Calculate new direction
            # # ----


            β = sum(r1.^2) / sum(r0.^2)
            d1 = r1 + β*d0


            itcnt +=1
            ridx  +=1
            r0 = r1
            d0 = d1
            x0 = x

        end
        elapsed = time()-start
        @printf("\nFinished %.2e iterations in %s",itcnt,elapsed)
        return x_out,itcnt,res
    end

    """
    -----------------------------------------------------------
    ddx_1d(ϕ,Δx,method)
    -----------------------------------------------------------

        Takes derivative of property ϕ where
        ϕ is given at the edges of the cell and
        Δx is the cell width. For each cell using one of the methods:

        Method (1) - Forward Differencing

            ϕ_x = [ϕ(i+1)-ϕ(i)] / Δx(i)

        Method (2) - Backward Differencing

            ϕ_x = [ϕ(i)-ϕ(i-1)] / Δx(i)

        Method (3) - Central Differencing

            ϕ_x = [ϕ(i+1)-ϕ(i-1)] / 2Δx(i)

        Inputs
          1) ϕ      = matrix of size [i]
          2) Δx     = cell spacing of size [1,imax]
          3) method = (1,2,or 3), see above
        Outputs
          1) ϕ_x    = Derivative of ϕ
          2) x_new  = New x indices
          3) Δx_new = Distances of new cell
    """
    function ddx_1d(ϕ,Δx,method)

        # Determine final size based on method
        imax = length(Δx)

        # For forward differencing, loop(1:n-1)
        if method == 1
            istart = 1
            imax_new = imax-1
            isize  = imax_new
        # For backward differencing, loop(2:n)
        elseif method == 2
            istart = 2
            imax_new = imax
            isize = imax-1
        # For central differencing, loop(2:n-1)
        elseif method == 3
            imax_new = imax-1
            isize    = imax-2
            istart   = 2
        end


        # Preallocate
        ϕ_x  = zeros(Float64,isize)
        x_new = [istart:imax_new;]
        Δx_new = Δx[istart:imax_new]

        icnt = 1
        for i = istart:imax_new

            # Forward Differencing
            if method == 1
                ϕ_x[icnt] = (ϕ[i+1] - ϕ[i]) / Δx[i]

            # Backward Differencing
            elseif method == 2
                ϕ_x[icnt] = (ϕ[i] - ϕ[i-1]) / Δx[i]

            # Central Differencing
            elseif method == 3
                ϕ_x[icnt] = (ϕ[i+1] - ϕ[i-1]) / (2*Δx[i])
            end
            icnt +=1
        end
        return ϕ_x,x_new,Δx_new
    end

    """
    -----------------------------------------------------------
    CN_make_matrix_2d
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
    function CN_make_matrix_2d(Δt,θ,IC,Cx,Cy,Bx,By,fr,fl,kprint,chk_per)

        start = time()

        # Determine LHS and RHS multipliers
        # LHS - For timestep (n+1), multiply by θ
        l_mult =  -Δt*(θ)
        # RHS - For timestep (n)  , multiply by 1-θ
        r_mult =  Δt*(1-θ)


        # Multiply variables by time and theta factors
        Bx      = Bx  .* l_mult
        By      = By  .* l_mult
        fl      = fl .* l_mult
        Cx      = Cx  .* r_mult
        Cy      = Cy  .* r_mult
        fr      = fr .* r_mult

        # Meth2: Add single digit post-multiplication
        # To diagonal (Ck,Bk)
        Cx[2,:] = Cx[2,:] .+ 1
        #Cy[2,:] = Cy[2,:] .+ 1
        By[2,:] = By[2,:] .+ 1
        #Bx[2,:] = Bx[2,:] .+ 1

        # Now combine terms
        Cx = Ax_2D(Cx,Cy,IC,chk_per)
        b     = Cx + fr - fl # Bring Sl from left side

        Ax     = Bx
        Ay     = By

        k = kprint
        #print_exeq(k,length(fr),b,C,IC,fr,fl,A)

        elapsed = time() - start
        #@printf("Completed matrix prep in %fs\n",elapsed)
        return Ax,Ay,b
    end

    """
    -----------------------------------------------------------
    PV Inversion Script (invPV_2d)
    -----------------------------------------------------------
        Aggregation of HW3 into a script.
        Given a vorticity field, computes the streamfunction

        Dependencies (from ocnmod module):
            FD_calc_coeff_2D
            cgk_2d

        Inputs

            1)  ζ    - [ixj] vorticity field

            2)  x_f  - [1xi] cell widths in x-dim
            3)  y_f  - [1xj] cell widths in y-dim
            4)  x_c  - [1xi] midpoint distance in x-dim
            5)  y_c  - [1xj] midpoint distance in y-dim
            6)  x_c0 - [1]   first midpoint distance in x-dim
            7)  y_c0 - [1]   first midpoint distance in y-dim
            8)  κx   - [1xi] diffusivity in x-dim
            9)  κy   - [1xj] diffusivity in y-dim
            10) κx0  - [1]   first diffusivity in x-dim
            11) κy0  - [1]   first diffusivity in y-dim

            12) NBC     - Type of BC for North (j = jmax)
            13) nb_val  - [ix1] values along nb
            14) SBC     - Type of BC for South (j = 1)
            15) sb_val  - [ix1] values along sb
            16) EBC     - Type of BC for East  (i = imax)
            17) eb_val  - [jx1] values along eb
            18) WBC     - Type of BC for West  (i = 1)
            19) wb_val  - [jx1] values along wb

            20) ug       - [ixj] Array of initial guess for the value
            21) tol      - Minimum tolerance for the residual
            22) max_iter - Maximum number of interations permitted

        Outputs
            1) ψ     - Streamfunction
            2) itcnt - iterations to convergence
            3) r     - residual at each iteration
    """
    function invPV_2d(ζ,
                      x_f,y_f,x_c,y_c,x_c0,y_c0,
                      κx,κy,κx0,κy0,
                      NBC,nb_val,SBC,sb_val,EBC,eb_val,WBC,wb_val,
                      ug,tol,max_iter)

        # Get # of cells along each dimension
        xmax = size(x_f,2)
        ymax = size(y_f,2)

        # Create chk_per vector based on BCs
        # [N, S, E, W]
        chk_per = zeros(Int,4)
        if NBC == 3
            chk_per[1] = 1
        end
        if SBC == 3
            chk_per[2] = 1
        end
        if EBC == 3
            chk_per[3] = 1
        end
        if WBC == 3
            chk_per[4] = 1
        end

        # Calculate vectors to store coefficients for iteration
        Cx,Bx = ocnmod.FD_calc_coeff_2D(xmax,x_f,x_c,κx,EBC,eb_val,WBC,wb_val,x_c0,κx0)
        Cy,By = ocnmod.FD_calc_coeff_2D(ymax,y_f,y_c,κy,NBC,nb_val,SBC,sb_val,y_c0,κy0)

        # Modify source term with BCs
        ζ[1,:]    -= Bx[1,:] # South BC
        ζ[xmax,:] -= Bx[2,:] # North BC
        ζ[:,1]    -= By[1,:] # West BC
        ζ[:,ymax] -= By[2,:] # East BC

        # Use Conjugate Gradient Krylov to solve for the streamfunction
        ψ,itcnt,r = ocnmod.cgk_2d(Cx,Cy,ζ,ug,chk_per,tol,max_iter)

        return ψ,itcnt,r
    end

    """
    -----------------------------------------------------------
    Quiverplot 2D
    -----------------------------------------------------------
        Setup for quiverplot

        Spaces points by x_sp, y_sp, and scales variables

        Inputs
            1) x        = vector: x points
            2) y        = vector: y points
            3) u        = array [i x j]
            4) v        = array [i x j]
            5) x_sp     = x-spacing (plot every x_sp points)
            6) y_sp     = y-spacing (plot every y_sp points)
            5) qscale   = number value to scale up vectors

        Outputs
            1) pts      = vector [i*j] of tuples, coordinate pairs
            2) uv       = vector [i*j] of tuples, vector component values
    """
    function quiverprep_2d(x,y,u,v,x_sp,y_sp,qscale)
        # Space variables
        us  = u[1:x_sp:end-1,1:y_sp:end-1] .* qscale
        vs  = v[1:x_sp:end-1,1:y_sp:end-1] .* qscale
        xp  = x[1:x_sp:end-1]; xpm = length(xp)
        yp  = y[1:y_sp:end-1]; ypm = length(yp)

        # Combine
        pts = vec([(xp[i], yp[j]) for i=1:xpm, j=1:ypm])
        uv = vec([(us[i,j],vs[i,j]) for i=1:xpm, j=1:ypm])
        return pts, uv
    end

    """
    -----------------------------------------------------------
    Streamfunction to uv: psi2uv(ψ,Δx,Δy,x,y,dropdim)
    -----------------------------------------------------------
    Using the relations:
        dψ/dx = v
        dψ/dy = -u

    Recover the velocity and wind field from the input ψ.
    Use central differencing.

    Dependencies: ddx_1d

    Inputs
        1) ψ       - streamfuncion [i x j] array
        2) Δx      - x cell size   [i]
        3) Δy      - y cell size   [j]
        4) x       - x coordinates [i]
        5) y       - y coordinates [j]
        6) dropdim - boolean, 1= drop outside dimensions for central differencing


    """
    function psi2uv(ψ,Δx,Δy,x,y,dropdim)
        xmax = length(x)
        ymax = length(y)

        # Calculate dψ/dx along each row
        v = zeros(xmax-2,ymax)
        for j = 1:ymax
            ψj = ψ[:,j]
            v[:,j],~,~ = ocnmod.ddx_1d(ψj,Δx,3)
        end

        # Calculate dψ/dy along each column
        u = zeros(xmax,ymax-2)
        for i = 1:xmax
            ψi = ψ[i,:]
            u[i,:],~,~ = ddx_1d(ψi,Δy,3)
        end
        u *= -1

        # Drop extra dimensions for the match
        if dropdim == 1
            u = u[2:end-1,:]
            v = v[:,2:end-1]
        end
        # Do the same to the input grids
        nx = x[2:end-1]
        ny = y[2:end-1]

        return u,v,nx,ny
    end

    """
    -----------------------------------------------------------
    ddx_2d()
    -----------------------------------------------------------
    Using the relations:
        dψ/dx = v
        dψ/dy = -u

    ddx_1d but along a 2D matrix

    Dependencies: ddx_1d

    Inputs
        1) ψ       - streamfuncion [i x j] array
        2) Δx      - x cell size   [i]
        4) x       - x coordinates [i]
        5) y       - y coordinates [j]
        6) dim     - 1 = calculate along row (j), 2 = calculate along column (i)


    """
    function ddx_2d(ψ,Δx,x,y,dim)
        xmax = length(x)
        ymax = length(y)

        # Calculate dψ/dx, looping for each j
        if dim == 1
            u = zeros(xmax-2,ymax)
            for j = 1:ymax
                ψj = ψ[:,j]'
                u[:,j],b,c = ocnmod.ddx_1d(ψj,Δx,3)
            end

            u = u[:,2:end-1]
            nx = x[2:end-1]
            ny = y
        # Calculate dψ/dy, looping along 2nd dim (row)
        elseif dim == 2
            u = zeros(xmax,ymax-2)
            for i = 1:xmax
                ψi = ψ[i,:]'
                u[i,:],~,~ = ddx_1d(ψi,Δx,3)
            end
            u *= -1
            u = u[2:end-1,:]
            nx = x
            ny = y[2:end-1]
        end
        return u,nx,ny
    end

















# Module End
end
