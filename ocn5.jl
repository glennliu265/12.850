module ocn5
    # Test Module for assignment 5, 1D Advection
    using Plots
    using Printf
    using LinearAlgebra


    """
    -----------------------------------------------------------
    UW_calc_coeff_1D
    -----------------------------------------------------------
        Compuee both C (Matrix of Coefficients) and modification to
        corresponding forcing term

        NOTE: Assumes Forcing term is on the same side of the equation
        as the diffusion term.

        Inpues:
            ||~ Box Geometry and Indices ~||
            z_f    = cell spacing
            zmax   = size of cells in target direction

            ||~ Initial Profiles         ~||
            u      = velocity field

            ||~ Boundary Conditions      ~||
            BC_top  = 1 for Dirichlet, 2 for Neumann
            val_top = Inpue value for top (matches size of other dim)
            BC_bot  = 1 for Dirichlet, 2 for Neumann
            val_bot = Inpue value for bot (matches size of other dim)

        Ouepues:
            1) C - Matrix of Coefficients [1:3,kmax]
                1 = C(k,k-1); 2 = C(k,k); 3 = C(k,k+1)
            2) B - Column vector of modification to Source/Sink Terms, with
                first and last entries modified for BCs along the other dimension
        """
    function UW_calc_coeff_1D(x_f,xmax,          # Cell face width [i]
        u,                                  # velocity field [i+1]
        BC_top,val_top,                     # Bottom BCs [1]
        BC_bot,val_bot)                     # Top BCs [1]


        C = zeros(Float64,3,xmax)
        B = zeros(Float64,xmax) # 1st dim, 1=bot,2=top,

        for i = 1:xmax
            # # For the bottom boundary
            # if i == 1
            #
            #     ip1 = i+1
            #     # For bottom, index from the top
            #     im1 = xmax
            #
            #     # Top boundary
            # elseif i == xmax
            #
            #         # Index top points back to bottom
            #         ip1 = 1
            #         im1 = i-1
            # else
            #
            #     ip1 = i+1
            #     im1 = i-1
            #
            # end

            # -------------
            # Get the winds
            uw = u[i]
            ue = u[i+1]

            # Calculate u+ and u- for top and bottom
            uep = (ue + abs(ue))/2 # u-top plus (ue)
            uem = (ue - abs(ue))/2 # u-bot minus (ue)

            uwp = (uw + abs(uw))/2 # u-bottom plus (uw)
            uwm = (uw - abs(uw))/2 # u-bottom minus (uw)

            dx  = -1*x_f[i]

            # ----------------------
            # Calculate Coefficients


            # Begin by indexing, assuming as if everything was periodic
            C[1,i] = -uwp#uwp / dx
            C[2,i]   = uep-uwm#-1* (uep - uwm) / dx
            C[3,i] = uem #uem / dx

            # Bottom Boundaries (assuming both im1 and im2 draw from same pool)
            if i == 1
                # Dirichlet BC
                if  BC_bot == 1

                    # Set i-1 cell to zero and move to forcing
                    C[1,i] = 0
                    B[i]     = 2*uwp*val_bot#2 * uwp * val_bot / dx

                    # Alter Interior Cell
                    C[2,i]   = 2*uep-uwm

                # Neumann BC
                elseif BC_bot == 2

                    # Set i-1 cell to zero and move to forcing
                    C[1,i] = 0
                    B[i]     = -val_bot#val_bot / dx

                    # Alter Interior Cell
                    C[2,i]   = uep#-1*uep / dx
                end
            end

            # Top Boundary
            if i == xmax

                # Set the i-1 cell (same as interior)
                #C[1,im1,j] = uwp / dz

                if BC_top == 1

                    C[3,i] = 0
                    B[i]     = 2*uem*val_top#-2*uem / dx * val_top

                    C[2,i]   = uep-2*uwm#-1*(uep-2*uwm) / dx

                elseif BC_top == 2

                    C[3,i] = 0
                    B[i]     = val_top*uem#-1*val_top / dx

                    C[2,i]   = uep-uwm#uwm / dx

                end


            end
            C[1,i] /= dx
            C[3,i] /= dx
            C[2,i] /= dx
            B[i] /= dx

        end

        return C,B
    end


    """
    -----------------------------------------------------------
    CN_make_matrix_1D
    -----------------------------------------------------------
    """

    function CN_make_matrix_1D(Δt,θ,IC,C,B,fr,fl,x_f,chk_per)

        start = time()

        # Determine LHS and RHS multipliers
        # LHS - For timestep (n+1), multiply by θ
        l_mult =  -Δt*(θ)
        # RHS - For timestep (n)  , multiply by 1-θ
        r_mult =  Δt*(1-θ)


        # Multiply variables by time and theta factors
        B      = B  .* l_mult
        fl     = fl .* l_mult
        C      = C  .* r_mult
        fr     = fr .* r_mult

        # Meth2: Add single digit post-multiplication
        # To diagonal (Ck,Bk)
        C[2,:] = C[2,:] .+ 1
        B[2,:] = B[2,:] .+ 1

        # Now combine terms
        Cx = Ax_1D(C,IC,chk_per)

        b     = Cx + fr - fl # Bring Sl from left side

        Ax     = B

        elapsed = time() - start
        #@printf("Completed matrix prep in %fs\n",elapsed)
        return Ax,b
    end

    """
    -----------------------------------------------------------
    Ax_1D(A,x,chk_per)
    -----------------------------------------------------------
        # Function to compuee Ax = b iteratively where
            x       = some initial guess [i x j]
            A       = Matrix of coefficients in x and y directions (Cx [5 x i], Cy [5 x j])
            chk_per = boolean (N-per,S-per,E-per,W-per)

        Inpues:
            1) Cx      = Coefficients in x-direction [5 x i]
            2) Cy      = Coefficients in y-direction [5 x j]
            4) x       = Guess                       [i x j]
            5) chk_per = Boolean for periodic BC     [4 (N,S,E,W)], 1 = period, 0 = not
    """
    function Ax_1D(A,x,chk_per)
        xmax = length(x)

        nper = chk_per[1]
        sper = chk_per[2]

        b = zeros(Float64,xmax)

        for i = 1:xmax

            B1 = A[1,i]
            B2 = A[2,i]
            B3 = A[3,i]


            # First, assume periodic boudaries
            # Make i indices periodic
            im1 = i-1
            ip1 = i+1

            if i == 1
                im1 = xmax
            end
            if i == xmax
                ip1 = 1
            end

            # Interior Points, Periodic
            x1 = x[im1]
            x2 = x[i]
            x3 = x[ip1]

            # BCs
            if sper != 1 && i == 1
                x1 = 0 # All j-1 terms = 0
            end

            if nper != 1 && i == ymax
                x3 = 0 # All j+1 terms = 0
            end

            b[i] = B1*x1 + B2*x2 + B3*x3
        end
        return b
    end

end
