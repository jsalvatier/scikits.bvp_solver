!! this file makes up the interface which is wrapped by python

Module BVP
		use BVP_M,
		implicit none

	! define some member parameters that python can use to get the solution parameters
	! this is necessary because f2py cannot pass 2-d arrays back from fortran

	!these are all the solution out parameters and they are
		INTEGER :: NODE,NPAR ! number of ODEs; number of unknown parameters.
   		INTEGER :: LEFTBC ! number of left boundary conditions.
   		INTEGER :: NPTS,INFO ! number of points in current mesh; indicator of
                         ! success (INFO=0) or failure (INFO=-1) of computation.
    	INTEGER :: MXNSUB ! maximum number of subintervals allowed.
  		DOUBLE precision, allocatable, DIMENSION(:), target :: X ! current mesh.

    	DOUBLE PRECISION, allocatable, DIMENSION(:,:), target :: Y ! current solution - Ith
	!   column is solution approximation at Ith mesh point in first NODE locations;
	!   last NPAR locations contain approximation to unknown parameters.

    	DOUBLE PRECISION, allocatable, DIMENSION(:),target  :: PARAMETERS ! Unknown parameters.

    	DOUBLE PRECISION, allocatable, DIMENSION(:) , target :: WORK ! double precision workspace.
    	INTEGER, allocatable, DIMENSION(:) :: IWORK ! integer workspace.

	!   The workspaces are used to save information that can be accessed by SOL_EVAL
	!   for evaluation of the numerical solution.

	! output arguments:

		DOUBLE PRECISION :: YERROR

	! member parameters for returning the results of an evaluation call
	    DOUBLE PRECISION, allocatable, DIMENSION(:,:) :: EVALUATED
	   	DOUBLE PRECISION, allocatable, DIMENSION(:,:) :: EVALUATED_D

PUBLIC :: GUESS_1_WRAP, GUESS_2_WRAP, GUESS_3_WRAP,   BVP_SOLVER_WRAP
!This is a wrapper for BVP_M so that the relevant BVP subroutines do not pass back
!derived data types but normal data types
contains

!!this function takes a set of parameters that would normally make up a SOL data type in fortran and turns it into a SOL data type
!! this is necessary because C does not like fortran derived data types
function sol_from_params  (NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, &
 MXNSUB_in ,X_in,  Y_in, PARAMETERS_in, work_in, iwork_in) result(sol_out)
		INTEGER :: NODE_in,NPAR_in ! number of ODEs; number of unknown parameters.
   		INTEGER :: LEFTBC_in ! number of left boundary conditions.
   		INTEGER :: NPTS_in,INFO_in ! number of points in current mesh; indicator of
                         ! success (INFO=0) or failure (INFO=-1) of computation.
    	INTEGER :: MXNSUB_in ! maximum number of subintervals allowed.

  		DOUBLE PRECISION, DIMENSION(:), target :: X_in ! current mesh.

    	DOUBLE PRECISION, DIMENSION(:,:), target :: Y_in ! current solution - Ith
	!   column is solution approximation at Ith mesh point in first NODE locations;
	!   last NPAR locations contain approximation to unknown parameters.

    	DOUBLE PRECISION, DIMENSION(:), target :: PARAMETERS_in ! Unknown parameters.

    	INTEGER, DIMENSION(:), target :: IWORK_in ! integer workspace.
    	DOUBLE PRECISION, DIMENSION(:), target :: WORK_in ! double precision workspace.
	!   The workspaces are used to save information that can be accessed by SOL_EVAL
	!   for evaluation of the numerical solution.

		TYPE(BVP_SOL) :: sol_out

	! copy all the passed parameters to the sol_out structure


		sol_out%NODE = NODE_in
    	sol_out%NPAR = NPAR_in
    	sol_out%LEFTBC = LEFTBC_in
		sol_out%MXNSUB = MXNSUB_in
		sol_out%NPTS = NPTS_in
		sol_out%INFO = INFO

		sol_out%X => X_in
		sol_out%Y => Y_in

		IF (NPAR_in > 0) THEN
		  sol_out%PARAMETERS => PARAMETERS_in
		ELSE
		  NULLIFY(sol_out%PARAMETERS)
		END IF

		if (size(WORK_in) > 0) then

			sol_out%WORK => WORK_in
		else
		  NULLIFY(sol_out%WORK)
		 end if


		if (size(iWORK_in) > 0) then

			sol_out%iWORK => iWORK_in
		else
		  NULLIFY(sol_out%iWORK)
		 end if

end function


subroutine mparams_f_sol( sol_in)
		TYPE(BVP_SOL), intent(in) :: sol_in


	! copy all  the sol_out structure parameters to the output parameters
		NODE = sol_in%NODE
		NPAR = sol_in%NPAR
		NPTS = sol_in%NPTS
		INFO = sol_in%INFO
		LEFTBC = sol_in%LEFTBC
    	MXNSUB = sol_in%MXNSUB

    	IF( ALLOCATED(X) ) DEALLOCATE( X )
    	IF( ALLOCATED(Y) ) DEALLOCATE( Y )
    	IF( ALLOCATED(PARAMETERS) ) DEALLOCATE( PARAMETERS )
    	IF( ALLOCATED(WORK) ) DEALLOCATE( WORK)
    	IF( ALLOCATED(IWORK) ) DEALLOCATE( IWORK)


    	!make sure there's something here to copy
    	!then allocate and copy the appropriate arrays


    	IF (ASSOCIATED(sol_in%X)) ALLOCATE(X(size(sol_in%X)))
    	IF (ASSOCIATED(sol_in%Y)) ALLOCATE(Y(size(sol_in%Y,1), size(sol_in%Y,2))  )
    	IF (ASSOCIATED(sol_in%PARAMETERS)) ALLOCATE(PARAMETERS(size(sol_in%PARAMETERS)))
    	IF (ASSOCIATED(sol_in%WORK)) ALLOCATE(WORK(size(sol_in%WORK)))
    	IF (ASSOCIATED(sol_in%IWORK)) ALLOCATE(IWORK(size(sol_in%IWORK)))

		IF (ASSOCIATED(sol_in%X)) X = sol_in%X
		IF (ASSOCIATED(sol_in%Y)) Y = sol_in%Y
		IF (ASSOCIATED(sol_in%PARAMETERS)) PARAMETERS = sol_in%PARAMETERS
		IF (ASSOCIATED(sol_in%WORK)) WORK = sol_in%WORK
		IF (ASSOCIATED(sol_in%IWORK)) IWORK = sol_in%IWORK
end subroutine



!! this wraps part of the BVP_INIT interface in fortran (I don't know how to use interface like structures in C and I would guess you cant, so this is necessary)
!! this part is for when Y guess is given as just constant values
subroutine GUESS_1_WRAP(NODE_in,NParam_in, LEFTBC_in, n, X_in, Y_in, np, &
PARAMETERS_in, MXNSUB_in)
	!   Input arguments:
    	INTEGER :: NODE_in, NParam_in, LEFTBC_in, n , np
    	DOUBLE PRECISION, DIMENSION(n) :: X_in
!f2py 	integer intent (hide), depend (X_in) :: n = len(X_in)
    	double precision, dimension (NODE_in) :: Y_in
    	DOUBLE PRECISION, DIMENSION(np) :: PARAMETERS_in
!f2py 	integer intent (hide), depend (PARAMETERS_in) :: np = len(PARAMETERS_in)
    	INTEGER :: MXNSUB_in

	    ! intermediate variable
		TYPE(BVP_SOL) :: init_sol

		init_sol = BVP_INIT(NODE_in,NParam_in, LEFTBC_in,X_in,Y_in,PARAMETERS_in, MXNSUB_in)

		! this gets the variables from the solution structure so we can return them
		call mparams_f_sol(init_sol)

		! this is necessary because all of the parameters in init_sol are copied


end subroutine


!! this wraps part of the BVP_INIT interface in fortran (I don't know how to use interface like structures in C and I would guess you cant, so this is necessary)
!! this part is for when Y guess is given as a list of points
subroutine GUESS_2_WRAP(NODE_in,NParam_in, LEFTBC_in,n, X_in,Y_in,&
 np, PARAMETERS_in, MXNSUB_in)
	!   Input arguments:
    	INTEGER :: NODE_in,NParam_in, LEFTBC_in , n, np
    	DOUBLE PRECISION, DIMENSION(n) :: X_in
!f2py 	integer intent (hide), depend (X_in) :: n = len(X_in)

    	DOUBLE PRECISION, DIMENSION(NODE_in,n) :: Y_in
    	DOUBLE PRECISION, DIMENSION(np), OPTIONAL :: PARAMETERS_in
!f2py 	integer intent (hide), depend (PARAMETERS_in) :: np = len(PARAMETERS_in)
    	INTEGER, OPTIONAL :: MXNSUB_in

	! intermediate variable
		TYPE(BVP_SOL) :: init_sol

		init_sol =  BVP_INIT(NODE_in,NParam_in, LEFTBC_in,X_in,Y_in,PARAMETERS_in, MXNSUB_in)

	! this gets the variables from the solution structure so we can return them
		call mparams_f_sol(init_sol)
end subroutine



!! this wraps part of the BVP_INIT interface in fortran (I don't know how to use interface like structures in C and I would guess you cant, so this is necessary)
!! this part is for when Y guess is given as a function
subroutine GUESS_3_WRAP(NODE_in,NParam_in, LEFTBC_in, n ,&
X_in,FCN, np, PARAMETERS_in, MXNSUB_in)
	!   Input arguments:
    	INTEGER :: NODE_in,NParam_in, LEFTBC_in , n, np
    	DOUBLE PRECISION, DIMENSION(n) :: X_in
!f2py 	integer intent (hide), depend (X_in) :: n = len(X_in)

!f2py 	intent(callback) fcn
    	EXTERNAL FCN
! 		The following lines define the callback signature of fcn
!f2py 	double precision intent(in) :: x
!f2py   double precision dimension(:), intent(out) :: y
!f2py 	call fcn(x ,y)

    	DOUBLE PRECISION, DIMENSION(np), OPTIONAL :: PARAMETERS_in
!f2py 	integer intent (hide), depend (PARAMETERS_in) :: np = len(PARAMETERS_in)
    	INTEGER, OPTIONAL :: MXNSUB_in

	! intermediate variable
	TYPE(BVP_SOL) :: init_sol

	init_sol = BVP_INIT(NODE = NODE_in,NParam = NParam_in,LEFTBC = LEFTBC_in, &
					 X = X_in,FCN = FCN,P = PARAMETERS_in,MAX_NUM_SUBINTERVALS= MXNSUB_in)


	! this gets the variables from the solution structure so we can return them
	call mparams_f_sol(init_sol)

end subroutine

!! this subroutine wraps the actual BVP solver functionality
!! boy does this function have a lot of arguments
subroutine BVP_SOLVER_WRAP(NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, MXNSUB_in &
, n ,X_in , ny1, ny2, Y_in, np,  PARAMETERS_in, work_in , nwork, iwork_in, niwork, &
FSUB,FSUBP,BCSUB,BCSUBP, singular, ns1, ns2, SINGULARTERM,METHOD, &
TOL,hasDFDY, DFDY,DFDYP, hasDBCDY, DBCDY,DBCDYP,TRACE,STOP_ON_FAIL)
	!   Input arguments:

	! replacement for SOL input
    	! output arguments:
		INTEGER :: NODE_in,NPAR_in ! number of ODEs; number of unknown parameters.
   		INTEGER :: LEFTBC_in ! number of left boundary conditions.
   		INTEGER :: NPTS_in,INFO_in ! number of points in current mesh; indicator of
                         ! success (INFO=0) or failure (INFO=-1) of computation.
    	INTEGER :: MXNSUB_in ! maximum number of subintervals allowed.

	    integer :: n, ny1, ny2, np, nwork, niwork

  		DOUBLE PRECISION, DIMENSION(n) :: X_in ! current mesh.
!f2py 	integer intent (hide), depend (X_in) :: n = len(X_in)

    	DOUBLE PRECISION, DIMENSION(ny1,ny2) :: Y_in ! current solution - Ith
!f2py 	integer intent (hide), depend (Y_in) :: ny1 = shape(Y_in, 0)
!f2py 	integer intent (hide), depend (Y_in) :: ny2 = shape(Y_in, 1)
	!   column is solution approximation at Ith mesh point in first NODE locations;
	!   last NPAR locations contain approximation to unknown parameters.

    	DOUBLE PRECISION, DIMENSION(np) :: PARAMETERS_in ! Unknown parameters.
!f2py 	integer intent (hide), depend (PARAMETERS_in) :: np = len(PARAMETERS_in)

    	INTEGER, DIMENSION(niwork) :: IWORK_in ! integer workspace.


    	DOUBLE PRECISION, DIMENSION(nwork) :: WORK_in ! double precision workspace.

	!   The workspaces are used to save information that can be accessed by SOL_EVAL
	!   for evaluation of the numerical solution.



	! input arguments for solving

		logical :: singular
		integer :: ns1, ns2
		DOUBLE PRECISION, DIMENSION(ns1,ns2), OPTIONAL :: SINGULARTERM
!f2py 	integer intent (hide), depend (singularterm) :: ns1 = shape(singularterm, 0)
!f2py 	integer intent (hide), depend (singularterm) :: ns2 = shape(singularterm, 1)
    	INTEGER, OPTIONAL :: METHOD,TRACE
    	LOGICAL, OPTIONAL :: STOP_ON_FAIL
    	DOUBLE PRECISION, OPTIONAL :: TOL

	!   User-supplied subroutines:

		logical :: hasDFDY, hasDBCDY

!	function at Y and T
!f2py 	intent(callback) FSUB
		external FSUB,FSUBP

! 		the following lines define the callback structure for FSUB
!f2py   double precision intent(in) :: T
!f2py   double precision dimension(:), intent(in) :: Y
!f2py   double precision dimension(:), intent(out) :: FTY
!f2py	call FSUB(T, Y, FTY)

! boundary conditions at left and right ends
!f2py 	intent(callback) BCSUB
		external BCSUB, BCSUBP
! 		the following lines define the callback structure for BCSUB
!f2py	double precision dimension(:), intent(in) :: YA
!f2py	double precision dimension(:), intent(in) :: YB
!f2py	double precision dimension(:), intent(out) :: BCA
!f2py	double precision dimension(:), intent(out) :: BCB
!f2py	call BCSUB(YA, YB, BCA, BCB)

! derivative of functions with respect to y
!f2py 	intent(callback) DFDY
		external DFDY,DFDYP
! 		the following lines define the callback structure for DFDY
!f2py	double precision intent(in) :: T
!f2py	double precision dimension(:), intent(in) :: Y
!f2py	double precision dimension(:,:), intent(out) :: PD
!f2py	call DFDY(T, Y, PD)

! derivative of boundary conditions with respect to y
!f2py 	intent(callback) DBCDY
    	EXTERNAL DBCDY,DBCDYP

! 		the following lines define the callback structure for DBCDY
!f2py 	double precision dimension(:), intent(in) :: YA
!f2py 	double precision dimension(:), intent(in) :: YB
!f2py 	double precision dimension(:,:), intent(out) :: PDYA
!f2py 	double precision dimension(:,:), intent(out) :: PDYB
!f2py 	call DBCDY( YA, TB, PDYA, PDYB)


	! intermediate solution variables
		TYPE(BVP_SOL) :: init_sol
		TYPE(BVP_SOL) :: output_sol


		init_sol = sol_from_params(NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, &
		MXNSUB_in ,X_in , Y_in, PARAMETERS_in, work_in, iwork_in)

		output_sol = BVP_SOLVER(init_sol,FSUB,FSUBP, BCSUB,BCSUBP, singular, SINGULARTERM,METHOD,TOL,hasDFDY, &
		DFDY,DFDYP, hasDBCDY, DBCDY,DBCDYP, TRACE,STOP_ON_FAIL,YERROR)



		call mparams_f_sol(output_sol)


	   ! calling terminate should not be necessary because python should take care of all variables that have been passed
end subroutine

subroutine BVP_EVAL_WRAP(eval_derivative,points,npoints,NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, MXNSUB_in &
, n ,X_in , ny1, ny2, Y_in, np,  PARAMETERS_in, work_in , nwork, iwork_in, niwork)
	!   Input arguments:

        LOGICAL :: eval_derivative
        DOUBLE PRECISION, DIMENSION(npoints) :: points ! current mesh.
        integer :: npoints
	! replacement for SOL input
    	! output arguments:
		INTEGER :: NODE_in,NPAR_in ! number of ODEs; number of unknown parameters.
   		INTEGER :: LEFTBC_in ! number of left boundary conditions.
   		INTEGER :: NPTS_in,INFO_in ! number of points in current mesh; indicator of
                         ! success (INFO=0) or failure (INFO=-1) of computation.
    	INTEGER :: MXNSUB_in ! maximum number of subintervals allowed.

	    integer :: n, ny1, ny2, np, nwork, niwork

  		DOUBLE PRECISION, DIMENSION(n) :: X_in ! current mesh.
!f2py 	integer intent (hide), depend (X_in) :: n = len(X_in)

    	DOUBLE PRECISION, DIMENSION(ny1,ny2) :: Y_in ! current solution - Ith
!f2py 	integer intent (hide), depend (Y_in) :: ny1 = shape(Y_in, 0)
!f2py 	integer intent (hide), depend (Y_in) :: ny2 = shape(Y_in, 1)
	!   column is solution approximation at Ith mesh point in first NODE locations;
	!   last NPAR locations contain approximation to unknown parameters.

    	DOUBLE PRECISION, DIMENSION(np) :: PARAMETERS_in ! Unknown parameters.
!f2py 	integer intent (hide), depend (PARAMETERS_in) :: np = len(PARAMETERS_in)

    	INTEGER, DIMENSION(niwork) :: IWORK_in ! integer workspace.


    	DOUBLE PRECISION, DIMENSION(nwork) :: WORK_in ! double precision workspace.

	!   The workspaces are used to save information that can be accessed by SOL_EVAL
	!   for evaluation of the numerical solution.


	! intermediate solution variables
		TYPE(BVP_SOL) :: input_sol

		input_sol = sol_from_params(NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, &
		MXNSUB_in ,X_in , Y_in, PARAMETERS_in, work_in, iwork_in)

        if (ALLOCATED(EVALUATED)) DEALLOCATE(EVALUATED)
        ALLOCATE(EVALUATED(input_sol%NODE, npoints))

        if (eval_derivative) THEN
                if (ALLOCATED(EVALUATED_D)) DEALLOCATE(EVALUATED_D)
                ALLOCATE(EVALUATED_D(input_sol%NODE, npoints))

                CALL BVP_EVAL(input_sol, points, EVALUATED, EVALUATED_D)
    	ELSE
    	        CALL BVP_EVAL(input_sol, points, EVALUATED)
    	END IF

	   ! calling terminate should not be necessary because python should take care of all variables that have been passed
end subroutine


subroutine bvp_extend_wrap(NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, MXNSUB_in &
, n ,X_in , ny1, ny2, Y_in, np,  PARAMETERS_in, work_in , nwork, iwork_in, niwork, ANEW,BNEW,ORDER,P,MAX_NUM_SUBINTERVALS)
	!   Input arguments:

	! replacement for SOL input
    	! output arguments:
		INTEGER :: NODE_in,NPAR_in ! number of ODEs; number of unknown parameters.
   		INTEGER :: LEFTBC_in ! number of left boundary conditions.
   		INTEGER :: NPTS_in,INFO_in ! number of points in current mesh; indicator of
                         ! success (INFO=0) or failure (INFO=-1) of computation.
    	INTEGER :: MXNSUB_in ! maximum number of subintervals allowed.

	    integer :: n, ny1, ny2, np, nwork, niwork

  		DOUBLE PRECISION, DIMENSION(n) :: X_in ! current mesh.
!f2py 	integer intent (hide), depend (X_in) :: n = len(X_in)

    	DOUBLE PRECISION, DIMENSION(ny1,ny2) :: Y_in ! current solution - Ith
!f2py 	integer intent (hide), depend (Y_in) :: ny1 = shape(Y_in, 0)
!f2py 	integer intent (hide), depend (Y_in) :: ny2 = shape(Y_in, 1)
	!   column is solution approximation at Ith mesh point in first NODE locations;
	!   last NPAR locations contain approximation to unknown parameters.

    	DOUBLE PRECISION, DIMENSION(np) :: PARAMETERS_in ! Unknown parameters.
!f2py 	integer intent (hide), depend (PARAMETERS_in) :: np = len(PARAMETERS_in)

    	INTEGER, DIMENSION(niwork) :: IWORK_in ! integer workspace.


    	DOUBLE PRECISION, DIMENSION(nwork) :: WORK_in ! double precision workspace.

	!   The workspaces are used to save information that can be accessed by SOL_EVAL
	!   for evaluation of the numerical solution.
	   DOUBLE PRECISION :: ANEW, BNEW
	   INTEGER :: order, MAX_NUM_SUBINTERVALS
	   DOUBLE PRECISION, DIMENSION(np) :: P ! double precision workspace.

	! intermediate solution variables
		TYPE(BVP_SOL) :: input_sol, output_sol

	    NULLIFY(input_sol%X)
	    NULLIFY(input_sol%Y)
	    NULLIFY(input_sol%PARAMETERS)
	    NULLIFY(input_sol%WORK)
	    NULLIFY(input_sol%IWORK)

	    NULLIFY(output_sol%X)
	    NULLIFY(output_sol%Y)
	    NULLIFY(output_sol%PARAMETERS)
	    NULLIFY(output_sol%WORK)
	    NULLIFY(output_sol%IWORK)

		input_sol = sol_from_params(NODE_in,NPAR_in, LEFTBC_in ,NPTS_in,INFO_in, &
		MXNSUB_in ,X_in , Y_in, PARAMETERS_in, work_in, iwork_in)

		output_sol = BVP_EXTEND(input_sol, ANEW, BNEW, order, P, MAX_NUM_SUBINTERVALS)
		IF (.NOT. (NPAR_in > 0)) THEN
			NULLIFY (output_sol%PARAMETERS)
		END IF
		CALL mparams_f_sol(output_sol)

end subroutine
end module