*
*  This file separates some special functions from the lapack code
*  to avoid overeager intraprocedural optimization from the intel
*  compiler.   The point of these functions was to force their arguments
*  to be trimmed to the size they would have in memory.   Floating point
*  registers actually carry more bits, which confuses LAPACK's initialization
*  where it tests floating point limits.    
*
* By putting these functions in a separate file, the clever opts that
* could keep values in registers are defeated.

************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
