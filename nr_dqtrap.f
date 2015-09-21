C> \file
C> Trapezoidal integrator from Numerical Recipes

C> Returns as `s` the integral of the function `func` from `a` to `b`. 
C>
C> The parameters `EPS` can be set to the desired fractional accuracy and 
C> `JMAX` so that 2 to the power `JMAX-1` is the maximum allowed number 
C> of steps. Integration is performed by the trapezoidal rule.
C>
C> @author W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery: 
C> Numerical Recipes in Fortran 77
C>
      SUBROUTINE qtrap(func,a,b,s)
      INTEGER JMAX
      DOUBLE PRECISION a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-8, JMAX=20)
CU    USES trapzd
      INTEGER j
      DOUBLE PRECISION olds
      olds=-1.d30
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (j.gt.5) then
          if (abs(s-olds).lt.EPS*abs(olds).or.
     *(s.eq.0.d0.and.olds.eq.0.d0)) 
     *return
        endif
        olds=s
11    continue
      STOP 'too many steps in qtrap'
      END
