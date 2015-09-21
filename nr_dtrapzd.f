C> \file
C> Trapezoidal rule from Numerical Recipes

C> This routine computes the `n`th stage of refinement of an extended 
C> trapezoidal rule.
C>
C> `func` is input as the name of the function to be integrated between 
C> limits `a` and `b`, also input. When called with `n=1`, the routine returns 
C> as `s` the crudest estimate of int(f(x),x=a..b) Subsequent calls with 
C> `n=2,3,...` (in that sequential order) will improve the accuracy of `s` by 
C> adding 2^(n-2) additional interior points. `s` should not be modified 
C> between sequential calls.
C>
C> @author W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery: 
C> Numerical Recipes in Fortran 77
C>
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
