      subroutine gauleg(x1,x2,x,w,n)
c Routine from Numerical Recipes (Press, Flannery, Teukolsky and Vetterling,
c    1986, Cambridge University Press, ISBN 0 521 30811 9)
c
c Given the lower and upper limits of integration (X1, X2) and the number
c  of intervals (N), this routine returns arrays X and W of length N,
c  containing the abscissas and weights of the Gauss-Legendre N-point
c  quadrature formula.
c
      implicit real*4 (a-h,o-z)
      REAL*4        X1,X2,X(N),W(N)     
      parameter (eps=3.e-7)
C
      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)
      do i=1,m
        z=cos(3.141592654*(i-.25)/(n+.5))
        z1 = 0.0
        do while (abs(z-z1).gt.eps)
          p1=1.0
          p2=0.0
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
          end do
          pp=n*(z*p1-p2)/(z*z-1.0)
          z1=z
          z=z1-p1/pp
        end do
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.0*xl/((1.0-z*z)*pp*pp)
        w(n+1-i)=w(i)
      end do
      return
      end
      
