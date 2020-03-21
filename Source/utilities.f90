function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
  real ( kind = 8 ) alnorm
  real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
  real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
  real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
  real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
  real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
  real ( kind = 8 ), parameter :: con = 1.28D+00
  real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
  real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
  real ( kind = 8 ), parameter :: ltone = 7.0D+00
  real ( kind = 8 ), parameter :: p = 0.398942280444D+00
  real ( kind = 8 ), parameter :: q = 0.39990348504D+00
  real ( kind = 8 ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = 8 ), parameter :: utzero = 18.66D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
  end if

  return
end



SUBROUTINE ran1sub(idum,x) 
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL(8) :: x 
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647,IQ=127773,IR=2836
REAL(8), SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
if (idum <= 0 .or. iy < 0) then
     am=nearest(1.0d-0,-1.0d-0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1
end if
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))
k=iy/IQ
iy=IA*(iy-k*IQ)-IR*k
if (iy<0) iy=iy+IM
x=am*ior(iand(IM,ieor(ix,iy)),1)
END SUBROUTINE ran1sub


SUBROUTINE ran1(idum,vv,n)
 
implicit none

external ran1sub

REAL(8) vv(n),x
INTEGER kk,n,idum

vvloop: do kk = 1,n
     call ran1sub(idum,x)
     vv(kk) = x
end do vvloop

idum = idum+2

END SUBROUTINE ran1
 


!------------------------


       SUBROUTINE GASDEV(idum,vv,n)

! This function returns a normally distributed deviate with
! zero mean and unit variance, using ran1(IDUM) as the
! source of uniform deviates

       implicit none

       integer iset,n,k,gotoc,idum
       real(8) v1,v2,r,fac,gset,gasdev1
       real(8) temp1(2),vv(n),ix1,ix2

       external ran1

!       write(*,*) 'want',n,'random numbers'

       vvloop: do k = 1,n

       iset = 0
       v1 = 0.0d-0
       v2 = 0.0d-0
       r = 0.0d-0
       fac = 0.0d-0
       gset = 0.0d-0
       gasdev1 = 0.0d-0
       ix1=0.0d-0
       ix2=0.0d-0
       gotoc = 0

       if (iset.eq.0) then


1         call ran1(idum,temp1,2)

          v1 = (2.0d-0)*temp1(1)-1.0d-0 
          v2 = (2.0d-0)*temp1(2)-1.0d-0

!          write(*,*) 'v1',v1
!          write(*,*) 'v2',v2
!          write(*,*) 'idum',idum

          r = v1**2.0d-0+v2**2.0d-0
          if (r.ge.1.0d-0.or.r.eq.0.0d-0) then 
             gotoc = gotoc + 1
             if (gotoc.gt.n) then
                !write(*,*) 'error in gasdev'
             end if
             go to 1
          end if

          fac = dsqrt(-2.0d-0*dlog(r)/r)
          gset = v1*fac
          gasdev1 = v2*fac
          iset = 1
       else
          gasdev1 = gset
          iset = 0
       end if
       
       vv(k)=gasdev1
!       write(*,*) 'random normal number',k,vv(k)

       end do vvloop

       return
       end SUBROUTINE GASDEV





         








 
         








 
