module parameters
   integer, parameter :: L=200, d=2*L
   integer, parameter :: events=int(10*L)!5000
   integer  :: bc=1   ! 0=open bc     1=abc   !!!-1=pbc
   double precision :: discr=0.2!0.001*L !0.0006*L
   double precision :: delta=-0.005, mu, muc=-2., u
   double precision :: z=1., ymu=1., kappa1=0., kappa2=3., gamm=1.
end module


!  Bogoliubov transformation in homogenous system with Aperiodic boundary condition


program kitaev

   use parameters

   double complex, dimension(:,:), allocatable :: Ham, corr
   double complex, dimension(:), allocatable :: gs

   integer :: i, j, m, limit
   double precision :: t0=0.
   double precision :: energy(d), deltaL, obs, obs0

   obs=0.
   

   open(1, file='ddissip.dat', status='old')
   !open(2, file='data-kitaev', status='old')
   
   !write(2,*) '#L     delta     mu     u     Eo'


   allocate (gs(d))
   allocate (Ham(d,d)) 



   u = gamm * L**(-z)
!   u = 0.1
!   theta = t * L**(-z)
   mu = kappa1/( L**(ymu) ) + muc 
!   mu = -5.

   call definition(Ham)
   call diagonalize(Ham, energy, gs, deltaL)

   deallocate (gs)

!   do i=1, L
!      if (energy(i)<0.) obs = obs + energy(i)                 !Ham(:,1)
!   enddo

   allocate (corr(d,d))
   corr=Ham
  
   do i=1, d
      do j=1, d
         Ham(i,j)=dconjg(corr(j,i))
      enddo
   enddo
   


   do i=1, L
      obs=obs+energy(L+i)
   enddo
   obs = -mu/2.*L - obs/2. 

!write(*,*) energy, obs
!stop

   !write(2,*)  L, delta, mu, u, obs
   !close(2)



   corr=(0.d0,0.d0)

   ! c_j^\dagger c_m       corr(j,m)
   do i=1, L
      do j=i, L
         do m=1, L
            corr(j,i) = corr(j,i) + dconjg(Ham(m,j))*Ham(m,i)
         enddo
         if (i/=j) corr(i,j) = dconjg(corr(j,i))    
      enddo
   enddo
   !write(*,*) corr(4,3), corr(3,4), corr(2,2), corr(4,4)
   !write(*,*) limit, obs, corr

   ! c_j^\dagger c_m^\dagger       corr(j,m)
   do i=1+L, d
      do j=i, d
         do m=1, L
            corr(j,i) = corr(j,i) + dconjg(Ham(m,j-L))*dconjg(Ham(m,i))
         enddo 
         if(i.eq.j) corr(i,i) = 0.
         corr(i,j) = - corr(j,i)
      enddo
   enddo
   !write(*,*) corr(3+L,4+L), corr(4+L,3+L), corr(7+L,7+L)
   !write(*,*) limit, obs, corr

!stop

   deallocate (Ham)

!  QUENCH:   kappa1   --->   kappa2
   mu = kappa2/( L**(ymu) ) + muc 

   !obs0=0.91260989414656368/L

   !obs0=corr(L/2,L/2)  
   !obs0=corr(L/3,L/3+L/3) + corr(L/3+L/3,L/3)
   obs0=corr(L/3+L/6, L/3-L/6) + corr(L/3-L/6, L/3+L/6)
   write(*,*) obs0
   
   obs0=corr(L/2+L/8   +L, L/2-L/8   +L)    +    dconjg( corr(L/2+L/8   +L, L/2-L/8   +L) )
   write(*,*) obs0
   stop
   !obs0=corr(L/3+L/3   +L, L/3   +L)    +    dconjg( corr(L/3+L/3   +L, L/3   +L) )
   obs=0.

   do i=1,events
!write(*,*) i
      if (mod(i,events/20).eq.0) write(*,*) i
      call rungekutta(t0, corr, obs) 
!      write(1,*) t0*L**(-z), L**(1.)*obs !/ obs0
      write(1,*) t0*L**(-z), L**(0.)*(obs - obs0)
!      write(1,*) t0  *L**(-z), log(real(L))*L**(1.)*obs!/obs0 
!     if (t0*L**(-z) >   0.75274437665939331) exit  !  theta_max=0.75274437665939331
   enddo



stop
end program





subroutine definition(Ham)

   use parameters
  implicit none

  double complex, dimension( d, d) :: Ham
  integer :: i, j

  Ham = (0.d0,0.d0)

  do i=1, L-1
     Ham(i,i) = Ham(i,i) -  mu 
     Ham(i,i+1) = Ham(i,i+1) - 1.
     Ham(i+1,i) = Ham(i+1,i) - 1.

     Ham(i,i+L+1) = Ham(i,i+L+1) + delta
     Ham(i+1,i+L) = Ham(i+1,i+L) - delta  

     Ham(i+L,i+1) = Ham(i+L,i+1) - delta
     Ham(i+L+1,i) = Ham(i+L+1,i) + delta    

     Ham(i+L,i+L) = Ham(i+L,i+L) +  mu !dconjg(mu) 
     Ham(i+L,i+L+1) = Ham(i+L,i+L+1) + 1. !dconjg(1.)
     Ham(i+L+1,i+L) = Ham(i+1+L,i+L) + 1.
     
  enddo
   
  Ham(L,L) = Ham(L,L) -  mu 
  Ham(2*L,2*L) = Ham(2*L,2*L) +  mu

  if (abs(bc).eq.1) then

     Ham(L,1) = Ham(L,1) + bc*1.
     Ham(1,L) = Ham(1,L) + bc*1.

     Ham(L,1+L) = Ham(L,1+L) - bc*delta
     Ham(1,L+L) = Ham(1,L+L) + bc*delta

     Ham(L+L,1) = Ham(L+L,1) + bc*delta
     Ham(1+L,L) = Ham(1+L,L) - bc*delta 

     Ham(L+L,1+L) = Ham(L+L,1+L) - bc*1. !dconjg(1.)
     Ham(1+L,L+L) = Ham(1+L,L+L) - bc*1.    

  endif
  

end subroutine



subroutine diagonalize( Ham, energy, gs, deltaL)

  use parameters
  implicit none

  double precision :: energy( d), deltaL
  double complex :: Ham( d, d), gs( d)
  
  integer :: LWork, LRWork, Info
  double complex, dimension(:), allocatable :: Work
  double precision, dimension(:), allocatable :: RWork
    
  
  LWork  = 22*(2* d-1)
  LRWork = 3* d-2
  allocate (Work(Lwork), RWork(LRWork))
  Work   = 0.d0
  RWork  = 0

  call zheev('V','U', d, Ham,  d, energy, Work, LWork, RWork, Info)
  if (Info /= 0) stop 'Failed diagonalization'

  deallocate (Work, RWork)

  gs = Ham(:,1)
  deltaL=energy(2) - energy(1)
 
end subroutine







subroutine rungekutta(t0, corr, obs)

    use parameters
    implicit none

    integer :: point
    double complex, dimension(d,d) :: corr, k1, k2, k3, k4
    double precision ::  obs
    double precision :: t0

    call derivate( t0           , corr              , k1)
    call derivate( t0 + discr/2., corr + k1*discr/2., k2) 
    call derivate( t0 + discr/2., corr + k2*discr/2., k3)
    call derivate( t0 + discr   , corr + k3*discr   , k4)

    corr= corr + discr/6.*( k1 + 2.*k2 + 2.*k3 + k4 )
    t0 = t0 + discr
    obs=0.
    point=L/2
    !do point=1, d
       !obs=obs + corr(point,point) 
    !enddo
    obs = corr(point,point)
    !obs=corr(point,point+L/3) + corr(point+L/3,point)
    !obs=corr(point+L/8, point-L/8) + corr(point-L/8, point+L/8)
    !obs=corr(point+L/8   +L, point-L/8   +L)    +    dconjg( corr(point+L/8   +L, point-L/8   +L) )
    !obs=corr(point+L/3   +L, point   +L)    +    dconjg( corr(point+L/3   +L, point   +L) )
    !obs = (0.,1.)*(corr(point+1,point) - corr(point,point+1) ) 
    !obs=1./6.*( k1(point,point) + 2.*k2(point,point) + 2.*k3(point,point) + k4(point,point) )

return
end subroutine



subroutine derivate( t, corr, k)

   use parameters
   implicit none

   integer :: i, j, m!, pump, loss
   double complex, dimension( d, d) :: corr, k
   double precision :: t

   !k=(0.d0,0.d0)

   do j=1,  L
      do m=j,  L
         k(j,m) = 0.
         k(j +L,m +L) = 0.


         !j=pump   and   m=pump
         !k(j,m) = k(j,m) - u * corr(j,m) 
         !if(j.eq.m) k(j,m) = k(j,m) + u 
         !k(j +L,m +L) = k(j +L,m +L) - u * corr(j +L,m +L)
         
         !j=loss   and   m=loss
         k(j,m) = k(j,m) - u * corr(j,m)
         k(j +L,m +L) = k(j +L,m +L) - u * corr(j +L,m +L)


         if (j < L) then

            k(j,m) = k(j,m) - (0.,1.) * corr(j+1,m) + (0.,1.) * delta * dconjg(corr(m +L,j+1 +L))
            k(j +L,m +L) = k(j +L,m +L) - (0.,1.) * corr(j+1 +L,m +L) - (0.,1.) * delta * corr(m,j+1)
            if ( m.eq.(j+1) ) k(j +L,m +L) = k(j +L,m +L) + (0.,1.) * delta * 1.

         else 

             if (abs(bc).eq.1) then
                k(j,m) = k(j,m) + bc*(0.,1.) * corr(1,m) - bc*(0.,1.) * delta * dconjg(corr(m +L,1 +L))
                k(j +L,m +L) = k(j +L,m +L) + bc*(0.,1.)*corr(1 +L,m +L) + bc*(0.,1.)* delta * corr(m,1)
                if ( m.eq.(1) ) k(j +L,m +L) = k(j +L,m +L) - bc*(0.,1.) * delta * 1.
             endif

         endif



         if (m < L) then

            k(j,m) = k(j,m) + (0.,1.) * corr(j,m+1) - (0.,1.) * delta * corr(j +L,m+1 +L)
            k(j +L,m +L) = k(j +L,m +L) - (0.,1.) * corr(j +L,m+1 +L) + (0.,1.) * delta * corr(j,m+1)

         else 

             if (abs(bc).eq.1) then
                k(j,m) = k(j,m) - bc*(0.,1.) * corr(j,1) + bc*(0.,1.) * delta * corr(j +L,1 +L)
                k(j +L,m +L) = k(j +L,m +L) + bc*(0.,1.)*corr(j +L,1 +L) - bc*(0.,1.)* delta * corr(j,1)
             endif

         endif



         if (j > 1) then

            k(j,m) = k(j,m) - (0.,1.) * corr(j-1,m) - (0.,1.) * delta * dconjg(corr(m +L,j-1 +L))
            k(j +L,m +L) = k(j +L,m +L) - (0.,1.) * corr(j-1 +L,m +L) + (0.,1.) * delta * corr(m,j-1)
            if ( m.eq.(j-1) ) k(j +L,m +L) = k(j +L,m +L) - (0.,1.) * delta * 1.

         else 

             if (abs(bc).eq.1) then
                k(j,m) = k(j,m) + bc*(0.,1.) * corr(L,m) + bc*(0.,1.) * delta * dconjg(corr(m +L,L +L))
                k(j +L,m +L) = k(j +L,m +L) + bc*(0.,1.)*corr(L +L,m +L) - bc*(0.,1.)* delta * corr(m,L)
                if ( m.eq.(L) ) k(j +L,m +L) = k(j +L,m +L) + bc*(0.,1.) * delta * 1.             
             endif

         endif



         if (m > 1) then

            k(j,m) = k(j,m) + (0.,1.) * corr(j,m-1) + (0.,1.) * delta * corr(j +L,m-1 +L)
            k(j +L,m +L) = k(j +L,m +L) - (0.,1.) * corr(j +L,m-1 +L) - (0.,1.) * delta * corr(j,m-1)

         else 

             if (abs(bc).eq.1) then
                k(j,m) = k(j,m) - bc*(0.,1.) * corr(j,L) - bc*(0.,1.) * delta * corr(j +L,L +L)
                k(j +L,m +L) = k(j +L,m +L) + bc*(0.,1.)*corr(j +L,L +L) + bc*(0.,1.)*delta * corr(j,L)
             endif

         endif

 
         k(j +L,m +L) = k(j +L,m +L) - 2. * (0.,1.) * mu * corr(j +L,m +L)

         
      
         k(m +L,j +L) = - k(j +L,m +L)

         if (m.eq.j) then
            k(j +L,m +L) = 0.
            !k(j,m) = 0.5 * ( k(j,m) + dconjg(k(j,m)) )
         else
            k(m,j) = dconjg(k(j,m))
         endif
           

      enddo
   enddo 

   
     

!write(*,*) k(7,3), dconjg(k(3,7)), k(8,8), k(2,2)
!write(*,*) k(7+L,3+L), dconjg(k(3+L,7+L)), k(8+L,8+L), k(2+L,2+L)

return
end subroutine
