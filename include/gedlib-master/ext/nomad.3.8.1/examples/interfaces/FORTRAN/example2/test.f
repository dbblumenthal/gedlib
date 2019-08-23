      program test

      real*8 x(5)
      real*8 lb(5)
      real*8 ub(5)

      x(1) = 0
      x(2) = 0
      x(3) = 0
      x(4) = 0
      x(5) = 0

      lb(1) = -6
      lb(2) = -6
      lb(3) = -6
      lb(4) = -6
      lb(5) = -6

      ub(1) = 5
      ub(2) = 6
      ub(3) = 7
      ub(4) = 1e20
      ub(5) = 1e20

      call nomad( 5 , 3 , x , lb , ub , -1 , 0 )

      write(*,*) 'solution:'
      do 10 i = 1, 5
         write(*,*) '   ', x(i)
  10  continue

      stop
      end




      subroutine bb(x,fx)
      real*8  x(5)
      real*8  fx(3)
      integer i

      fx(1) = x(5)
      fx(2) = 0.0
      fx(3) = 0.0

      do 10 i = 1, 5
         fx(2) = fx(2) + (x(i)-1)*(x(i)-1)
         fx(3) = fx(3) + (x(i)+1)*(x(i)+1)
  10  continue

      fx(2) = fx(2)-25
      fx(3) = 25-fx(3)

      return
      end
