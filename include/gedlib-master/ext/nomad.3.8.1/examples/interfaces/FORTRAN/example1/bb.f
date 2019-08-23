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
