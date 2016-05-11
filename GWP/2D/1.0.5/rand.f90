          program test_random_number
            REAL*8 :: r(5,5),x(100),y(100)
         integer*4 i 
!            CALL init_random_seed()         ! see example of RANDOM_SEED
      
          call RANDom_number(x)
          call random_number(y)
          do i =1,100
          print *,x(i),y(i)
          enddo 
    
          end program
