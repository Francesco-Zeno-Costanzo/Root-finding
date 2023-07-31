      program testsolve
      implicit real*8(a-h,o-z)
      real*8, dimension(:,:), allocatable :: a
      real*8, dimension(:), allocatable :: x, b, y
      real*8, dimension(:), allocatable :: d1, d2, d3, x1, b1, y1
      
      N = 1000 !dimensioni matrice quadrata
      
C=========================================================================
C Test successive over-relaxation
C=========================================================================
      
      allocate(a(N,N), x(N), b(N), y(N))
      
      omega = 0.5	!fattore di rilassamento
      tol = 1E-10	!tolleranza
      
      do i = 1, N
          !creo matrice sparsa
          a(i, i) = dcos(i*1.d0)
          a(i, (i-1)+1) = dsin(i*1.d0)
          a((i-1)+1, i) = dsin(2.6d0*i)+ dcos(5.d0*i)
          a(i, (i-90)+90) = 2*dsin(i*1.d0)+ dcos(2.d0*i)
          a((i-150)+150, i) = dsin(9.6d0*i)+3.1* dcos(5.d0*i)

          x(i) = dcos(5.d0*i)*dcos(i*1.d0)	!soluzione del sistema

      enddo
      
      do i = 1, N
          b(i) = sum(a(i,:)*x) !termine noto creato a posteriori
      enddo
      
      !se va tutto bene y sarà vicino a x entro tol
      call SOR(N, a, b, y, tol, omega)
      
      print*, "test successive over-relaxation: "
      print*, sqrt(sum((x-y)**2))
      
C========================================================================
C Test matrice tridiagonale
C========================================================================
      
      allocate(d1(N), d2(N), d3(N), x1(N), b1(N), y1(N))
      
      do i = 1, N
         d1(i) = dcos(1.d0*i)			!diagonale
         d2(i) = dsin(2.d0*i)+ 2.4*cos(3.4*i)	!diagonale superiore
         d3(i) = 0.5*dsin(3.d0*i) + 2*sin(5.d-1*i)	!diagonale inferiore
         x1(i) = dcos(5.d0*i)*cos(i*1.d0)		!soluzione del sistema
      enddo
      
      !termine noto creato a posteriori
      b1(1) = d1(1)*x1(1)+d2(1)*x1(2)
      do i = 2, N-1
         b1(i) = d3(i-1)*x1(i-1)+d1(i)*x1(i)+d2(i)*x1(i+1)
      enddo
      b1(n) = d3(n-1)*x1(n-1)+d1(n)*x1(n)
      
      !se va tutto bene y1 sarà abbastanza vicino a x1
      call tri(N, d1, d2, d3, b1, y1)
      
      print*, "test matrice tridiagonale:"
      print*, sqrt(sum((x1-y1)**2))

      deallocate(d1, d2, d3, x1, b1, y1, x, b, y, a)
      
      end program testsolve
      
C========================================================================
C soubroutine per la risoluzione di sistemi con matrici tridiagonali
C========================================================================

      subroutine tri(N, diag, sup_diag, inf_diag, b, x)
	implicit real*8 (a-h,o-z)
	real*8, dimension(N) :: diag, sup_diag, inf_diag
	real*8, dimension(N) :: d_ii, d_up, d_lo, b, x
	
	! we write the information on these arrays
	! so that the subroutine doesn't modify the old ones
	d_ii = diag       !d_ii is the diagonal of the matrix
	d_up = sup_diag   !d_up is the upper diagonal
	d_lo = inf_diag   !d_lo is the lower diagonal
	
	! make the matrix an upper triangular by canceling d_lo
	do i = 2, N
	    if(d_ii(i-1) == 0.0) then
              print *,"Division by zero, non invertible matrix"
              return
          endif
          a = d_lo(i-1)/d_ii(i-1)
          d_ii(i) = d_ii(i) - a*d_up(i-1)
          b(i) = b(i) - a*b(i-1)
	enddo

	! at this point it is easy to find x_n given that
	! in the last row of the matrix, there is only
	! the element left on the diagonal which is non-zero
	if(d_ii(N) == 0.0) then
          print *,"Division by zero, non invertible matrix"
          return
      endif
      
      x(N) = b(N)/d_ii(N)
      	
      !and now going backwards I can find all the other x_i
      do i = N-1, 1, -1
          b(i) = b(i) - d_up(i)*x(i+1)
          if(d_ii(i) == 0.0) then
              print *,"Division by zero, non invertible matrix"
              return
          endif
          x(i) = b(i)/d_ii(i)
      enddo

	return
	end


C========================================================================
C soubroutine algoritmo Test successive over-relaxation
C========================================================================
   
      subroutine SOR(N, A, b, x, tol, omega)
      real*8, dimension(N, N) :: A
      real*8, dimension(N) :: b, x, prod
      real*8 :: tol, omega, sigma, res
      
      !calcolo residuo iniziale
      res = 0    
      do i = 1, N
          prod(i) = sum(A(i,:)*x)
      enddo
      res = sqrt(sum((prod-b)**2))
      
      !iter = 0
      do while (res>tol)
          do i = 1, N
              sigma = 0
              do j = 1, N
                  if (j/=i) then
                      sigma = sigma + A(i, j)*x(j)
                  endif 
              enddo
              x(i) = (1-omega)*x(i) + (omega/A(i, i))*(b(i)-sigma)
          enddo
          
          do i = 1, N
              prod(i) = sum(A(i,:)*x)
          enddo
          res = sqrt(sum((prod-b)**2))
          !iter = iter + 1
          !print*, iter, res
      enddo
      
      
      return
      end
