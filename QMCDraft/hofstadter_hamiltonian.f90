program generate_hofstadter_dispersion_matrix_for_w2dynamics

  implicit none

  integer, parameter :: p=1
  integer, parameter :: q=3
  integer, parameter :: ksteps=200
  real(kind=8), parameter :: pi=3.141592653589793238440d0
  real(kind=8), parameter :: t=0.250d0
  real(kind=8), parameter :: t1= 0.0d0
  real(kind=8), parameter :: t2= 0.0d0

  complex(kind=8), parameter :: Xi=dcmplx(0,1.0d0)

   !epsilon(k) Matrices
  complex(kind=8), dimension(:,:), allocatable :: epsmat 
  complex(kind=8), dimension(:,:), allocatable :: TMatrix 
  complex(kind=8), dimension(:,:), allocatable :: TPrimeMatrix 
  complex(kind=8), dimension(:,:), allocatable :: TPrimePrimeMatrix 

  character(len=1000) :: FormatString
  integer :: ikx,iky,l,lp, i, j
  integer :: kstepsX, kstepsY
  real(kind=8) :: kx,ky
  real(kind=8) :: dkx,dky
  real(kind=8) :: B

  !Allocate and initialize matrices
  allocate(epsmat(1:q,1:q))
  allocate(TMatrix(1:q,1:q))
  allocate(TPrimeMatrix(1:q,1:q))
  allocate(TPrimePrimeMatrix(1:q,1:q))
  do l=1,q
     do lp=1,q
        epsmat(l,lp)=dcmplx(0.0d0,0.0d0)
        TMatrix(l,lp)=dcmplx(0.0d0,0.0d0)
        TPrimeMatrix(l,lp)=dcmplx(0.0d0,0.0d0)
        TPrimePrimeMatrix(l,lp)=dcmplx(0.0d0,0.0d0)
     enddo
  enddo


  write(FormatString,*)2*q
  FormatString='('//trim(adjustl(FormatString))//'f22.17)'
  !initialize output data of epsilon matrix
  open(10,file="Hk_Hofstadter.dat",form="formatted",status="replace")
  write(10,'(3I6)')ksteps*(ksteps/q),q,q


  B = dfloat(p)/dfloat(q)


  !start Sampling of kx and ky in the better way (14.3.23)
  kstepsX = ksteps
  kstepsY = ksteps                                                         !Maybe change this .........  if kstepsY = ksteps/q then have same Genauigkeit, but less steps. 
                                                                           !maybe better to have better Genauigkeit by kstepsY = ksteps since exp(k*q) oscillates more
  dkx = (2 * Pi)/dfloat(kstepsX)
  dky = (2 * Pi)/dfloat(kstepsY * q)                                 !Make sure to have float division 

  do ikx = 1, kstepsX                                                !get kx from the running variable ikx
    kx = -Pi + dfloat(ikx - 1) * dkx
    do iky = 1, kstepsY
     ky = -Pi/dfloat(q) + dfloat(iky - 1) * dky



   !Construct dispersion matrix for given values of kx and ky
      TMatrix(:,:) = 0.0													!Fill the complete matrix with 0 to initialize it
      do i = 1, q									!last iterate over i and j, so over the matrix indices from 0 to (q-1), i is row, j is column
            do j = 1, q
               if (i == j) then
                  TMatrix(i,j) = TMatrix(i,j)  + 2 * t * cos(kx + (i-1) * (2 * pi * B))		!diagonale
               END IF
               if ((i == q) .and. (j == 1)) then
                  TMatrix(i,j) = TMatrix(i,j) + t*exp(- Xi*ky*q)									!bottom left
               END IF
               if ((i == 1) .and. (j == q)) then
                  TMatrix(i,j) = TMatrix(i,j) + t*exp(Xi*ky*q)									!top right
               END IF
               !(q > 1) since do not have this for q == 1, but in q = 2 have this in corner with other term, else only that term
               if (((i == (j + 1)) .or. (i == (j - 1))) .and. (q > 1)) then									
                  TMatrix(i,j) = TMatrix(i,j) + t 										!next to diagonale
               END IF
            end do
      end do	


      !Now build the TPrimeMatrix 
      TPrimeMatrix(:,:) = 0.0													!Fill the complete matrix with 0 to initialize it
      do i = 1, q									!last iterate over i and j, so over the matrix indices from 0 to (q-1), i is row, j is column
            do j = 1, q
               if ((i == q) .and. (j == 1)) then
                  TPrimeMatrix(i,j) = TPrimeMatrix(i,j) + 2 * t1*exp(-Xi*ky*q) * cos(kx + (2 * pi * B * (q - 0.5)))						!bottom left
               END IF
               if ((i == 1) .and. (j == q)) then
                  TPrimeMatrix(i,j) = TPrimeMatrix(i,j) + 2 * t1*exp(Xi*ky*q)	* cos(kx - pi * B)								!top right
               END IF
               if (i == (j - 1)) then									
                  TPrimeMatrix(i,j) = TPrimeMatrix(i,j) + 2 * t1 * cos(kx + (2 * pi * B * ((i-1) + 0.5))) 		 								!upper next to diagonale
               END IF
               if (i == (j + 1)) then									
                  TPrimeMatrix(i,j) = TPrimeMatrix(i,j) + 2 * t1 * cos(kx + (2 * pi * B * ((i-2) + 0.5))) 		 								!lower next to diagonale
               END IF
            end do
      end do	

      
      !Now build the TPrimePrimeMatrix 
      TPrimePrimeMatrix(:,:) = 0.0													!Fill the complete matrix with 0 to initialize it
      do i = 1, q									!last iterate over i and j, so over the matrix indices from 0 to (q-1), i is row, j is column
            do j = 1, q
               if (i == j) then
                  TPrimePrimeMatrix(i,j) = TPrimePrimeMatrix(i,j)  + 2 * t2 * cos(2*kx + (i-1) * (4 * pi * B))		!diagonale
               END IF
               if (((i == q) .and. (j == 2)) .or. ((i == q - 1) .and. (j == 1))) then
                  TPrimePrimeMatrix(i,j) = TPrimePrimeMatrix(i,j) + t2*exp(- Xi*ky*q)									!bottom left
               END IF
               if (((i == 1) .and. (j == q - 1)) .or. ((i == 2) .and. (j == q))) then
                  TPrimePrimeMatrix(i,j) = TPrimePrimeMatrix(i,j) + t2*exp(Xi*ky*q)									!top right
               END IF
               !(q > 2) since do not have this for q == 2, but in q = 3 have this in corner with other term
               if (((i == (j + 2)) .or. (i == (j - 2))) .and. (q > 2)) then									
                  TPrimePrimeMatrix(i,j) = TPrimePrimeMatrix(i,j) + t2 										!next to next to diagonale
               END IF
            end do
         end do


         epsmat = TMatrix + TPrimeMatrix + TPrimePrimeMatrix


         write(10,'(3f17.10)')kx, &              !writes kx
         ky, &                             !writes ky
         0.0d0                                                                      !writes kz = 0 since in 2D
         !Now finally write epsmat into the file (real and imag part one after another)
         do l=1,q
            write(10,FormatString)(dreal(epsmat(l,lp)),dimag(epsmat(l,lp)),lp=1,q)                 !a do in one row to enforce that one row in matrix is one row in file
         end do
      end do
   end do




!OLD FORMAT to describe epsilon matrix
!  do ikx=-ksteps+1,ksteps
!     kx=Pi*dfloat(ikx)/dfloat(ksteps)
!     do iky=-ksteps/q+1,ksteps/q
!        ky=Pi*dfloat(iky)/dfloat(ksteps)
        !Construct dispersion matrix for given valuex of kx and ky
!        if (q.gt.2) then
!           epsmat(1,1)=dcmplx(-2.0d0*t*dcos(kx),0.0d0)
!           epsmat(1,2)=dcmplx(-t,0.0d0)
!           epsmat(1,q)=dcmplx(-t*dcos(dfloat(q)*ky),0.0d0)+ &
!                dcmplx(0.0d0,-t*dsin(dfloat(q)*ky))
!           epsmat(q,q)=dcmplx(-2.0d0*t* &
!                dcos(kx+dfloat(q-1)*B*2.0d0*pi),0.0d0)
!           epsmat(q,q-1)=dcmplx(-t,0.0d0)
!           epsmat(q,1)=dcmplx(-t*dcos(dfloat(q)*ky),0.0d0)- &
!                dcmplx(0.0d0,-t*dsin(dfloat(q)*ky))               
!           do l=2,q-1
!              epsmat(l,l)=dcmplx(-2.0d0*t* &
!                   dcos(kx+dfloat(l-1)*B*2.0d0*Pi),0.0d0)
!              epsmat(l,l-1)=dcmplx(-t,0.0d0)
!              epsmat(l,l+1)=dcmplx(-t,0.0d0)
!           enddo
!        elseif (q.eq.2) then
!           epsmat(1,1)=dcmplx(-2.0d0*t*dcos(kx),0.0d0)
!           epsmat(1,2)=dcmplx(-t-t*dcos(dfloat(L)*ky),0.0d0)+ &
!                dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
!           epsmat(2,2)=dcmplx(-2.0d0*t*dcos(kx+B*2.0d0*Pi),0.0d0)
!           epsmat(2,1)=dcmplx(-t-t*dcos(dfloat(2)*ky),0.0d0)- &
!                dcmplx(0.0d0,-t*dsin(dfloat(2)*ky))
!        else
!           epsmat(1,1)=-2*t*(dcos(kx)+dcos(ky))
!        endif
!        write(10,'(3f17.10)')dfloat(ikx+ksteps-1)/dfloat(2*ksteps), &
!             dfloat(iky+ksteps/q-1)/dfloat(2*(ksteps/q)),0.0d0
!        do l=1,q
!           write(10,FormatString)(dreal(epsmat(l,lp)),dimag(epsmat(l,lp)),lp=1,q)
!        enddo
!     enddo
!  enddo

end program generate_hofstadter_dispersion_matrix_for_w2dynamics
