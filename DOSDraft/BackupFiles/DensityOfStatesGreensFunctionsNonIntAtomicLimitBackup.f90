program DensityOfStates
  
  use iso_fortran_env
  use mpi_f08

  implicit none

  !---------------------------------------------Constants and parameters---------------------------------------------
  
  !\pi.
  real(real64), parameter :: pi=3.141592653589793238462643383280d0
  
  !Nearest neighbor hoping parameter t. This is set as constant parameter which fixes the energy scale.
  real(real64), parameter :: t=0.250d0

  !Logical variables to determine what should be calculated.
  logical, parameter :: calcdos=.TRUE.
  logical, parameter :: calcdoslowrankupdate=.TRUE.
  logical, parameter :: calcgretarded=.TRUE.
  logical, parameter :: calcgretardedlowrankupdate=.TRUE.
  logical, parameter :: calcgmatsubara=.TRUE.
  logical, parameter :: calcgmatsubaralowrankupdate=.TRUE.
  
  !Thresholds for equality of eigenvalues and size of weights to be considered in the low rank update.
  !The value of tolerance is determined by the intrinsic epsilon(x) function.
  real(real64) :: tolerance,tolerance2
  
  !------------------------------------------------------------------------------------------------------------------

  
  !---------------------------------------Input model and calculation parameters--------------------------------------

  !Total number of sets of input model and calculation parameters for which a calculation should be performed.
  integer :: number_input_sets

  !Model parameters for one data set:

  !Magnetic field: B = a^2 * \Phi_0 * p/q where the lattice constant "a" and the flux quantum "\Phi_0" are set to 1.
  integer :: p,q

  !(Inverse) temperature \beta=1/k_BT. Required for the calculation of Matsubara Green's function.
  real(real64) :: beta

  !Hubbard interaction if Matsubara and retarded Green's function should be calculated with the self-energy of the atomic limit.
  !For a non-interacting calculation set uhub to 0.0d0.
  real(real64) :: uhub
  
  !Calculation parameters for one data set.

  !Numbers Nkx and Nky of grid points in kx and ky direction, respectively,  in the interval [0,\pi] for the momentum integration for eigenvalues.
  integer :: Nkx,Nky
  
  !Interval I=[epsmin,epsmax] in which the density of states should be calculated.
  real(real64) :: epsmin,epsmax
  
  !Number of (equally spaced) sampling points Neps in the interval I=[epsmin,epsmax] and/or number of positive Matsubara frequencies.
  integer :: Neps
  
  !Broadening delta of the Lorentzian curve.
  real(real64) :: delta
  
  !------------------------------------------------------------------------------------------------------------------
  

  !------------------------------------------------------Output------------------------------------------------------

  !Density of states (DOS). "dospart" denotes the contribution to the DOS from a set of kx from one MPI process.
  real(real64), dimension(:), allocatable :: dospart,dos,dospartlowrankupdate,doslowrankupdate

  !Local (i.e., momentum summed) part of of the Green's function:
  !-) gretardedlocpart, gretardedloc: Local retarded Green's function calculated by using the Woodbury formula for the matrix inversion.
  !-) gretardedpart,gretarded: Local retarded Green's function calculated with plain matrix inversion (i.e., without the Woodbury formula).
  !-) gmatsubaralocpart,gmatsubaraloc: Local Matsubara Green's function calculated by using the Woodbury formula for the matrix inversion.
  !-)gmatsubarapart,gmatsubara: Local Matsubara Green's function calculated with plain matrix inversion (i.e., without the Woodbury formula).
  complex(real64), dimension(:,:,:), allocatable :: gretardedloc,gretarded,gmatsubaraloc,gmatsubara
  complex(real64), dimension(:,:,:), allocatable :: gretardedlocpart,gretardedpart,gmatsubaralocpart,gmatsubarapart
  
  !------------------------------------------------------------------------------------------------------------------


  !------------------------------------------------Internal variables------------------------------------------------

  !Diagonal and subdiagonal/supperdiagonal parts of the partial dispersion matrix [without the element (1,q) and (q,1) in the corners].
  !!!Note that later we define the dispersion with an additional "-" sign!!!
  real(real64), dimension(:), allocatable :: eps_diag,eps_subdiag
  complex(real64), dimension(:), allocatable :: eps_diag_complex,eps_subdiag_complex,eps_superdiag_complex

  !Full dispersion matrix.
  complex(real64), dimension(:,:), allocatable :: eps_full
  
  !Frequency dependent part of the diagonal of the inverse Green's function [z-\Sigma(z)].
  !z=\omega+i\delta or z=i\nu for the retarded and Matsubara Green's function, respectively.
  !\Sigma(z)=0 (non-interacting system) or \Sigma(z)=U/(4z) (for the strong coupling atomic limit case), respectively.
  real(real64), dimension(:), allocatable :: omega,nu
  complex(real64), dimension(:), allocatable :: zetaretarded,zetamatsubara
  
  !Diagonal part of the inverse of the retarded/Matsubara Green's function:
  ![\omega-\varepsilon_k-\Sigma(\omega)]_{ll} and [i\nu-\varepsilon_k-\Sigma(i\nu)]_{ll}.
  !\Sigma(z)=i\delta / 0 (non-interacting case) or \Sigma(z)=U/2+U^2(4z) (strong coupling atomic limit self-energy).
  complex(real64), dimension(:), allocatable :: ginvretarded_diag,ginvmatsubara_diag

  !Inverse of the tridiagonal part of [\omega-\varepsilon_k-\Sigma(\omega)] or [i\nu-\varepsilon_k-\Sigma(i\nu)], respectively.
  complex(real64), dimension(:,:), allocatable :: gretardedtridiag,gmatsubaratridiag

  !Inverse of the tridiagonal part of [z-\varepsilon_k-\Sigma(z)] multiplied with the vectors v=(1,...,e^(+/-I*q*k_y)).
  !The vectors v correspond to the rank one-update of the tridiagonal matrix to obtain the full dispersion matrix.
  complex(real64), dimension(:), allocatable :: gtridiagv,vdaggtridiag

  !Full inverse Green's function at a given momentum.
  complex(real64), dimension(:,:), allocatable :: ginv_full
  
  !Denominator of rank-one update of inverse of tridiagonal matrix.
  complex(real64) :: denomlowrankupdate
  
  !Eigenvectors of the partial dispersion matrix. The corresponding eigenvalues are stored in eps_diag on ouput.
  real(real64), dimension(:,:), allocatable :: eigenvectors

  !Non-degenerate eigenvalues of the partial dispersion matrix.
  real(real64), dimension(:), allocatable :: eigenvalues_not_equal
  
  !Eigenvalues and weights (i.e., squared components of the update vector v in the eigenbasis of the partial dispersion matrix) which are used for the rank-1 update.
  real(real64), dimension(:), allocatable :: eigenvalues_used,weight_used

  !Updated eigenvalues of the full (negative) dispersion matrix (i.e., after the rank-1 update).
  real(real64), dimension(:), allocatable :: eigenvalues

  !Array for precalculating the weights w_l=|(a_l,....,b_l)*(1,...,e^(I*q*ky))|^2=a_l^2+b_l^2+2*a_l*b_l*cos(q*ky).
  !(a_l,...,b_l) denotes the l-th eigenvector of the partial dispersion matrix.
  !The ky-independent parts of these term, i.e., a_l^2+b_l^2 and 2*a_l*b_l are stored in the array.
  real(real64), dimension(:,:), allocatable :: weight_precalc

  !Array for precalculating cos(q*ky).
  real(real64), dimension(:), allocatable :: precalcky

  !Array for precalculating 2*cos(q*ky) and exp(+/-i*q*k_y).
  complex(real64), dimension(:), allocatable :: precalcky2
  complex(real64), dimension(:,:), allocatable :: precalckyexp
  
  !Auxiliary variables.
  integer :: i,j,k,l,lp,m,n,count_no_update,count_equal,Nky2,isign,kxsign
  real(real64) :: kx,ky,dkx,dky,Bmag,dpq,dpNeps,normkxsum,normkysum,deltasquared,weight,new_eigenvalue,sum
  real(real64), dimension(:), allocatable :: Bmagshift,multipkx,multipky

  !Lapack variables.
  integer :: info
  integer, dimension(:), allocatable :: ipiv
  real(real64), dimension(:), allocatable :: work,eigendiff,rwork
  complex(real64), dimension(:), allocatable :: workinv,lwork

  !Variables for distributing the Nkx grid points in kx direction among the processes.
  integer :: knumber,krest,k_min,k_max,koffset

  !MPI variables.
  integer :: nprocs,myid,ierror
  
  !MPI initialization.
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)

  !Set tolerances for eigenvalues and eigenvectors of the partial dispersion matrix to be updated by the rank-1 update.
  tolerance=epsilon(0.0d0)*0.50d0
  tolerance2=tolerance*tolerance
  
  !Read model and calculation parameters from input file "Input_parameters_DOS.dat" (root process only).
  if (myid.eq.0) then
     open(10,file="InputParametersDOS.dat",form="formatted",status="unknown")
     !Read total number of input parameter sets and logical variables from the second line of the input file.
     read(10,*)
     read(10,*)number_input_sets
  endif
  !Broadcast the number of input parameter sets and logical variables which determine what should be calculated to all processes.
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call MPI_BCAST(number_input_sets,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  !Open files for output.
  if (myid.eq.0) then
     open(11,file="DensityOfStatesLowRankUpdate.dat",form="formatted",status="replace")
     open(12,file="DensityOfStates.dat",form="formatted",status="replace")
     open(13,file="RetardedGreensFunctionLowRankUpdate.dat", form="formatted",status="replace")
     open(14,file="MatsubaraGreensFunctionLowRankUpdate.dat", form="formatted",status="replace")
     open(15,file="RetardedGreensFunction.dat", form="formatted",status="replace")
     open(16,file="MatsubaraGreensFunction.dat", form="formatted",status="replace")     
  endif
  
  !Loop over all sets of input parameters.
  do n=1,number_input_sets

     if (myid.eq.0) then
        write(*,*)"Calculation for input set",n
     endif
     
     !Read n-th set of model and calculation parameters (root process only).
     if (myid.eq.0) then
        read(10,*)p,q,uhub,beta,Nkx,Nky,epsmin,epsmax,Neps,delta
     endif
     
     !Broadcast n-th set of model and calculation parameters to all processes.
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
     call MPI_BCAST(p,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(q,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(uhub,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(Nkx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(Nky,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(epsmin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(epsmax,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(Neps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     call MPI_BCAST(delta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)

     !Number of grid points for ky in the interval [-\pi,\pi).
     Nky2=2*(Nky-1)
     
     !Allocate arrays.
     allocate(Bmagshift(1:q))
     allocate(eps_diag(1:q))
     allocate (eps_diag_complex(1:q))
     allocate(eps_subdiag(1:q-1))
     allocate(eps_subdiag_complex(1:q-1))
     allocate(eps_superdiag_complex(1:q-1))
     allocate(eps_full(1:q,1:q))
     allocate(ginvretarded_diag(1:q))
     allocate(ginvmatsubara_diag(1:q))
     allocate(gretardedtridiag(1:q,1:q))
     allocate(gmatsubaratridiag(1:q,1:q))
     allocate(gtridiagv(1:q))
     allocate(vdaggtridiag(1:q))
     allocate(ginv_full(1:q,1:q))
     allocate(eigenvectors(1:q,1:q))
     allocate(eigenvalues_not_equal(1:q))
     allocate(weight_precalc(1:2,1:q))
     allocate(precalcky(1:Nky))
     allocate(precalcky2(1:Nky2))
     allocate(precalckyexp(1:2,1:Nky2))
     allocate(eigenvalues_used(1:q))
     allocate(weight_used(1:q))
     allocate(eigenvalues(1:q))
     allocate(omega(0:Neps))
     allocate(nu(0:Neps))
     allocate(zetaretarded(0:Neps))
     allocate(zetamatsubara(0:Neps))
     allocate(multipkx(1:Nkx))
     allocate(multipky(1:Nky))
     allocate(ipiv(1:q))
     allocate(work(1:2*q-2))
     allocate(lwork(1:10*q))
     allocate(rwork(1:3*q-2))
     allocate(workinv(1:10*q))
     allocate(eigendiff(1:q))
     if (calcgretardedlowrankupdate) allocate(gretardedlocpart(0:Neps,1:q,1:q))
     if (calcgretardedlowrankupdate) allocate(gretardedloc(0:Neps,1:q,1:q))
     if (calcgretarded) allocate(gretardedpart(0:Neps,1:q,1:q))
     if (calcgretarded) allocate(gretarded(0:Neps,1:q,1:q))
     if (calcgmatsubaralowrankupdate) allocate(gmatsubaralocpart(0:Neps,1:q,1:q))
     if (calcgmatsubaralowrankupdate) allocate(gmatsubaraloc(0:Neps,1:q,1:q))
     if (calcgmatsubara) allocate(gmatsubarapart(0:Neps,1:q,1:q))
     if (calcgmatsubara) allocate(gmatsubara(0:Neps,1:q,1:q))
     if (calcdoslowrankupdate) allocate(dospartlowrankupdate(0:Neps))
     if (calcdoslowrankupdate) allocate(doslowrankupdate(0:Neps))
     if (calcdos) allocate(dospart(0:Neps))
     if (calcdos) allocate(dos(0:Neps))

     
     !Define double precision variables from integers to avoid repeated type casting.
     dpq=dble(q)
     Bmag=dble(p)/dpq
     dpNeps=dble(Neps)
     normkxsum=dble(Nkx-1)
     normkysum=dble(Nky-1)

     !Calculate square of Lorentzian boradening.
     deltasquared=delta*delta
     
     !Precalculate shifts in cosine functions in dispersion realtion due to magnetic field Bmagshift=2*\pi*p*(l-1)/q.
     do l=1,q
        Bmagshift(l)=2.0d0*pi*Bmag*dble(l-1)
     enddo

     !Distribute kx grid points among the processes.
     knumber=Nkx/nprocs
     krest=mod(Nkx,nprocs)
     if (myid.lt.krest) then
        knumber=knumber+1
     endif
     if (myid.lt.krest) then
        koffset=myid*knumber
     else
        koffset=krest+myid*knumber
     endif
     k_min=koffset+1
     k_max=koffset+knumber

     !Set grid resolution (interval).
     dkx=pi/normkxsum
     dky=pi/dble(q*(Nky-1))
  
     !Precalculate cos(q*k_y) and and exp(+/-iqk_y)
     do i=1,Nky
        ky=dble(i-1)*dky
        precalcky(i)=cos(ky*dpq)
     enddo
     do i=1,Nky2
        ky=-pi/dpq+dble(i-1)*dky
        precalcky2(i)=cmplx(2.0d0*cos(ky*dpq),0.0d0,kind=real64)
        precalckyexp(1,i)=exp(cmplx(0.0d0,dpq*ky,kind=real64))
        precalckyexp(2,i)=exp(cmplx(0.0d0,-dpq*ky,kind=real64))
     enddo

     !Precalculate multiplicities of points in the kx and ky momentum grids.
     !This is required if one integrates only over the interval [0:\pi] instead of [-\pi,\pi] which is only possible for the DOS.
     multipkx(1)=0.50d0
     multipkx(Nkx)=0.50d0
     do i=2,Nkx-1
        multipkx(i)=1.0d0
     enddo
     multipky(1)=0.50d0
     multipky(Nky)=0.50d0
     do i=2,Nky-1
        multipky(i)=1.0d0
     enddo
     
     !Initialize density of states.
     do i=0,Neps
        if (calcdoslowrankupdate) dospartlowrankupdate(i)=0.0d0
        if (calcdoslowrankupdate) doslowrankupdate(i)=0.0d0
        if (calcdos) dospart(i)=0.0d0
        if (calcdos) dos(i)=0.0d0
     enddo

     !Initialize the local retarded/Matsubara Green's function G_R_loc(\omega)/G_loc(i\nu).
     do k=0,Neps
        do l=1,q
           do lp=1,q
              if (calcgretardedlowrankupdate) gretardedloc(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgretardedlowrankupdate) gretardedlocpart(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgmatsubaralowrankupdate) gmatsubaraloc(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgmatsubaralowrankupdate) gmatsubaralocpart(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgretarded) gretarded(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgretarded) gretardedpart(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgmatsubara) gmatsubara(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
              if (calcgmatsubara) gmatsubarapart(k,l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
           enddo
        enddo
     enddo
     
     !Precalculate frequency grid for z and zeta(z)=U^2/(4z) for z=\omega+i\delta (retarded Green's function) and z=i\nu (Matsubara Green's function), respectively.
     do i=0,Neps
        omega(i)=epsmin+dble(i)*(epsmax-epsmin)/dpNeps
        nu(i)=dble(2*i+1)*pi/beta
        zetaretarded(i)=cmplx(omega(i)*(1.0d0-uhub**2/(4.0d0*(omega(i)**2+delta**2))), &
             delta*(1.0d0+uhub**2/(4.0d0*(omega(i)**2+delta**2))),kind=real64)
        zetamatsubara(i)=cmplx(0.0d0,nu(i)+uhub**2/(4.0d0*nu(i)),kind=real64)                
     enddo

     !Loop of all kx values assigned to the actual MPI process.
     do i=k_min,k_max

        !For the calculation of the Green's functions the symmetries kx<->-kx and ky<->-ky cannot be used.
        !Therefore and additional sum over +/-kx for each kx-point has to be performed [apart from kx=0 (i=1) and kx=\pi (i=Nkx)].
        if ((i.eq.1).or.(i.eq.Nkx)) then
           isign=0
        else
           isign=-1
        endif

        !Sum over +/-kx.
        do kxsign=isign,0
        
           kx=dble((i-1)*(1+2*kxsign))*dkx

           !Construct diagonal part of the (negative) dispersion matrix at the given kx-point.
           call construct_diag_dispersion(q,Bmagshift,t,kx,eps_diag)
        
           !Cast the diagonal part of the (negative) dispersion matrix to type complex.
           !write(*,*)"Retarded Green's function started."
           do l=1,q
              eps_diag_complex(l)=cmplx(eps_diag(l),0.0d0,kind=real64)
           enddo

           !Loop over all frequencies \omega or i\nu, respectively.
           do k=0,Neps
              
              !The following block is evaluated if the retarded Green's function should be calculated via a low rank update.
              if (calcgretardedlowrankupdate) then
                 !Calculate the diagonal part of [zeta(\omega+i\delta)-\varepsilon_k], zeta[z]=z or zeta[z]=z-U/(4z) for the non-interacting system and the atomic limit, respectively.
                 !Notice the "+" in front of eps_diag_complex since eps_diag is defined here with an overall negative sign.
                 do l=1,q
                    ginvretarded_diag(l)=zetaretarded(k)+eps_diag_complex(l)
                 enddo                 
                 !Create identity Matrix for the right hand side of the equation [zeta(\omega+i\delta)-\varepsilon]*gtridiag=1 (=identitymatrix).
                 do l=1,q
                    do lp=1,q
                       gretardedtridiag(l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
                    enddo
                    gretardedtridiag(l,l)=cmplx(1.0d0,0.0d0,kind=real64)              
                 enddo                 
                 !Set (negative of the) first sub- and super-diagonal parts of inverse Green's function matrix.
                 do l=1,q-1
                    eps_subdiag_complex(l)=cmplx(t,0.0d0,kind=real64)
                    eps_superdiag_complex(l)=cmplx(t,0.0d0,kind=real64)
                 enddo

                 !Invert the tridiagonal part of (\omega-\varepsilon+i\delta).
                 call zgtsv(q,q,eps_subdiag_complex,ginvretarded_diag,eps_superdiag_complex,gretardedtridiag,q,info)                 
              endif
              
              !The following block is evaluated if the Matsubara Green's function should be calculated via a low rank update.
              if (calcgmatsubaralowrankupdate) then              
                 !Calculate diagonal part of [zeta(i\nu)-\varepsilon_k], zeta[z]=z or zeta[z]=z-U/(4z) for the non-interacting system and the atomic limit, respectively.
                 !Notice the "+" in front of eps_diag_complex since eps_diag is defined here with an overall negative sign.
                 do l=1,q
                    ginvmatsubara_diag(l)=zetamatsubara(k)+eps_diag_complex(l)
                 enddo                 
                 !Create identity Matrix for the right hand side of the equation [zeta(i\nu)-\varepsilon]*gtridiag=1 (=identitymatrix).
                 do l=1,q
                    do lp=1,q
                       gmatsubaratridiag(l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
                    enddo
                    gmatsubaratridiag(l,l)=cmplx(1.0d0,0.0d0,kind=real64)              
                 enddo
                 !Set (negative of the) first sub- and super-diagonal parts of inverse Green's function matrix.
                 do l=1,q-1
                    eps_subdiag_complex(l)=cmplx(t,0.0d0,kind=real64)
                    eps_superdiag_complex(l)=cmplx(t,0.0d0,kind=real64)
                 enddo                 
                 !Invert the tridiagonal part of (i\nu-\varepsilon).
                 call zgtsv(q,q,eps_subdiag_complex,ginvmatsubara_diag,eps_superdiag_complex,gmatsubaratridiag,q,info)
              endif
                            
              !Add the rank-one update due to non-tridiagonal term in the dispersion matrix.
              do j=1,Nky2
                 
                 !The following block is evaluated if the retarded Green's function should be calculated via a low rank update.                 
                 if (calcgretardedlowrankupdate) then                                        
                    !Calculate v^+*gretardedtridiag*v for v^+=(1,0,...,0,e^(I*q*ky)).
                    denomlowrankupdate=cmplx(1.0d0,0.0d0,kind=real64)+cmplx(t,0.0d0,kind=real64)* &
                         (gretardedtridiag(1,1)+gretardedtridiag(q,q)+precalcky2(j)*gretardedtridiag(1,q))
                    !Calculate the vectors gretardedtridiag*v and v^+*gretardedtridiag for v^+=(1,0,...,0,e^(I*q*ky)).
                    !Note: the first vector is multiplied with t and divided by dnomlowrankupdate.
                    do l=1,q
                       gtridiagv(l)=cmplx(t,0.0d0,kind=real64)*(gretardedtridiag(l,1)+ &
                            gretardedtridiag(l,q)*precalckyexp(2,j))/denomlowrankupdate
                       vdaggtridiag(l)=gretardedtridiag(1,l)+gretardedtridiag(q,l)*precalckyexp(1,j)
                       !write(300,*)gtridiagv(l),vdaggtridiag(l)
                    enddo
                    !Calculate the full contribution of the low-rank-update to the local retarded Green's function.
                    do l=1,q
                       do lp=1,q
                          gretardedlocpart(k,l,lp)=gretardedlocpart(k,l,lp)+(gretardedtridiag(l,lp)- &
                               gtridiagv(l)*vdaggtridiag(lp))
                       enddo
                    enddo
                 endif
                 
                 !The following block is evaluated if the Matsubara Green's function should be calculated via a low rank update.
                 if (calcgmatsubaralowrankupdate) then              
                    !Calculate v^+*gretardedtridiag*v for v^+=(1,0,...,0,e^(I*q*ky)).
                    denomlowrankupdate=cmplx(1.0d0,0.0d0,kind=real64)+cmplx(t,0.0d0,kind=real64)* &
                         (gmatsubaratridiag(1,1)+gmatsubaratridiag(q,q)+precalcky2(j)*gmatsubaratridiag(1,q))
                    !Calculate the vectors gretardedtridiag*v and v^+*gretardedtridiag for v^+=(1,0,...,0,e^(I*q*ky)).
                    !Note: the first vector is multiplied with t and divided by dnomlowrankupdate.
                    do l=1,q
                       gtridiagv(l)=cmplx(t,0.0d0,kind=real64)*(gmatsubaratridiag(l,1)+ &
                            gmatsubaratridiag(l,q)*precalckyexp(2,j))/denomlowrankupdate
                       vdaggtridiag(l)=gmatsubaratridiag(1,l)+gmatsubaratridiag(q,l)*precalckyexp(1,j)
                    enddo
                    !Calculate the full contribution of the low-rank-update to the local retarded Green's function.
                    do l=1,q
                       do lp=1,q
                          gmatsubaralocpart(k,l,lp)=gmatsubaralocpart(k,l,lp)+(gmatsubaratridiag(l,lp)- &
                               gtridiagv(l)*vdaggtridiag(lp))!*cmultipkx(i)*cmultipky(j)
                       enddo
                    enddo
                 endif

                 !The following block is evaluated if the retarded Green's function should be calculated without a low rank update.
                 if (calcgretarded) then
                    !Construct full inverse retarded Green's function.
                    do l=1,q
                       do lp=1,q
                          ginv_full(l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
                       enddo
                    enddo
                    do l=1,q-1
                       ginv_full(l,l)=zetaretarded(k)+eps_diag_complex(l)
                       ginv_full(l,l+1)=cmplx(t,0.0d0,kind=real64)
                       ginv_full(l+1,l)=cmplx(t,0.0d0,kind=real64)
                    enddo
                    ginv_full(1,1)=ginv_full(1,1)+cmplx(t,0.0d0,kind=real64)
                    ginv_full(q,q)=zetaretarded(k)+eps_diag_complex(q)+cmplx(t,0.0d0,kind=real64)
                    ginv_full(1,q)=cmplx(t,0.0d0,kind=real64)*precalckyexp(1,j)
                    ginv_full(q,1)=cmplx(t,0.0d0,kind=real64)*precalckyexp(2,j)                 
                    call zgetrf(q,q,ginv_full,q,ipiv,info)
                    call zgetri(q,ginv_full,q,ipiv,workinv,10*q,info)
                    do l=1,q
                       do lp=1,q
                          gretardedpart(k,l,lp)=gretardedpart(k,l,lp)+ginv_full(l,lp)!*cmultipkx(i)*cmultipky(j)
                       enddo
                    enddo
                 endif

                 !The following block is evaluated if the Matsubara Green's function should be calculated without a low rank update.
                 if (calcgmatsubara) then
                    !Construct full inverse retarded Green's function.
                    do l=1,q
                       do lp=1,q
                          ginv_full(l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
                       enddo
                    enddo
                    do l=1,q-1
                       ginv_full(l,l)=zetamatsubara(k)+eps_diag_complex(l)
                       ginv_full(l,l+1)=cmplx(t,0.0d0,kind=real64)
                       ginv_full(l+1,l)=cmplx(t,0.0d0,kind=real64)
                    enddo
                    ginv_full(1,1)=ginv_full(1,1)+cmplx(t,0.0d0,kind=real64)
                    ginv_full(q,q)=zetamatsubara(k)+eps_diag_complex(q)+cmplx(t,0.0d0,kind=real64)
                    ginv_full(1,q)=cmplx(t,0.0d0,kind=real64)*precalckyexp(1,j)
                    ginv_full(q,1)=cmplx(t,0.0d0,kind=real64)*precalckyexp(2,j)                 
                    call zgetrf(q,q,ginv_full,q,ipiv,info)
                    call zgetri(q,ginv_full,q,ipiv,workinv,10*q,info)
                    do l=1,q
                       do lp=1,q
                          gmatsubarapart(k,l,lp)=gmatsubarapart(k,l,lp)+ginv_full(l,lp)!*cmultipkx(i)*cmultipky(j)
                       enddo
                    enddo
                 endif
                 
              enddo     !End loop over ky (i.e., loop over index j).
           enddo     !End loop over omega/i\nu (i.e., loop over index k).

        enddo     !End loop over sign of kx (i.e., loop over index kxsign)
        !write(*,*)"Retarded Green's function finished."

        !The following block is evaluated if the DOS should be calculated with a low-rank-update.
        if (calcdoslowrankupdate) then
           !Calculate eigenvalues and eigenvectors of the partial (negative tridiagonal) dispersion matrix (eps_subdiag,eps_diag,eps_subdiag) at the given kx-point.
           !Construct subdiagonal part of the partial (negative tridiagonal) dispersion matrix at the given kx-point.
           do l=1,q-1
              eps_subdiag(l)=t
           enddo
           call dstev('V',q,eps_diag,eps_subdiag,eigenvectors,q,work,info)       
           !Precalculate a_j=v_1_j^2+v_q_j^2 and b_j=2*v_1_j*v_q_j for obtaining the weights w_j=a_j+b_j*cos(q*ky) and sum them for equal eigenvalues.
           do l=1,q
              weight_precalc(1,l)=0.0d0
              weight_precalc(2,l)=0.0d0
           enddo
           weight_precalc(1,1)=eigenvectors(1,1)*eigenvectors(1,1)+eigenvectors(q,1)*eigenvectors(q,1)
           weight_precalc(2,1)=2.0d0*eigenvectors(1,1)*eigenvectors(q,1)
           eigenvalues_not_equal(1)=eps_diag(1)
           count_equal=0
           do l=2,q
              if (eps_diag(l)-eps_diag(l-1).le.tolerance) then
                 count_equal=count_equal+1
                 eigenvalues(count_equal)=eps_diag(l)
              endif
              weight_precalc(1,l-count_equal)=weight_precalc(1,l-count_equal)+ &
                   eigenvectors(1,l)*eigenvectors(1,l)+eigenvectors(q,l)*eigenvectors(q,l)
              weight_precalc(2,l-count_equal)=weight_precalc(2,l-count_equal)+ &
                   2.0d0*eigenvectors(1,l)*eigenvectors(q,l)
              eigenvalues_not_equal(l-count_equal)=eps_diag(l)
           enddo     
           do j=1,Nky
              !Calculate weights for the given ky-point.
              count_no_update=0
              do l=1,q-count_equal
                 weight=sqrt(abs(weight_precalc(1,l)+weight_precalc(2,l)*precalcky(j))/2.0d0)
                 if (weight.le.tolerance2) then
                    count_no_update=count_no_update+1
                    eigenvalues(count_no_update+count_equal)=eigenvalues_not_equal(l)
                 else
                    weight_used(l-count_no_update)=weight
                    eigenvalues_used(l-count_no_update)=eigenvalues_not_equal(l)
                 endif
              enddo
              !Use lapack routine to solve secular equation.
              m=q-count_equal-count_no_update
              do l=1,m
                 call dlaed4(m,l,eigenvalues_used(1:m),weight_used(1:m),eigendiff(1:m),2.0d0*t,new_eigenvalue,info)
                 eigenvalues(l+count_equal+count_no_update)=new_eigenvalue
              enddo           
              !Calculate density of states.
              do k=0,Neps
                 sum=0.0d0
                 do l=1,q
                    sum=sum+1.0d0/(deltasquared+(omega(k)+eigenvalues(l))**2)
                 enddo
                 dospartlowrankupdate(k)=dospartlowrankupdate(k)+multipkx(i)*multipky(j)*sum
              enddo
           enddo
        endif

        !The following block is evaluated if the DOS should be calculated without a low rank update.
        if (calcdos) then
           !Construct (negative) full dispersion matrix.
           do j=1,Nky
              do l=1,q
                 do lp=1,q
                    eps_full(l,lp)=cmplx(0.0d0,0.0d0,kind=real64)
                 enddo
              enddo
              do l=1,q-1
                 eps_full(l,l)=eps_diag_complex(l)
                 eps_full(l,l+1)=cmplx(t,0.0d0,kind=real64)
                 eps_full(l+1,l)=cmplx(t,0.0d0,kind=real64)
              enddo
              eps_full(1,1)=eps_full(1,1)+cmplx(t,0.0d0,kind=real64)
              eps_full(q,q)=eps_diag_complex(q)+cmplx(t,0.0d0,kind=real64)
              eps_full(1,q)=cmplx(t,0.0d0,kind=real64)*precalckyexp(1,j)
              eps_full(q,1)=cmplx(t,0.0d0,kind=real64)*precalckyexp(2,j)
              !Calculate eigenvalues of the (negative) dispersion matrix eps_full for the given k-point (kx,ky).
              call zheev('N','U',q,eps_full,q,eigenvalues,lwork,10*q,rwork,info)
              !Calculate density of states.
              do k=0,Neps
                 sum=0.0d0
                 do l=1,q
                    sum=sum+1.0d0/(deltasquared+(omega(k)+eigenvalues(l))**2)
                 enddo
                 dospart(k)=dospart(k)+multipkx(i)*multipky(j)*sum
              enddo
           enddo
        endif
        
     enddo     !End loop over kx (i.e., loop over index i).
     
     !Sum up contributions from different kx to the total retarded/Matsubara Green's function and DOS.
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
     if (calcgretardedlowrankupdate) call MPI_REDUCE(gretardedlocpart(0:Neps,1:q,1:q), &
          gretardedloc(0:Neps,1:q,1:q),q*q*(Neps+1),MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     if (calcgmatsubaralowrankupdate) call MPI_REDUCE(gmatsubaralocpart(0:Neps,1:q,1:q), &
          gmatsubaraloc(0:Neps,1:q,1:q),q*q*(Neps+1),MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     if (calcgretarded) call MPI_REDUCE(gretardedpart(0:Neps,1:q,1:q),gretarded(0:Neps,1:q,1:q),q*q*(Neps+1), &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     if (calcgmatsubara) call MPI_REDUCE(gmatsubarapart(0:Neps,1:q,1:q),gmatsubara(0:Neps,1:q,1:q),q*q*(Neps+1), &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     if (calcdoslowrankupdate) call MPI_REDUCE(dospartlowrankupdate(0:Neps),doslowrankupdate(0:Neps),Neps+1, &
          MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     if (calcdos) call MPI_REDUCE(dospart(0:Neps),dos(0:Neps),Neps+1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     
     !Normalize retarded Green's function and DOS.
     do k=0,Neps
        do l=1,q
           do lp=1,q
              if (calcgretardedlowrankupdate) gretardedloc(k,l,lp)=gretardedloc(k,l,lp)/(cmplx(4.0d0*normkxsum*normkysum,0.0d0,kind=real64))
              if (calcgmatsubaralowrankupdate) gmatsubaraloc(k,l,lp)=gmatsubaraloc(k,l,lp)/(cmplx(4.0d0*normkxsum*normkysum,0.0d0,kind=real64))
              if (calcgretarded) gretarded(k,l,lp)=gretarded(k,l,lp)/(cmplx(4.0d0*normkxsum*normkysum,0.0d0,kind=real64))
              if (calcgmatsubara) gmatsubara(k,l,lp)=gmatsubara(k,l,lp)/(cmplx(4.0d0*normkxsum*normkysum,0.0d0,kind=real64))
           enddo
        enddo
        if (calcdoslowrankupdate) doslowrankupdate(k)=delta*doslowrankupdate(k)/(pi*normkxsum*normkysum)        
        if (calcdos) dos(k)=delta*dos(k)/(pi*normkxsum*normkysum)
     enddo
     if (myid.eq.0) then
        !Write retarded/Matsubara Green's function to file.
        write(13,*)"#Retarded Green's function for parameter set number",n,"start."
        write(14,*)"#Matsubara Green's function for parameter set number",n,"start."
        write(15,*)"#Retarded Green's function for parameter set number",n,"start."
        write(16,*)"#Matsubara Green's function for parameter set number",n,"start."                
        do l=1,q
           do lp=1,q
              do i=0,Neps
                 if (calcgretardedlowrankupdate) write(13,'(2I6,3f30.12)')l,lp,omega(i), real(gretardedloc(i,l,lp),kind=real64), &
                      aimag(gretardedloc(i,l,lp))
                 if (calcgmatsubaralowrankupdate) write(14,'(2I6,3f30.12)')l,lp,nu(i), real(gmatsubaraloc(i,l,lp),kind=real64), &
                      aimag(gmatsubaraloc(i,l,lp))
                 if (calcgretarded) write(15,'(2I6,3f30.12)')l,lp,omega(i), real(gretarded(i,l,lp),kind=real64), &
                      aimag(gretarded(i,l,lp))
                 if (calcgmatsubara) write(16,'(2I6,3f30.12)')l,lp,nu(i), real(gmatsubara(i,l,lp),kind=real64), &
                      aimag(gmatsubara(i,l,lp))
              enddo
           enddo
        enddo
        write(13,*)"#Retarded Green's function for parameter set number",n,"end."
        write(14,*)"#Matsubara Green's function for parameter set number",n,"end."
        write(15,*)"#Retarded Green's function for parameter set number",n,"end."
        write(16,*)"#Matsubara Green's function for parameter set number",n,"end."
        !Write density of states to file.
        write(11,*)"#DOS for parameter set number",n,"start."
        write(12,*)"#DOS for parameter set number",n,"start."
        do i=0,Neps
           if (calcdoslowrankupdate) write(11,'(2f30.12)')omega(i),doslowrankupdate(i)
           if (calcdos) write(12,'(2f30.12)')omega(i),dos(i)
        enddo
        write(11,*)"#DOS for parameter set number",n,"end."
        write(12,*)"#DOS for parameter set number",n,"end."
     endif
     
     !Deallocate arrays.
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)  
     deallocate(Bmagshift)
     deallocate(eps_diag)
     deallocate(eps_diag_complex)
     deallocate(eps_subdiag)
     deallocate(eps_subdiag_complex)
     deallocate(eps_superdiag_complex)
     deallocate(eps_full)
     deallocate(ginvretarded_diag)
     deallocate(ginvmatsubara_diag)
     deallocate(gretardedtridiag)
     deallocate(gmatsubaratridiag)
     deallocate(gtridiagv)
     deallocate(vdaggtridiag)
     deallocate(ginv_full)
     deallocate(eigenvectors)
     deallocate(eigenvalues_not_equal)
     deallocate(weight_precalc)
     deallocate(precalcky)
     deallocate(precalcky2)
     deallocate(precalckyexp)
     deallocate(eigenvalues_used)
     deallocate(weight_used)
     deallocate(eigenvalues)
     deallocate(omega)
     deallocate(nu)
     deallocate(zetaretarded)
     deallocate(zetamatsubara)
     deallocate(multipkx)
     deallocate(multipky)
     deallocate(ipiv)
     deallocate(work)
     deallocate(lwork)
     deallocate(rwork)
     deallocate(workinv)
     deallocate(eigendiff)
     if (calcgretardedlowrankupdate) deallocate(gretardedlocpart)
     if (calcgretardedlowrankupdate) deallocate(gretardedloc)
     if (calcgretarded) deallocate(gretardedpart)
     if (calcgretarded) deallocate(gretarded)
     if (calcgmatsubaralowrankupdate) deallocate(gmatsubaralocpart)
     if (calcgmatsubaralowrankupdate) deallocate(gmatsubaraloc)
     if (calcgmatsubara) deallocate(gmatsubarapart)
     if (calcgmatsubara) deallocate(gmatsubara)
     if (calcdoslowrankupdate) deallocate(dospartlowrankupdate)
     if (calcdoslowrankupdate) deallocate(doslowrankupdate)
     if (calcdos) deallocate(dospart)
     if (calcdos) deallocate(dos)
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
     
  !End loop over all input sets.
  enddo

  !Close open files (root process only).
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  if (myid.eq.0) then
     close(10)
     close(11)
     close(12)
     close(13)
     close(14)
     close(15)
     close(16)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
  call MPI_FINALIZE(ierror)
  
  stop
  
contains
    
  subroutine construct_diag_dispersion(q,Bmagshift,t,kx,eps_diag)
    !Input
    integer, intent(in) :: q
    real(real64), intent(in) :: kx,t
    real(real64), dimension(:), intent(in) :: Bmagshift
    !Output
    real(real64), dimension(:), intent(out) :: eps_diag
    !Subroutine internal variables
    integer :: l
    !Construct diagonal part of dispersion matrix.
    do l=1,q
       eps_diag(l)=2.0d0*t*cos(kx+Bmagshift(l))
    enddo
    !Add correction term to first and last element for low rank updates.
    eps_diag(1)=eps_diag(1)-t
    eps_diag(q)=eps_diag(q)-t
  end subroutine construct_diag_dispersion

end program DensityOfStates
