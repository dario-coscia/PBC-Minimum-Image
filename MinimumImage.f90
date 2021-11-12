! ============ ============= ============= ============= =============
! FORTRAN90 program to calculate minimum distance, using minimum image 
! convention for 2D (1,1,1) lattice. The ideces (m,n,l) refers to the
! Miller indeces used in crystallography. The algorithm projects the 
! Euclidean distance between the two particles into the two basis vectors, 
! then classical minimum image is applied. 
!
! To compile and run the code simply type:
!
!       gfortran MinimumImage.f90 -o MinimumImage.x
!       ./MinimumImage.x < try.txt
!
! If everything is correct you should get as output:
!
!       Minimum distance between the two particles:    0.00000000
!
! ============ ============= ============= ============= =============

PROGRAM Minimum_Image 
    IMPLICIT NONE
    
    INTEGER,PARAMETER :: rk=selected_real_kind(6)
    REAL(kind=rk) ::x1, x2, y1, y2, sep, alpha, beta
    REAL(kind=rk),DIMENSION(2) :: dr, dr1, dr2, L1, L2, mod_1, mod_2
    
    !----INPUT VALUES----!
    
    ! (x1,y1) position particle 1 in cartesian coordinates
    ! (x2,y2) position particle 2 in cartesian coordinates
    ! L1, L2 two basis supercell vector for (m,l,l) = (111) cell
    
    READ*,  x1, y1, x2, y2, L1, L2      ! Reading input values from terminal
    
    !--------------------!
    
    dr(1) = x1 - x2
    dr(2) = y1 - y2
    
    beta = (L1(1)*dr(2)-L1(2)*dr(1))/(L1(1)*L2(2)-L1(2)*L2(1))
    alpha = (L2(2)*dr(1)-L2(1)*dr(2))/(L1(1)*L2(2)-L1(2)*L2(1))

    dr1 = alpha*L1
    dr2 = beta*L2 
    
    mod_1 = L1*anint(norm2(dr1)/norm2(L1))
    mod_2 = L2*anint(norm2(dr2)/norm2(L2))
    
    IF ( norm2(dr1 - mod_1)<= norm2(dr1 + mod_1)) THEN
        dr1 = dr1 - mod_1
    ELSE
        dr1 = dr1 + mod_1
    END IF
        
    IF ( norm2(dr2 - mod_2)<= norm2(dr2 + mod_2)) THEN
        dr2 = dr2 - mod_2
    ELSE
        dr2 = dr2 + mod_2
    END IF
    
    sep = norm2(dr1 + dr2)             ! Minimum distance between particle 1 and particle 2
    
    PRINT*, "Minimum distance between the two particles:",sep

END PROGRAM