MODULE grids
IMPLICIT NONE
CONTAINS
    FUNCTION must_finer(point, boxsize, min_spacing, max_spacing, alpha1, &
                        alpha2, threshold, coordinates, n)
!f2py   intent(hide) n
        REAL(8),DIMENSION(3) :: point
        REAL(8) :: boxsize, min_spacing, max_spacing, alpha1, alpha2, threshold
        REAL(8),DIMENSION(n,3) :: coordinates
        INTEGER :: n
        LOGICAL :: must_finer
        ! local vars
        INTEGER :: i
        REAL(8) :: distance, new_distance, s
        IF (boxsize .gt. max_spacing) THEN
            must_finer = .true.
        ELSE IF (boxsize .lt. min_spacing) THEN
            must_finer = .false.
        ELSE
            distance = SQRT(                          &
                (coordinates(1,1) - point(1))**2 +  &
                (coordinates(1,2) - point(2))**2 +  &
                (coordinates(1,3) - point(3))**2    &
            )
            DO i=2,n
                new_distance = SQRT(            &
                    (coordinates(i,1) - point(1))**2 +   &
                    (coordinates(i,2) - point(2))**2 +   &
                    (coordinates(i,3) - point(3))**2     &
                )
                IF (new_distance .lt. distance) distance = new_distance
            END DO
            !must_finer = (distance*alpha .lt. boxsize)
            !must_finer = (alpha*(distance - 0.3*exp(-(2*(distance-1))**2)) .lt. boxsize)
            !must_finer = ( &
            !    ((alpha1*distance*threshold/distance + alpha2*distance**5) / &
            !     (threshold/distance + distance**4)) &
            !    .lt. boxsize &
            !)
            s = (1 - EXP(-distance/threshold))**2
            must_finer = (alpha1*(1-s) + alpha2*s)*distance .lt. boxsize
        END IF
    END FUNCTION

    SUBROUTINE recursive_grid(center, boxsize, div, min_spacing, max_spacing, &
                              alpha1, alpha2, threshold, coordinates, n, cb)
!f2py   intent(hide) n

        REAL(8),DIMENSION(3) :: center
        REAL(8) :: boxsize
        INTEGER :: div
        REAL(8) :: min_spacing, max_spacing, alpha1, alpha2, threshold
        REAL(8),DIMENSION(n,3) :: coordinates
        INTEGER :: n
        EXTERNAL :: cb


        REAL(8) :: subsize
        REAL(8),DIMENSION(3) :: point
        INTEGER :: i,j,k
        subsize = boxsize/div
        DO i=1,div
            DO j=1,div
                DO k=1,div
                    point(1) = (i-0.5-0.5*div)*subsize + center(1)
                    point(2) = (j-0.5-0.5*div)*subsize + center(2)
                    point(3) = (k-0.5-0.5*div)*subsize + center(3)
                    IF (must_finer(point, subsize, min_spacing, max_spacing, &
                        alpha1, alpha2, threshold, coordinates, n)) THEN
                        CALL recursive_grid(point, subsize, div, min_spacing, &
                                            max_spacing, alpha1, alpha2,      &
                                            threshold, coordinates, n, cb)
                    ELSE
                        CALL cb(point, subsize**3)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE pseudo_log_grid(low, sizes, delta, beta, numu, coordinates, n, cb)
!f2py   intent(hide) n

        REAL(8),DIMENSION(3) :: low
        INTEGER,DIMENSION(3) :: sizes
        REAL(8) :: delta, beta
        REAL(8),DIMENSION(n,3) :: coordinates
        INTEGER :: n, numu
        EXTERNAL :: cb


        REAL(8) :: distance,new_distance,dp,dm,uxmin,uxmax,uymin,uymax,uzmin, &
                   uzmax,wx,wy,wz,ux,uy,uz,w,dux,duy,duz
        REAL(8),DIMENSION(3) :: point, gridp

        INTEGER :: i,j,k,l,closest,nn,in,jn,kn

        PRINT *, "low", low
        PRINT *, "sizes", sizes
        PRINT *, "delta", delta
        PRINT *, "coordinates", coordinates

        DO i=1,sizes(1)
            DO j=1,sizes(2)
                DO k=1,sizes(3)
                    !PRINT *, "indices", i, j, k
                    point(1) = (i-0.5)*delta + low(1)
                    point(2) = (j-0.5)*delta + low(2)
                    point(3) = (k-0.5)*delta + low(3)
                    !PRINT *, "point", point

                    closest = 1
                    distance = SQRT(                          &
                        (coordinates(1,1) - point(1))**2 +  &
                        (coordinates(1,2) - point(2))**2 +  &
                        (coordinates(1,3) - point(3))**2    &
                    )
                    DO l=2,n
                        new_distance = SQRT(            &
                            (coordinates(l,1) - point(1))**2 +   &
                            (coordinates(l,2) - point(2))**2 +   &
                            (coordinates(l,3) - point(3))**2     &
                        )
                        IF (new_distance .lt. distance) THEN
                            distance = new_distance
                            closest = l
                        END IF
                    END DO

                    !PRINT *, "distance=", distance

                    dp = (distance+delta*1.73/2)
                    dm = (distance-delta*1.73/2)

                    !PRINT *, "dp=", dp, "dm=", dm

                    ! numu = aantal gridpunten op u-as
                    dp = numu/beta*asinh(dp*sinh(beta)/8D0)
                    dm = numu/beta*asinh(dm*sinh(beta)/8D0)

                    !PRINT *, "dp=", dp, "dm=", dm

                    nn = CEILING(abs(dp-dm))
                    !PRINT *, "nn", nn

                    uxmin = asinh(                                     &
                        (point(1) - delta/2 - coordinates(closest,1))* &
                        sinh(beta)/8D0                                  &
                    )/beta
                    uxmax = asinh(                                     &
                        (point(1) + delta/2 - coordinates(closest,1))* &
                        sinh(beta)/8D0                                  &
                    )/beta

                    uymin = asinh(                                     &
                        (point(2) - delta/2 - coordinates(closest,2))* &
                        sinh(beta)/8D0                                  &
                    )/beta
                    uymax = asinh(                                     &
                        (point(2) + delta/2 - coordinates(closest,2))* &
                        sinh(beta)/8D0                                  &
                    )/beta

                    uzmin = asinh(                                     &
                        (point(3) - delta/2 - coordinates(closest,3))* &
                        sinh(beta)/8D0                                  &
                    )/beta
                    uzmax = asinh(                                     &
                        (point(3) + delta/2 - coordinates(closest,3))* &
                        sinh(beta)/8D0                                  &
                    )/beta
                    !PRINT *, i, uzmax, uzmin

                    DO in=1,nn
                        dux = (uxmax - uxmin)/nn
                        ux = uxmin + (in - 0.5)*dux
                        gridp(1) = coordinates(closest,1) + &
                                   8D0*sinh(beta*ux)/sinh(beta)
                        !wx = (uxmax - uxmin)/nn*(beta*8.0/sinh(beta)*cosh(beta*ux))
                        wx = 8D0*sinh(beta*(ux+0.5*dux))/sinh(beta) - &
                             8D0*sinh(beta*(ux-0.5*dux))/sinh(beta)
                        DO jn=1,nn
                            duy = (uymax - uymin)/nn
                            uy = uymin + (jn - 0.5)*duy
                            gridp(2) = coordinates(closest,2) + &
                                       8D0*sinh(beta*uy)/sinh(beta)
                            !wy = (uymax - uymin)/nn*(beta*8.0/sinh(beta)*cosh(beta*uy))
                            wy = 8D0*sinh(beta*(uy+0.5*duy))/sinh(beta) - &
                                 8D0*sinh(beta*(uy-0.5*duy))/sinh(beta)
                            DO kn=1,nn
                                duz = (uzmax - uzmin)/nn
                                uz = uzmin + (kn - 0.5)*duz
                                gridp(3) = coordinates(closest,3) + &
                                           8D0*sinh(beta*uz)/sinh(beta)
                                !wz = (uzmax - uzmin)/nn*(beta*8.0/sinh(beta)*cosh(beta*uz))
                                wz = 8D0*sinh(beta*(uz+0.5*duz))/sinh(beta) - &
                                     8D0*sinh(beta*(uz-0.5*duz))/sinh(beta)
                                !PRINT *, gridp, wx*wy*wz
                                w = wx*wy*wz
                                CALL cb(gridp, w)
                            END DO
                        END DO
                    END DO

                END DO
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE list_neighbors(grid, volumes, callback, ngrid)
        INTEGER :: ngrid
        REAL(8),DIMENSION(ngrid,3),INTENT(IN) :: grid
        REAL(8),DIMENSION(ngrid),INTENT(IN) :: volumes
        EXTERNAL callback
!f2py   intent(hide) :: ngrid

        INTEGER :: i,j
        REAL(8),DIMENSION(3) :: delta
        REAL(8) :: distance
        REAL(8),DIMENSION(ngrid) :: sizes

        sizes = volumes**(1.0/3.0)
        DO i=1,ngrid
            DO j=1,i-1
                delta = grid(i,:) - grid(j,:)
                distance = MAXVAL(delta)
                IF (distance < sizes(i) + sizes(j) + 0.5*MIN(sizes(i),sizes(j))) THEN
                    CALL callback(i,j)
                END IF
            END DO
        END DO
    END SUBROUTINE
END MODULE

MODULE potential
IMPLICIT NONE
    REAL(8),PARAMETER :: pi = 3.141592653589793238462643d0
CONTAINS
    FUNCTION gridv(observer, grid, densities, volumes, ion_charges, ion_coordinates, ngrid, nion)
        INTEGER,INTENT(IN) :: ngrid,nion
        REAL(8),DIMENSION(3),INTENT(IN) :: observer
        REAL(8),DIMENSION(ngrid,3),INTENT(IN) :: grid
        REAL(8),DIMENSION(ngrid),INTENT(IN) :: densities, volumes
        REAL(8),DIMENSION(nion),INTENT(IN) :: ion_charges
        REAL(8),DIMENSION(nion,3),INTENT(IN) :: ion_coordinates
        REAL(8) :: gridv
!f2py   intent(hide) :: ndata,nions
        INTEGER :: i
        REAL(8) :: distance
        gridv = 0
        DO i=1,ngrid
            distance = SQRT(SUM( (observer - grid(i,:))**2 ))
            gridv = gridv - volumes(i)*densities(i)/distance
        END DO
        DO i=1,nion
            distance = SQRT(SUM( (observer - ion_coordinates(i,:))**2 ))
            gridv = gridv + ion_charges(i)/distance
        END DO
     END FUNCTION

     FUNCTION switch_cos(distance, low, high)
        REAL(8),INTENT(IN) :: distance, low, high
        REAL(8) :: switch_cos
        IF (distance < low) THEN
            !PRINT *, "ERONDER", distance, " < ", low
            switch_cos = 0D0
        ELSE IF (distance > high) THEN
            !PRINT *, "EROP", distance, " > ", high
            switch_cos = 1D0
        ELSE
            !PRINT *, "ERTUSSEN", low, " < ", distance, " < ", high
            switch_cos = 0.5*(1-COS(pi*(distance-low)/(high-low)))
        END IF
    END FUNCTION

    FUNCTION measure_ortho(point, observer, direction, length)
        REAL(8),DIMENSION(3),INTENT(IN) :: point,observer,direction
        REAL(8),INTENT(IN) :: length
        REAL(8) :: measure_ortho
        measure_ortho = SUM( (point-observer)*direction ) / length
    END FUNCTION

    SUBROUTINE spherical_switchv(observer,rho,vols,ionval,ioncor, &
                                 switch_low,switch_high,closest,V, &
                                 qedep,qidep,ionic_s,ndata,nions)
        INTEGER,INTENT(IN) :: ndata,nions
        REAL(8),DIMENSION(3),INTENT(IN) :: observer
        REAL(8),DIMENSION(ndata,4),INTENT(IN) :: rho
        REAL(8),DIMENSION(ndata),INTENT(IN) :: vols
        REAL(8),DIMENSION(nions),INTENT(IN) :: ionval
        REAL(8),DIMENSION(nions,3),INTENT(IN) :: ioncor
        REAL(8),INTENT(IN) :: switch_low,switch_high
        INTEGER,DIMENSION(ndata),INTENT(IN) :: closest
        REAL(8),INTENT(OUT) :: V,qedep,qidep
        REAL(8),DIMENSION(nions),INTENT(OUT) :: ionic_s
!f2py   intent(hide) :: ndata,nions
        INTEGER :: i,n
        !INTEGER,DIMENSION(1) :: m
        REAL(8) :: d, s
        REAL(8),DIMENSION(3) :: point
        REAL(8),DIMENSION(nions) :: ionic_w
        V = 0
        qedep = 0
        qidep = 0
        ionic_s = 0
        ionic_w = 0
        DO i=1,ndata
            point = rho(i,1:3)
            d = SQRT(SUM( (observer - point)**2 ))
            s = switch_cos(d, switch_low, switch_high)
            IF (s > 0) THEN
                V = V - s*vols(i)*rho(i,4)/d
                n = closest(i)+1
                ionic_s(n) = ionic_s(n) + s*vols(i)*rho(i,4)!/d
                ionic_w(n) = ionic_w(n) + vols(i)*rho(i,4)!/d
            END IF
            qedep = qedep - (1-s)*vols(i)*rho(i,4)
        END DO
        WHERE (ionic_w > 0)
            ionic_s = ionic_s / ionic_w
        ELSEWHERE
            ionic_s = 0.0
        END WHERE
        DO i=1,nions
            point = ioncor(i,:)
            d = SQRT(SUM( (observer - point)**2 ))
            V = V + ionic_s(i)*ionval(i)/d
            qidep = qidep + (1-ionic_s(i))*ionval(i)
        END DO
    END SUBROUTINE

    SUBROUTINE linear_switchv(observer,direction,length,rho,vols,ionval,ioncor,switch_low,switch_high,V,qedep,qidep,ndata,nions)
        INTEGER,INTENT(IN) :: ndata,nions
        REAL(8),DIMENSION(3),INTENT(IN) :: observer,direction
        REAL(8),INTENT(IN) :: length,switch_low,switch_high
        REAL(8),DIMENSION(ndata,4),INTENT(IN) :: rho
        REAL(8),DIMENSION(ndata),INTENT(IN) :: vols
        REAL(8),DIMENSION(nions),INTENT(IN) :: ionval
        REAL(8),DIMENSION(nions,3),INTENT(IN) :: ioncor
        REAL(8),INTENT(OUT) :: V,qedep,qidep
!f2py   intent(hide) :: ndata,nions
        INTEGER :: i
        REAL(8) :: d, s!, ssum, counter
        REAL(8),DIMENSION(3) :: point
        !PRINT *, "OBSERVER", observer
        !PRINT *, "DIRECTION", direction
        !PRINT *, "LENGTH", length
        V = 0
        qedep = 0
        qidep = 0
        !ssum = 0
        !counter = 0
        DO i=1,ndata
            point = rho(i,1:3)
            d = SQRT(SUM( (observer - point)**2 ))
            s = measure_ortho(point, observer, direction, length)
            s = switch_cos(s, switch_low, switch_high)
            !PRINT *, "S", s
            !ssum = ssum + s
            !counter = counter + 1
            IF (s > 0) THEN
                V = V - s*vols(i)*rho(i,4)/d
            END IF
            qedep = qedep - (1-s)*vols(i)*rho(i,4)
        END DO
        DO i=1,nions
            point = ioncor(i,:)
            d = SQRT(SUM( (observer - point)**2 ))
            s = measure_ortho(point, observer, direction, length)
            s = switch_cos(s, switch_low, switch_high)
            !PRINT *, "S", s
            !ssum = ssum + s
            !counter = counter + 1
            IF (s > 0) THEN
                V = V + s*ionval(i)/d
            END IF
            qidep = qidep + (1-s)*ionval(i)
        END DO
        !PRINT *, "SSUM", ssum
        !PRINT *, "COUNTER - SSUM", (counter - ssum)
        !PRINT *, "COUNTER", counter
    END SUBROUTINE

    SUBROUTINE all_distances(ioncor, grid, nion, ngrid, distances)
        INTEGER,INTENT(IN) :: nion, ngrid
        REAL(8),DIMENSION(nion,3),INTENT(IN) :: ioncor
        REAL(8),DIMENSION(ngrid,4),INTENT(IN) :: grid
        REAL(8),DIMENSION(ngrid,nion),INTENT(OUT) :: distances
!f2py   intent(hide) :: nion,ngrid
        INTEGER :: i,j
        DO j=1,ngrid
            DO i=1,nion
                print *,i,j
                distances(i,j) = SQRT(SUM( (ioncor(i,:) - grid(j,1:3))**2 ))
            END DO
        END DO
    END SUBROUTINE


    FUNCTION neg_switch_cos(distance, low, high)
        REAL(8),INTENT(IN) :: distance, low, high
        REAL(8) :: neg_switch_cos
        IF (distance < low) THEN
            !PRINT *, "ERONDER", distance, " < ", low
            neg_switch_cos = 1D0
        ELSE IF (distance > high) THEN
            !PRINT *, "EROP", distance, " > ", high
            neg_switch_cos = 0D0
        ELSE
            !PRINT *, "ERTUSSEN", low, " < ", distance, " < ", high
            neg_switch_cos = 0.5*(1+COS(pi*(distance-low)/(high-low)))
        END IF
    END FUNCTION

    SUBROUTINE neg_spherical_switchv(observer,rho,vols,ionval,ioncor, &
                                     switch_low,switch_high,closest,V, &
                                     qedep,qidep,ionic_s,ndata,nions)
        INTEGER,INTENT(IN) :: ndata,nions
        REAL(8),DIMENSION(3),INTENT(IN) :: observer
        REAL(8),DIMENSION(ndata,4),INTENT(IN) :: rho
        REAL(8),DIMENSION(ndata),INTENT(IN) :: vols
        REAL(8),DIMENSION(nions),INTENT(IN) :: ionval
        REAL(8),DIMENSION(nions,3),INTENT(IN) :: ioncor
        REAL(8),INTENT(IN) :: switch_low,switch_high
        INTEGER,DIMENSION(ndata),INTENT(IN) :: closest
        REAL(8),INTENT(OUT) :: V,qedep,qidep
        REAL(8),DIMENSION(nions),INTENT(OUT) :: ionic_s
!f2py   intent(hide) :: ndata,nions
        INTEGER :: i,n
        !INTEGER,DIMENSION(1) :: m
        REAL(8) :: d, s
        REAL(8),DIMENSION(3) :: point
        REAL(8),DIMENSION(nions) :: ionic_w
        V = 0
        qedep = 0
        qidep = 0
        ionic_s = 0
        ionic_w = 0
        DO i=1,ndata
            point = rho(i,1:3)
            d = SQRT(SUM( (observer - point)**2 ))
            s = neg_switch_cos(d, switch_low, switch_high)
            IF (s > 0) THEN
                V = V - s*vols(i)*rho(i,4)/d
                n = closest(i)+1
                ionic_s(n) = ionic_s(n) + s*vols(i)*rho(i,4)!/d
                ionic_w(n) = ionic_w(n) + vols(i)*rho(i,4)!/d
            END IF
            qedep = qedep - (1-s)*vols(i)*rho(i,4)
        END DO
        WHERE (ionic_w > 0)
            ionic_s = ionic_s / ionic_w
        ELSEWHERE
            ionic_s = 0.0
        END WHERE
        DO i=1,nions
            point = ioncor(i,:)
            d = SQRT(SUM( (observer - point)**2 ))
            V = V + ionic_s(i)*ionval(i)/d
            qidep = qidep + (1-ionic_s(i))*ionval(i)
        END DO
    END SUBROUTINE

END MODULE

MODULE partitioning
IMPLICIT NONE
CONTAINS
    FUNCTION lin_interpol(x, xs, ys, n, lefty, righty)
        INTEGER,INTENT(IN) :: n
        REAL(8),INTENT(IN) :: x, lefty, righty
        REAL(8),INTENT(IN),DIMENSION(n) :: xs, ys
        REAL(8) :: lin_interpol

        INTEGER :: b,e,j
        IF (x <= xs(1)) THEN
            lin_interpol = lefty
        ELSE IF (x >= xs(n)) THEN
            lin_interpol = righty
        ELSE
            b = 1
            e = n
            DO
                j = (b+e)/2
                IF (xs(j) > x) THEN
                    e = j
                ELSEIF (xs(j+1) < x) THEN
                    b = j+1
                ELSE
                    b = j
                    e = j+1
                    EXIT
                END IF
            END DO
            lin_interpol = (x - xs(b))/(xs(e) - xs(b))*(ys(e) - ys(b)) + ys(b)
        END IF
    END FUNCTION

    ! Calculate hirshfeld weights
    SUBROUTINE hirshfeld_weights(grid, weights, atom_densities, atom_indices, &
                                 atom_coordinates, ngrid, nrad, natden, natom)
        INTEGER,INTENT(IN) :: ngrid, nrad, natden, natom
        REAL(8),INTENT(IN),DIMENSION(ngrid,3) :: grid
        REAL(8),INTENT(IN),DIMENSION(natden,2,nrad) :: atom_densities
        INTEGER,INTENT(IN),DIMENSION(natom) :: atom_indices
        REAL(8),INTENT(IN),DIMENSION(natom,3) :: atom_coordinates
        REAL(8),INTENT(OUT),DIMENSION(ngrid,natom) :: weights
!cf2py  intent(hide) :: ngrid, nrad, natden, natom
        INTEGER :: i,j
        REAL(8) :: distance, total_weight
        DO i=1,ngrid
            DO j=1,natom
                distance = SQRT(SUM( (grid(i,:) - atom_coordinates(j,:))**2 ))
                weights(i,j) = lin_interpol(               &
                    distance,                              &
                    atom_densities(atom_indices(j)+1,1,:),   &
                    atom_densities(atom_indices(j)+1,2,:),   &
                    nrad,                                  &
                    atom_densities(atom_indices(j)+1,2,1),   &
                    0D0                                    &
                )
            END DO
            total_weight = SUM(weights(i,:))
            IF (total_weight > 0) THEN
                weights(i,:) = weights(i,:) / total_weight
            END IF
        END DO
    END SUBROUTINE

    ! Calculate voronoi weights
    SUBROUTINE voronoi_weights(grid, weights, atom_coordinates, ngrid, natom)
        INTEGER,INTENT(IN) :: ngrid, natom
        REAL(8),INTENT(IN),DIMENSION(ngrid,3) :: grid
        REAL(8),INTENT(IN),DIMENSION(natom,3) :: atom_coordinates
        REAL(8),INTENT(OUT),DIMENSION(ngrid,natom) :: weights
!cf2py  intent(hide) :: ngrid, natom
        INTEGER :: i,j
        INTEGER,DIMENSION(1) :: tmp
        DO i=1,ngrid
            DO j=1,natom
                weights(i,j) = SQRT(SUM( (grid(i,:) - atom_coordinates(j,:))**2 ))
            END DO
            tmp = MINLOC(weights(i,:))
            weights(i,:) = 0.0
            weights(i,tmp(1)) = 1.0
        END DO
    END SUBROUTINE

    ! Calculate charges and dipoles based on densities and weights
    SUBROUTINE charges_and_dipoles(grid, volumes, densities, weights,  &
                                   atom_coordinates, cd, ngrid, natom)
        INTEGER,INTENT(IN) :: ngrid, natom
        REAL(8),INTENT(IN),DIMENSION(ngrid,3) :: grid
        REAL(8),INTENT(IN),DIMENSION(ngrid) :: volumes, densities
        REAL(8),INTENT(IN),DIMENSION(ngrid,natom) :: weights
        REAL(8),INTENT(IN),DIMENSION(natom,3) :: atom_coordinates
        REAL(8),INTENT(OUT),DIMENSION(natom,4) :: cd
!cf2py  intent(hide) :: ngrid, natom
        INTEGER :: i,j
        REAL(8) :: c
        DO i=1,ngrid
            c = volumes(i)*densities(i)
            cd(:,1) = cd(:,1) + c*weights(i,:)
            DO j=1,natom
                cd(j,2:4) = (                                            &
                    cd(j,2:4) +                                          &
                    c*weights(i,j)*(grid(i,:) - atom_coordinates(j,:))   &
                )
            END DO
        END DO
    END SUBROUTINE

    FUNCTION mfn(x,y,z)
        REAL(8),INTENT(IN) :: x,y,z
        REAL(8) :: r2
        REAL(8),DIMENSION(25) :: mfn
        REAL(8),PARAMETER :: c20  = 1.0/2.0
        REAL(8),PARAMETER :: c21  = SQRT(3.0)
        REAL(8),PARAMETER :: c22c = SQRT(3.0/4.0)
        REAL(8),PARAMETER :: c31  = SQRT(3.0/8.0)
        REAL(8),PARAMETER :: c32c = SQRT(15.0/4.0)
        REAL(8),PARAMETER :: c32s = SQRT(15.0)
        REAL(8),PARAMETER :: c33  = SQRT(5.0/8.0)
        REAL(8),PARAMETER :: c40  = 1.0/8.0
        REAL(8),PARAMETER :: c42c = SQRT(5.0/16.0)
        REAL(8),PARAMETER :: c42s = SQRT(5.0/4.0)
        REAL(8),PARAMETER :: c43  = SQRT(35.0/8.0)
        REAL(8),PARAMETER :: c44c = SQRT(35.0/64.0)
        REAL(8),PARAMETER :: c44s = SQRT(35.0/4.0)

        r2 = x*x+y*y+z*z
        ! monopole, 00
        mfn(1) = 1
        ! dipole, 10, 11c, 11s
        mfn(2) = z
        mfn(3) = x
        mfn(4) = y
        ! quadrupole, 20, 21c, 21s, 22c, 22s
        mfn(5) = c20 *(3*z*z-r2)
        mfn(6) = c21 *x*z
        mfn(7) = c21 *y*z
        mfn(8) = c22c*(x*x-y*y)
        mfn(9) = c21 *x*y
        ! octopole, 30, 31C, 31s, 32c, 32s, 33c, 33s
        mfn(10) = c20 *(5*z**3-3*z*r2)
        mfn(11) = c31 *x*(5*z*z-r2)
        mfn(12) = c31 *y*(5*z*z-r2)
        mfn(13) = c32c*z*(x*x-y*y)
        mfn(14) = c32s*x*y*z
        mfn(15) = c33 *x*(x*x-3*y*y)
        mfn(16) = c33 *y*(3*x*x-y*y)
        ! hexadecapole, 40, 41c, 41s, 42c, 42s, 43c, 43s, 44c, 44s
        mfn(17) = c40 *(35*z**4-30*z*z*r2+3*r2*r2)
        mfn(18) = c33 *x*z*(7*z*z-3*r2)
        mfn(19) = c33 *y*z*(7*z*z-3*r2)
        mfn(20) = c42c*(x*x-y*y)*(7*z*z-r2)
        mfn(21) = c42s*x*y*(7*z*z-r2)
        mfn(22) = c43 *z*x*(x*x-3*y*y)
        mfn(23) = c43 *z*y*(3*x*x-y*y)
        mfn(24) = c44c*(x**4-6*x*x*y*y+y**4)
        mfn(25) = c44s*x*y*(x*x-y*y)
    END FUNCTION

    SUBROUTINE multipoles(grid, volumes, densities, weights,  &
                          atom_coordinates, output, ngrid, natom)
        INTEGER,INTENT(IN) :: ngrid, natom
        REAL(8),INTENT(IN),DIMENSION(ngrid,3) :: grid
        REAL(8),INTENT(IN),DIMENSION(ngrid) :: volumes, densities
        REAL(8),INTENT(IN),DIMENSION(ngrid,natom) :: weights
        REAL(8),INTENT(IN),DIMENSION(natom,3) :: atom_coordinates
        REAL(8),INTENT(OUT),DIMENSION(natom,25) :: output
!cf2py  intent(hide) :: ngrid, natom
        INTEGER :: i,j
        REAL(8) :: c
        REAL(8),DIMENSION(3) :: delta
        output = 0
        DO i=1,ngrid
            DO j=1,natom
                c = volumes(i)*densities(i)*weights(i,j)
                delta = grid(i,:) - atom_coordinates(j,:)
                output(j,:) = output(j,:) + c*mfn(delta(1),delta(2),delta(3))
            END DO
        END DO
    END SUBROUTINE
END MODULE
