MODULE fsim
IMPLICIT NONE
CONTAINS
    SUBROUTINE similarity_rows(row1, row2, margin, reciproke, result, n1, n2)
!f2py   intent(hide) n1, n2
        REAL(8),DIMENSION(n1) :: row1
        REAL(8),DIMENSION(n2) :: row2
        REAL(8) :: margin, reciproke
        REAL(8),INTENT(out) :: result
        INTEGER :: n1, n2
        ! local vars
        INTEGER :: i1, i2
        REAL(8) :: delta, mean

        result = 0
        DO i1=1,n1
            DO i2=1,n2
                mean = 0.5*(row1(i1)+row2(i2))
                IF (1/mean .GT. reciproke) THEN
                    delta = row1(i1)-row2(i2)
                    result = result + EXP(-(delta/(mean*margin))**2)*(1-mean*reciproke)
                END IF
            END DO
        END DO

    END SUBROUTINE
END MODULE
