MODULE reserve_mod
  IMPLICIT NONE
  INTERFACE reserve
    MODULE PROCEDURE reserve_i_d1
    MODULE PROCEDURE reserve_i_d2
    MODULE PROCEDURE reserve_l_d1
    MODULE PROCEDURE reserve_r8_d1
    MODULE PROCEDURE reserve_r8_d2
    MODULE PROCEDURE reserve_r8_d3
    MODULE PROCEDURE reserve_r8_d4
    MODULE PROCEDURE reserve_r8_lb2
    MODULE PROCEDURE reserve_r8_lb3
  END INTERFACE
  
CONTAINS

  SUBROUTINE reserve_i_d1(a, n, offload)
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: a(:)
    INTEGER, INTENT(IN) :: n
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a, 1) < n) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_i_d2(a, n, offload)
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: a(:, :)
    INTEGER, INTENT(IN) :: n(2)
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a, 1) < n(1) .OR. SIZE(a, 2) < n(2)) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n(1), n(2)))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_l_d1(a, n, offload)
    IMPLICIT NONE
    LOGICAL, ALLOCATABLE :: a(:)
    INTEGER, INTENT(IN) :: n
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a, 1) < n) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_r8_d1(a, n, offload)
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: a(:)
    INTEGER, INTENT(IN) :: n
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a) < n) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_r8_d2(a, n, offload)
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: a(:, :)
    INTEGER, INTENT(IN) :: n(2)
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a, 1) < n(1) .OR. SIZE(a, 2) < n(2)) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n(1), n(2)))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_r8_d3(a, n, offload)
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: a(:, :, :)
    INTEGER, INTENT(IN) :: n(3)
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a, 1) < n(1) .OR. SIZE(a, 2) < n(2) &
        .OR. SIZE(a, 3) < n(3)) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n(1), n(2), n(3)))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_r8_d4(a, n, offload)
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: a(:, :, :, :)
    INTEGER, INTENT(IN) :: n(4)
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (SIZE(a, 1) < n(1) .OR. SIZE(a, 2) < n(2) &
        .OR. SIZE(a, 3) < n(3) .OR. SIZE(a, 4) < n(4)) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(n(1), n(2), n(3), n(4)))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_r8_lb2(a, lb, ub, offload)
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: a(:, :)
    INTEGER, INTENT(IN) :: lb(2), ub(2)
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (LBOUND(a, 1) /= lb(1) .OR. LBOUND(a, 2) /= lb(2) &
        .OR. UBOUND(a, 1) < ub(1) .OR. UBOUND(a, 2) < ub(2)) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(lb(1):ub(1), lb(2):ub(2)))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

  SUBROUTINE reserve_r8_lb3(a, lb, ub, offload)
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: a(:, :, :)
    INTEGER, INTENT(IN) :: lb(3), ub(3)
    LOGICAL, OPTIONAL, INTENT(IN) :: offload
    LOGICAL :: off
    off = .FALSE.
    IF (PRESENT(offload)) off = offload
    IF (ALLOCATED(a)) THEN
      IF (LBOUND(a, 1) /= lb(1) .OR. LBOUND(a, 2) /= lb(2) .OR. LBOUND(a, 3) /= lb(3) &
        .OR. UBOUND(a, 1) < ub(1) .OR. UBOUND(a, 2) < ub(2) .OR. UBOUND(a, 3) < ub(3)) THEN
        IF (off) THEN
          !$omp target exit data map(delete: a)
        END IF
        DEALLOCATE(a)
      END IF
    END IF
    IF (.NOT. ALLOCATED(a)) THEN
      ALLOCATE(a(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
      IF (off) THEN
        !$omp target enter data map(alloc: a)
      END IF
    END IF
  END SUBROUTINE

END MODULE

