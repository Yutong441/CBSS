subroutine lingress_glm(x, n, m, y, b)
    implicit none
    integer, intent(in) :: n, m
    real(8), dimension(n, m), intent(in) :: x
    real(8), dimension(n), intent(in) :: y
    real(8), dimension(m,1), intent(out) :: b

    real(8), dimension(m,m) :: xx
    real(8), dimension(m) :: sing
    real(8), dimension(m,m) :: sing_mat
    real(8), dimension(m,m) :: u, vt
    integer :: lwork, info
    real(8), dimension(:), allocatable :: work
    real(8), dimension(1,1) :: dummy

    real(8), dimension(m,m) :: inv_xx
    real(8), dimension(n,1) :: y_sub
    integer :: i
    external :: dgesvd

    xx = matmul(transpose(x), x)

    lwork = -1
    call dgesvd("a", "a", m, m, xx, m, sing, u, m, vt, m, dummy, lwork, info)
    lwork = max(133*n, nint(dummy(1,1)))

    allocate(work(lwork))
    call dgesvd("a", "a", m, m, xx, m, sing, u, m, vt, m, dummy, lwork, info)
    deallocate(work)

    sing_mat(:, :) = 0
    do i=1,m
        sing_mat(i, i) = 1/sing(i)
    end do
    inv_xx = matmul(matmul(transpose(vt), sing_mat), transpose(u))
    y_sub(:, 1) = y(:)
    b = matmul(matmul(inv_xx, transpose(x)), y_sub)
end subroutine lingress_glm


subroutine glm(one_sig, n, ref_sig, confound, m, outb)
    implicit none
    integer, intent(in) :: n, m
    real(8), dimension(n), intent(in) :: one_sig, ref_sig
    real(8), dimension(n, m), intent(in) :: confound
    real(8), intent(out) :: outb

    integer :: m3, i
    real(8), dimension(m, 1) :: b1
    real(8), dimension(n, m+3) :: x
    external :: lingress_glm

    m3 = m + 3
    x(:, 1) = ref_sig(:)
    do i=1,m
        x(:, i+1) = confound(:, i)
    end do

    ! add the signal time course as an additional regressor
    do i=1,n
        x(i, m+2) = i
        x(i, m+3) = 1
    end do
    call lingress_glm(x, n, m3, one_sig, b1)
    outb = b1(1, 1)  ! /b1(m3, 1)
end subroutine glm


subroutine get_mat(confound, n, m, outmat)
    implicit none
    integer, intent(in) :: n, m
    real(8), dimension(n, m), intent(in) :: confound
    real(8), dimension(m, n), intent(out) :: outmat

    integer :: i
    real(8), dimension(m,m) :: xx
    real(8), dimension(m) :: sing
    real(8), dimension(m,m) :: sing_mat
    real(8), dimension(m,m) :: u, vt
    integer :: lwork, info
    real(8), dimension(:), allocatable :: work
    real(8), dimension(1,1) :: dummy

    real(8), dimension(m,m) :: inv_xx
    external :: dgesvd

    ! ----------decompose matrix----------
    xx = matmul(transpose(confound), confound)
    ! write(*,*) xx

    lwork = -1
    call dgesvd("a", "a", m, m, xx, m, sing, u, m, vt, m, dummy, lwork, info)
    lwork = max(133*m, nint(dummy(1,1)))

    write(*,*) "step1"
    allocate(work(lwork))
    call dgesvd("a", "a", m, m, xx, m, sing, u, m, vt, m, dummy, lwork, info)
    deallocate(work)

    write(*,*) "step2"
    sing_mat(:, :) = 0
    do i=1,m
        sing_mat(i, i) = 1/sing(i)
    end do
    inv_xx = matmul(matmul(transpose(vt), sing_mat), transpose(u))
    outmat = matmul(inv_xx, transpose(confound))
end subroutine get_mat


subroutine get_refmat(ref_sig, n, confound, m, outmat)
    implicit none
    integer, intent(in) :: n, m
    real(8), dimension(n), intent(in) :: ref_sig
    real(8), dimension(n, m), intent(in) :: confound
    real(8), dimension(m+3, n), intent(out) :: outmat

    integer :: m3, i
    real(8), dimension(n, m+3) :: x
    real(8), dimension(m+3,m+3) :: xx
    real(8), dimension(m+3) :: sing
    real(8), dimension(m+3,m+3) :: sing_mat
    real(8), dimension(m+3,m+3) :: u, vt
    integer :: lwork, info
    real(8), dimension(:), allocatable :: work
    real(8), dimension(1,1) :: dummy

    real(8), dimension(m+3,m+3) :: inv_xx
    external :: dgesvd

    ! ----------compile the matrix----------
    m3 = m + 3
    x(:, 1) = ref_sig(:)
    do i=1,m
        x(:, i+1) = confound(:, i)
    end do

    ! add the signal time course as an additional regressor
    do i=1,n
        x(i, m+2) = i
        x(i, m+3) = 1
    end do

    ! ----------decompose matrix----------
    xx = matmul(transpose(x), x)

    lwork = -1
    call dgesvd("a", "a", m3, m3, xx, m3, sing, u, m3, vt, m3, dummy, lwork, info)
    lwork = max(133*m3, nint(dummy(1,1)))

    write(*,*) "step1"
    allocate(work(lwork))
    call dgesvd("a", "a", m3, m3, xx, m3, sing, u, m3, vt, m3, dummy, lwork, info)
    deallocate(work)

    write(*,*) "step2"
    sing_mat(:, :) = 0
    do i=1,m3
        sing_mat(i, i) = 1/sing(i)
    end do
    inv_xx = matmul(matmul(transpose(vt), sing_mat), transpose(u))
    outmat = matmul(inv_xx, transpose(x))
end subroutine get_refmat


subroutine get_b(glm_mat, m, n, y, outb)
    integer, intent(in) :: m, n
    real(8), dimension(m, n), intent(in) :: glm_mat
    real(8), dimension(n), intent(in) :: y
    real(8), intent(out) :: outb

    real(8), dimension(m, 1) :: bcoef
    real(8), dimension(n, 1) :: y_sub

    y_sub(:, 1) = y(:)
    bcoef = matmul(glm_mat, y_sub)

    outb = bcoef(1,1)
    ! if (bcoef(m,1) > 0) then
    !     outb = bcoef(1,1)/bcoef(m,1)
    ! else
    !     outb = 0
    ! end if
end subroutine get_b
