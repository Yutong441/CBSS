subroutine detrend_quad(img4d, t, x, y, z, mask, confound, m, img4d_de)
    ! quadratic detrending
    ! `mask`: mask in which detrending is carried out
    ! `confound`: matrix of confounding variables, e.g. motion parameters, CSF signals,
    ! WM signals
    implicit none
    integer, intent(in) :: t, x, y, z, m
    real(8), dimension(t, x, y, z), intent(in) :: img4d
    integer, dimension(x, y, z), intent(in) :: mask
    real(8), dimension(t, m), intent(in) :: confound
    real(8), dimension(t, x, y, z), intent(out) :: img4d_de

    integer :: i, j, k, m1
    real(8), dimension(t) :: residual
    real(8), dimension(t, 1) :: one_signal, y_pred
    real(8), dimension(m+1, 1) :: b
    real(8), dimension(t, m+1) :: confound_mat
    real(8), dimension(m+1, t) :: refmat
    external :: get_mat

    m1 = m + 1
    confound_mat(:, 1:m) = confound(:, :)
    do i=1,t
        confound_mat(i, m1) = 1
    end do

    call get_mat(confound_mat, t, m1, refmat)

    img4d_de(:, :, :, :) = 0
    do k=1,z
        do j=1,y
            do i=1,x
                if (mask(i, j, k) > 0) then
                    one_signal(:, 1) = img4d(:, i, j, k)
                    ! obtain coefficient
                    b = matmul(refmat, one_signal)
                    ! make prediction
                    y_pred = matmul(confound_mat, b)
                    ! add the baseline constant (last variable)
                    residual(:) = one_signal(:, 1) - y_pred(:, 1) + b(m1, 1)
                    img4d_de(:, i, j, k) = residual(:)
                end if
            end do
        end do
    end do
end subroutine detrend_quad
