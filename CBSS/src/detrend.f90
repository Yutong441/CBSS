subroutine detrend(x, n, y, y_de, b0)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: x, y
    real(8), dimension(n), intent(out) :: y_de
    real(8), intent(out) :: b0
    real(8) :: mu_x, mu_y, var_x, b1

    mu_x = sum(x)/n
    mu_y = sum(y)/n
    var_x = sum((X(:) - mu_x)**2)
    b1 = sum((X(:)-mu_x)*(Y(:)-mu_y))/var_x
    b0 = mu_y - mu_x*b1
    y_de = y - (b0 + b1*X)
end subroutine


subroutine detrend4d(img4d, t, x, y, z, mask, img4d_de)
    implicit none
    integer, intent(in) :: t, x, y, z
    real(8), dimension(t, x, y, z), intent(in) :: img4d
    integer, dimension(x, y, z), intent(in) :: mask
    real(8), dimension(t+1, x, y, z), intent(out) :: img4d_de

    integer :: i, j, k
    real(8), dimension(t) :: one_x, one_signal, y_de
    real(8) :: b0
    external :: detrend

    do i=1,t
        one_x(i)=i
    end do
    
    img4d_de(:, :, :, :) = 0
    do k=1,z
        do j=1,y
            do i=1,x
                if (mask(i, j, k) > 0) then
                    one_signal(:) = img4d(:, i, j, k)
                    call detrend(one_x, t, one_signal, y_de, b0)
                    img4d_de(1, i, j, k) = b0
                    img4d_de(2:(t + 1), i, j, k) = y_de(:)
                end if
            end do
        end do
    end do
end subroutine
