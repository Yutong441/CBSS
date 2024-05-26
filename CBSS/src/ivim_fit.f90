! Estimate ADC from IVIM data
subroutine lingress(x, n, y, slope)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: x, y
    real(8), intent(out) :: slope
    slope = sum(x(:)*y(:))/sum(x**2)

    ! I decided against adding an intercept estimate
    ! because ln(S/S0) =  -bD
    ! If I input only the S/S0 as y, there is no need for intercept
    ! 
    ! real(8) :: mu_x, mu_y, var_x
    ! mu_x = sum(x)/n
    ! mu_y = sum(y)/n
    ! var_x = sum((X(:) - mu_x)**2)
    ! slope = sum((X(:) - mu_x)*(Y(:) - mu_y))/var_x
end subroutine lingress


subroutine fit_ivim(img4d, t, h, l, d, bvals, s0, mask, diffu)
    implicit none
    integer, intent(in) :: t, h, l, d
    real(8), dimension(t, h, l, d), intent(in) :: img4d
    real(8), dimension(t), intent(in) :: bvals
    integer, dimension(h, l, d), intent(in) :: mask
    real(8), dimension(h, l, d), intent(in) :: s0
    real(8), dimension(h, l, d), intent(out) :: diffu

    real(8), dimension(t) :: y
    real(8) :: slope
    integer :: i, j, k
    external :: lingress

    diffu(:, :, :) = 0
    do k=1,d
        do j=1,l
            do i=1,h
                if (mask(i, j, k) > 0) then
                    y = img4d(:, i, j, k)
                    y(:) = -log(y(:)/s0(i, j, k))
                    call lingress(bvals, t, y, slope)
                    diffu(i, j, k) = slope
                end if
            end do
        end do
    end do
end subroutine
