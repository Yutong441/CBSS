subroutine wls_iter_dti(W, n, sig, min_signal, dti_params)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, 7), intent(in) :: W
    real(8), dimension(n), intent(in) :: sig
    real(8), intent(in) :: min_signal
    real(8), dimension(7), intent(out) :: dti_params
    external :: inv_mat

    real(8), dimension(7, n) :: invWTS2W_WTS2
    real(8), dimension(n, 1) :: sig_pos
    real(8), dimension(7, 1) :: params
    integer :: i

    call inv_mat(W, n, sig, invWTS2W_WTS2)

    do i=1,n
        if (sig(i) > 0) then
            sig_pos(i, 1) = sig(i)
        else
            sig_pos(i, 1) = min_signal
        end if
    end do

    params = matmul(invWTS2W_WTS2, log(sig_pos))
    dti_params(:) = params(:, 1)
end subroutine wls_iter_dti
 
 
subroutine wls_fit_dti(W, img, t, h, l, d, mask, min_signal, param4d)
    implicit none
    integer, intent(in) :: t, h, l, d
    real(8), dimension(t, 7), intent(in) :: W
    real(8), dimension(t, h, l, d), intent(in) :: img
    integer, dimension(h, l, d), intent(in) :: mask
    real(8), intent(in) :: min_signal
    real(8), dimension(7, h, l, d), intent(out) :: param4d

    real(8), dimension(t) :: sig
    integer :: i, j, k
    real(8), dimension(7) :: params
    external :: wls_iter_dti

    do k=1,d
        do j=1,l
            do i=1,h
                if (mask(i,j,k) > 0) then
                    sig = img(:,i,j,k)
                    call wls_iter_dti(W, t, sig, min_signal, params)
                    param4d(:,i,j,k) = params(:)
                else
                    param4d(:,i,j,k) = 0
                end if
            end do
        end do
    end do
end subroutine wls_fit_dti
