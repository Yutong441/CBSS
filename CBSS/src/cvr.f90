subroutine correlation(x, n, y, corr)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) ::  x, y
    real(8), intent(out) ::  corr

    real(8) :: mu_x, mu_y, numer, denom
    
    mu_x = sum(x)/n
    mu_y = sum(y)/n
    numer = sum((x(:) - mu_x)*(y(:) - mu_y))
    denom = sqrt(sum((x(:) - mu_x)**2)*sum((y(:) - mu_y)**2))
    corr = numer/denom
end subroutine correlation

subroutine time_shift(one_sig, n, ref_sig, max_delay, delay)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: one_sig, ref_sig
    integer, intent(in) :: max_delay
    integer, intent(out) :: delay

    integer :: i, len_vec
    real(8) :: corr, best_corr
    real(8), dimension(:), allocatable :: ref_i, one_i
    external :: correlation

    best_corr = 0
    delay = 0
    do i=-max_delay,max_delay
        len_vec = n - abs(i)
        allocate(ref_i(len_vec))
        allocate(one_i(len_vec))
        if (i>=0) then
            ref_i(:) = ref_sig(1:len_vec)
            one_i(:) = one_sig((1+i):(len_vec + i))
        else
            ref_i(:) = ref_sig((1+abs(i)):len_vec + abs(i))
            one_i(:) = one_sig(1:len_vec)
        end if
        call correlation(ref_i, len_vec, one_i, corr)
        deallocate(ref_i)
        deallocate(one_i)

        if (corr > best_corr) then
            best_corr = corr
            delay = i
        end if
    end do
end subroutine time_shift


subroutine lingress_beta1(x, n, y, b1)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: x, y
    real(8), intent(out) :: b1
    real(8) :: mu_x, mu_y, var_x

    mu_x = sum(x)/n
    mu_y = sum(y)/n
    var_x = sum((X(:) - mu_x)**2)
    b1 = sum((X(:)-mu_x)*(Y(:)-mu_y))/var_x
end subroutine lingress_beta1


subroutine lag_regression(one_sig, n, ref_sig, max_delay, delay, b1, b1_shift)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: one_sig, ref_sig
    integer, intent(in) :: max_delay
    integer, intent(out) :: delay
    real(8), intent(out) :: b1, b1_shift
    real(8), dimension(:), allocatable :: ref_i, one_i
    external :: lingress_beta1, time_shift
    integer :: len_vec

    call lingress_beta1(ref_sig, n, one_sig, b1)
    call time_shift(one_sig, n, ref_sig, max_delay, delay)

    len_vec = n - abs(delay)
    allocate(ref_i(len_vec))
    allocate(one_i(len_vec))
    if (delay>=0) then
        ref_i(:) = ref_sig(1:len_vec)
        one_i(:) = one_sig((1+delay):(len_vec + delay))
    else
        ref_i(:) = ref_sig((1+abs(delay)):len_vec + abs(delay))
        one_i(:) = one_sig(1:len_vec)
    end if
    call lingress_beta1(ref_i, len_vec, one_i, b1_shift)
end subroutine lag_regression


subroutine lag_regression4d(img4d, t, x, y, z, mask, ref_sig, max_delay, maps)
    implicit none
    integer, intent(in) :: t, x, y, z
    real(8), dimension(t, x, y, z), intent(in) :: img4d
    integer, dimension(x, y, z), intent(in) :: mask
    real(8), dimension(t), intent(in) :: ref_sig
    integer, intent(in) :: max_delay
    real(8), dimension(x, y, z, 3), intent(out) :: maps

    integer :: i, j, k
    real(8), dimension(t) :: one_signal
    integer :: delay
    real(8) :: b1, b1_shift
    external :: lag_regression

    maps(:, :, :, :) = 0
    do k=1,z
        do j=1,y
            do i=1,x
                if (mask(i, j, k) > 0) then
                    one_signal(:) = img4d(:, i, j, k)
                    call lag_regression(one_signal, t, ref_sig, max_delay, delay, b1, b1_shift)
                    maps(i, j, k, 1) = delay
                    maps(i, j, k, 2) = b1
                    maps(i, j, k, 3) = b1_shift
                end if
            end do
        end do
    end do
end subroutine lag_regression4d
