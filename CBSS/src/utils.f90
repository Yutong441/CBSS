subroutine linspace(low, high, num, out_arr)
    implicit none
    real(8), intent(in) :: low, high
    integer, intent(in) :: num
    real(8), dimension(num), intent(out) :: out_arr

    integer :: i
    do i=0,(num-1)
        out_arr(i+1) = i*(high - low)/(num-1) + low
    end do
end subroutine linspace


subroutine fa_md(eig_val, fa, md)
    implicit none
    real(8), dimension(3), intent(in) :: eig_val
    real(8), intent(out) :: fa, md
    real(8) :: L1, L2, L3, numer, denom

    L1 = eig_val(1)
    L2 = eig_val(2)
    L3 = eig_val(3)

    md = (L1 + L2 + L3)/3
    numer = sqrt((L1 - md)**2 + (L2 - md)**2 + (L3 - md)**2)
    denom = sqrt(L1**2 + L2**2 + L3**2)
    fa = sqrt(1.5)*numer/(denom + 1e-5)
end subroutine fa_md


subroutine inv_mat(W, n, sig, invWTS2W_WTS2)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n, 7), intent(in) :: W
    real(8), dimension(n), intent(in) :: sig
    real(8), dimension(7, n), intent(out) :: invWTS2W_WTS2

    real(8), dimension(n, n) :: S2  ! signal variance
    real(8), dimension(7, n) :: WTS2
    real(8), dimension(7, 7) :: WTS2W
    real(8), dimension(7, 7) :: inv_WT_S2_W
    integer :: lwork, info
    real(8), dimension(:), allocatable :: work
    real(8), dimension(1,1) :: dummy
    real(8), dimension(7) :: sing
    real(8) :: max_sing
    real(8), dimension(7,7) :: u, vt, sing_inv
    real(8) :: invert_thres  ! threshold of singular value for inversion
    integer :: i
    external :: dgesvd  ! SVD to invert matrix

    invert_thres = 1e-14
    S2(:, :) = 0
    do i=1,n
        S2(i,i) = sig(i)**2
    end do

    ! Defining matrix to solve fwDTI wls solution
    WTS2 = matmul(transpose(W), S2)
    WTS2W = matmul(WTS2, W)

    ! invert WTS2W
    ! A = USV^T
    ! A^-1 = VS^-1U^T
    lwork = -1
    call dgesvd("a", "a", 7, 7, WTS2W, 7, sing, u, 7, vt, 7, dummy, lwork, info)
    lwork = max(133*7, nint(dummy(1,1)))

    allocate(work(lwork))
    call dgesvd("a", "a", 7, 7, WTS2W, 7, sing, u, 7, vt, 7, work, lwork, info)
    deallocate(work)

    sing_inv(:,:) = 0
    max_sing = maxval(sing)
    do i=1,7
        if (sing(i)/max_sing < invert_thres) then
            sing_inv(i,i) = 0
        else
            sing_inv(i,i) = 1/sing(i)
        end if
    end do

    inv_WT_S2_W = matmul(sing_inv, transpose(u))
    inv_WT_S2_W = matmul(transpose(vt), inv_WT_S2_W)
    invWTS2W_WTS2 = matmul(inv_WT_S2_W, WTS2)
end subroutine inv_mat


subroutine eig_dti(D, evals)
    implicit none
    real(8), dimension(6), intent(in) :: D
    real(8), dimension(3), intent(out) :: evals
    real(8), dimension(9) :: dti_mat
    real(8), dimension(3, 3) :: dti_mat2

    integer :: lwork, info
    real(8), dimension(:), allocatable :: work
    real(8) :: zero
    integer :: i
    external :: dsyev

    zero = 0
    dti_mat = (/D(1), D(2), D(4), D(2), D(3), D(5), D(4), D(5), D(6)/)
    ! dti_mat = (/D(1), D(2), D(4), zero, D(3), D(5), zero, zero, D(6)/)
    dti_mat2 = reshape(dti_mat, (/3, 3/))

    lwork = 12
    allocate(work(lwork))
    call dsyev('N', 'U', 3, dti_mat2, 3, evals, work, lwork, info)
    deallocate(work)

    ! set minimum eigen values to 0
    do i=1,3
        if (evals(i) < 0) then
            evals(i) = 0
        end if
    end do
end subroutine eig_dti


subroutine min_index(arr, n, scalar, ind)
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: arr
    real(8), intent(in) :: scalar
    integer, intent(out) :: ind

    integer :: i
    integer, dimension(1) :: ind_num
    real(8), dimension(n) :: diff
    do i=1,n
        diff(i) = abs(arr(i) - scalar)
    end do
    ind_num = minloc(diff)
    ind = ind_num(1)
end subroutine min_index


subroutine tensor_to_scalar(param4d, t, h, l, d, maps)
    ! compute FA and MD from tensor
    integer, intent(in) :: t, h, l, d
    real(8), dimension(t, h, l, d), intent(in) :: param4d
    real(8), dimension(h, l, d, 2), intent(out) :: maps
    integer :: i, j, k
    real(8), dimension(3) :: evals
    real(8) :: fa_one, md_one
    external :: eig_dti, fa_md

    do k=1,d
        do j=1,l
            do i=1,h
                call eig_dti(param4d(1:6, i, j, k), evals)
                call fa_md(evals, fa_one, md_one)
                maps(i,j,k,1) = fa_one
                maps(i,j,k,2) = md_one
            end do
        end do
    end do
end subroutine tensor_to_scalar


subroutine gauss2d(im, t, h, l, d, sigma, im_out)
    integer, intent(in) :: t, h, l, d
    real(8), dimension(t, h, l, d), intent(in) :: im
    real(8) :: sigma
    real(8), dimension(t, h, l, d), intent(out) :: im_out

    real(8), dimension(:,:), allocatable :: kern
    integer :: i, j, k, k_i, k_j, ind_i, ind_j, r
    real(8), dimension(t) :: g

    ! define kernel
    r = int(sigma*2)
    allocate(kern(2*r+1, 2*r+1))
    do k_j=-r,r
        do k_i=-r,r
            kern(k_i+r+1, k_j+r+1) = exp(-(k_i**2 + k_j**2)/(2*sigma**2))
        end do
    end do
    kern(:,:) = kern(:,:)/sum(kern)

    do k=1,d
        do j=1,l
            do i=1,h
                g = 0
                do k_j=-r,r
                    do k_i=-r,r
                       ind_i = i+k_i 
                       ind_j = j+k_j 
                       if (ind_i>=1 .and. ind_i<=h .and. ind_j>=1 .and. ind_j<=l) then
                           g(:) = g(:) + kern(k_i+r+1, k_j+r+1)*im(:,ind_i,ind_j,k)
                       end if
                    end do
                end do
                im_out(:,i,j,k) = g(:)
            end do
        end do
    end do
    deallocate(kern)
end subroutine gauss2d
