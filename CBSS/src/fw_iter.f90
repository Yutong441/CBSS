subroutine process_fw(sig, n, W, MD, SS, Diso, min_signal, piterations, ns, MDm, fw_params)
    ! sig: diffusion signal from one voxel
    ! piterations: number of iterations
    ! W: design matrix, num x 7; Bxx, Bxy, Byy, Bxz, Byz, Bzz, -1
    ! Bxx = -b * gx * gx
    ! Bxy = -2b * gx * gy
    ! MDm: MD threshold
    ! SS: b0 signal
    ! Diso: diffusivity of water at 37^0C, which is 3e-3 mm^2/s
    implicit none
    integer, intent(in) :: n, piterations, ns
    real(8), dimension(n), intent(in) :: sig
    real(8), dimension(n, 7), intent(in) :: W
    real(8), intent(in) :: MD, SS, MDm, min_signal, Diso
    real(8), dimension(9), intent(out) :: fw_params

    ! intermediary variables
    real(8) :: df, flow, fhig
    real(8), dimension(ns) :: fs  ! each f value to test at each iteration
    real(8), dimension(ns) :: f2  ! store MSE of each iteration
    real(8) :: f2s1
    real(8), dimension(n) :: sig_pred ! predicted signal
    real(8), dimension(n) :: SA ! signal subtracted from free water diffusion 
    real(8), dimension(n, 1) :: y ! y variable in a least square problem
    real(8), dimension(7, 1) :: all_new_params

    real(8), dimension(7, n) :: invWTS2W_WTS2
    real(8), dimension(7) :: Dwater
    real(8), dimension(n) :: fwsig
    real(8), dimension(3) :: evals
    integer :: i, j, p

    real(8) :: fa_one, md_one
    real(8), dimension(ns) :: fa2, md2
    real(8) :: MDa, MDm2
    real(8), dimension(7, 1) :: params1
    real(8) :: f
    real(8) :: zero
    integer :: Mind1, Mind2
    external :: linspace, fa_md, inv_mat, min_index

    ! ------------------------------------------------------------
    call inv_mat(W, n, sig, invWTS2W_WTS2)
    ! exp(-bD)
    zero = 0.
    Dwater = (/Diso, zero, Diso, zero, zero, Diso, zero/)
    do i=1,n
        fwsig(i) = exp(dot_product(W(i, :), Dwater))
    end do

    ! ------------------------------------------------------------
    df = 1  ! initialize precision
    flow = 0  ! lower f evaluated
    fhig = 1  ! higher f evaluated
    do p=1,piterations
        df = df * 0.1
        call linspace(flow, fhig, ns, fs)  ! sampling f
        fs(ns) = 0.98

        do i=1,ns
            do j=1,n
                SA(j) = sig(j) - fs(i)*SS*fwsig(j)
                if (SA(j) < 0) then
                    SA(j) = min_signal  ! signal can only be positive
                end if
            end do

            y(:,1) = log(SA(:) / (1-fs(i)))
            all_new_params = matmul(invWTS2W_WTS2, y)

            ! Select params for lower f2 (RMSE)
            do j=1,n
                sig_pred(j) = (1-fs(i))*exp(dot_product(W(j, :), all_new_params(:, 1))) + fs(i)*SS*fwsig(j)
            end do
            f2(i) = sum((sig(:) - sig_pred(:))**2)

            ! perform eigen decomposition
            call eig_dti(all_new_params(1:6, 1), evals)
            call fa_md(evals, fa_one, md_one)
            fa2(i) = fa_one
            md2(i) = md_one
        end do

        if (p == 1) then
            MDa = MD2(1)
        end if

        ! select the optimal parameter
        if (MD > MDm) then
            ! Mind1 = minloc(abs(md2(:) - MDm))
            call min_index(md2, ns, MDm, Mind1)
            ! Mind2 = minloc(abs(fa2(:) - 3.0*fa2(1)))
            call min_index(fa2, ns, 3.0*fa2(1), Mind2)
            Mind1 = min(Mind1, Mind2)
        else
            MDm2 = 0.00042*MDa/MDm
            call min_index(md2, ns, MDm2, Mind1)
            call min_index(fa2, ns, 3.0*fa2(1), Mind2)
            Mind1 = min(Mind1, Mind2)
        end if

        f2s1 = f2(1) - f2(Mind1)
        ! params1 = all_new_params(:, Mind1)
        do j=1,n
            SA(j) = sig(j) - fs(Mind1)*SS*fwsig(j)
            if (SA(j) < 0) then
                SA(j) = min_signal  ! signal can only be positive
            end if
        end do
        y(:,1) = log(SA(:) / (1-fs(Mind1)))
        params1 = matmul(invWTS2W_WTS2, y)

        f = fs(Mind1)  ! Updated f

        ! search near the regions of optimal f
        flow = max(f - df, 0.)
        fhig = min(f + df, 0.98)
        fw_params(1:7) = params1(:, 1)
        fw_params(8) = f
        fw_params(9) = f2s1
    end do
end subroutine process_fw


subroutine wls_iter_fw(W, n, sig, MD, S0, Diso, mdreg, min_signal, piterations, MDm, fw_params)
    implicit none
    integer, intent(in) :: n, piterations
    real(8), dimension(n, 7), intent(in) :: W
    real(8), dimension(n), intent(in) :: sig
    real(8), intent(in) :: MD, S0, Diso, mdreg, min_signal, MDm
    real(8), dimension(9), intent(out) :: fw_params

    external :: process_fw 
    real(8) :: mean_sig
    integer :: ns = 21 ! number of samples tested per iteration

    mean_sig = sum(sig)/n
    if (MD < mdreg) then
        if (mean_sig > min_signal .and. S0 > min_signal) then
            ! Process voxel if it has significant signal from tissue
            call process_fw(sig, n, W, MD, S0, Diso, min_signal, piterations, ns, MDm, fw_params)
        else
            ! in regions of very low signals, free water fraction must be 0
            fw_params(:) = 0
        end if
    else
        ! in regions of very high MD, free water fraction must be 1
        fw_params(:) = 0
        fw_params(8) = 1.
    end if
end subroutine wls_iter_fw


subroutine wls_fit_tensor_fw(W, img, t, h, l, d, md_data, mask, S0, Diso, min_signal, piterations, mdreg, MDm, param4d)
    integer, intent(in) :: t, h, l, d
    real(8), dimension(t, 7), intent(in) :: W
    real(8), dimension(t, h, l, d), intent(in) :: img
    real(8), dimension(h, l, d), intent(in) :: md_data, S0
    integer, dimension(h, l, d), intent(in) :: mask
    real(8), intent(in) :: Diso, min_signal, mdreg, MDm
    integer, intent(in) :: piterations
    real(8), dimension(9, h, l, d), intent(out) :: param4d

    real(8), dimension(t) :: sig
    real(8) :: md_one, s0_one
    integer :: i, j, k
    real(8), dimension(9) :: fw_params
    external :: wls_iter_fw

    do k=1,d
        do j=1,l
            do i=1,h
                if (mask(i,j,k) > 0) then
                    sig = img(:,i,j,k)
                    md_one = md_data(i,j,k)
                    s0_one = S0(i,j,k)
                    call wls_iter_fw(W, t, sig, md_one, s0_one, Diso, mdreg, min_signal, piterations, MDm, fw_params)
                    param4d(:,i,j,k) = fw_params(:)
                else
                    param4d(:,i,j,k) = 0
                end if
            end do
        end do
    end do
end subroutine wls_fit_tensor_fw
! Diso=3e-3, mdreg=2e-3, min_signal=1e-6, piterations=2, MDm=0.0006
