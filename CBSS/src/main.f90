subroutine find_fw_sing_shell(W, img4d, t, h, l, d, mask, S0, maps)
    ! Args:
    ! W: design matrix, combining bvals and bvecs into a single matrix, 
    ! can be obtained from `dipy.reconst.dti.design_matrix`
    ! img4d: 4D DWI signal
    ! mask: brain mask
    ! S0: averaged b0 image
    ! 
    ! Return:
    ! maps: 3 maps stacked along the fourth axis in the order of FW, 
    ! FW-corrected FA, FW-corrected MD
    integer, intent(in) :: t, h, l, d
    real(8), dimension(t, 7), intent(in) :: W
    real(8), dimension(t, h, l, d), intent(in) :: img4d
    integer, dimension(h, l, d), intent(in) :: mask
    real(8), dimension(h, l, d), intent(in) :: S0
    real(8), dimension(h, l, d, 3), intent(out) :: maps

    real(8), dimension(7, h, l, d) :: dti_params
    real(8), dimension(h, l, d, 2) :: dti_maps, fw_maps
    real(8), dimension(h, l, d) :: md_map
    real(8), dimension(9, h, l, d) :: fw_params
    real(8) :: mdreg, min_signal, MDm, Diso
    real(8) :: one_md, mcsf, mdreg1
    integer :: piterations
    integer :: i, j, k
    integer :: num_csf
    external :: wls_fit_dti, tensor_to_scalar, wls_fit_tensor_fw

    ! define parameters
    mdreg = 2.0e-3
    min_signal = 1e-6
    piterations = 2
    MDm = 0.0006
    Diso = 3e-3

    call wls_fit_dti(W, img4d, t, h, l, d, mask, min_signal, dti_params)
    call tensor_to_scalar(dti_params, 7, h, l, d, dti_maps)

    md_map(:, :, :) = dti_maps(:, :, :, 2)
    ! mean CSF signal
    mcsf = 0
    num_csf = 0
    do k=1,d
        do j=1,l
            do i=1,h
                one_md = md_map(i, j, k)
                if (one_md > 0.002) then
                    mcsf = mcsf + one_md
                    num_csf = num_csf + 1
                end if
            end do
        end do
    end do

    mcsf = mcsf/num_csf
    mdreg1 = 0.002*mCSF/0.0025
    mdreg = min(mdreg, mdreg1)

    call wls_fit_tensor_fw(W, img4d, t, h, l, d, md_map, mask, S0, Diso, min_signal, piterations, mdreg, MDm, fw_params)
    call tensor_to_scalar(fw_params, 9, h, l, d, fw_maps)

    maps(:, :, :, 1) = fw_params(8, :, :, :)
    maps(:, :, :, 2) = fw_maps(:, :, :, 1)
    maps(:, :, :, 3) = fw_maps(:, :, :, 2)
end subroutine find_fw_sing_shell 
