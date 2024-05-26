! fill gray matter voxels with the annotation from adjacent voxels
! use gray matter mask derived from freesurfer
subroutine fill_gray(mask, h, w, d, gray, vox_h, vox_w, vox_d, dist_thres, out_mask)
    ! Method:
    ! 1. for each voxel in the gray matter mask, identify which sulcus it belongs to
    ! 2. this is the label of the nearest voxels
    ! 3. If the voxel is not adjacent to another voxels that have annotation, the 
    ! voxel is not labelled
    ! Args:
    ! `mask`: sulcus annotation mask
    ! `gray`: binary mask of gray matter
    ! `vox_h`, `vox_w`, `vox_d`: voxel dimension
    ! `dist_thres`: distance threshold
    implicit none
    integer, intent(in) :: h, w, d
    real(8), intent(in) :: vox_h, vox_w, vox_d, dist_thres
    integer, dimension(h, w, d), intent(in) :: mask, gray
    integer, dimension(h, w, d), intent(out) :: out_mask

    integer :: i, j, k, label
    external :: nearest_vox_label
    out_mask(:, :, :) = 0
    do k=1,d
        do j=1,w
            do i=1,h
                if (gray(i, j, k) > 0) then
                    call nearest_vox_label(mask, h, w, d, i, j, k, vox_h, vox_w, vox_d, dist_thres, label)
                    out_mask(i, j, k) = label
                end if
            end do
        end do
    end do
end subroutine fill_gray


subroutine nearest_vox_label(mask, h, w, d, h0, w0, d0, vox_h, vox_w, vox_d, dist_thres, label)
    implicit none
    integer, intent(in) :: h, w, d, h0, w0, d0
    real(8), intent(in) :: vox_h, vox_w, vox_d, dist_thres
    integer, dimension(h, w, d), intent(in) :: mask
    integer, intent(out) :: label

    real(8) :: min_dist, dist
    integer :: best_label, mask_label
    integer :: i, j, k, step

    integer :: start_h, end_h, start_w, end_w, start_d, end_d

    min_dist = dist_thres
    best_label = 3000
    step = int(dist_thres, 4)

    start_h = max(1, h0 - step)
    end_h = min(h, h0 + step)
    start_w = max(1, w0 - step)
    end_w = min(w, w0 + step)
    start_d = max(1, d0 - step)
    end_d = min(d, d0 + step)

    do k=start_d,end_d
        do j=start_w,end_w
            do i=start_h,end_h
                dist = (vox_h*real(h0 - i, 8))**2 + (vox_w*real(w0 - j, 8))**2 + (vox_d*real(d0 - k, 8))**2
                mask_label = mask(i, j, k)
                if (dist < min_dist**2 .and. mask_label > 0) then
                    min_dist = dist
                    best_label = mask_label
                end if
            end do
        end do
    end do
    label = best_label
end subroutine nearest_vox_label
