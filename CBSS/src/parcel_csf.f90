! parcellation of CSF
subroutine parcellate(csf_mask, h, w, d, cortex, x_res, y_res, z_res, dist_cutoff, csf_parc)
    ! cortex: a label image from 1 to n
    ! dist_cutoff: if above a certain distance, no neighbor is found
    ! then the point is classified as 0
    implicit none
    integer, intent(in) :: h, w, d, dist_cutoff
    integer, dimension(h, w, d), intent(in) :: csf_mask, cortex
    integer, dimension(h, w, d), intent(out) :: csf_parc
    real(8), intent(in) :: x_res, y_res, z_res

    integer :: i, j, k, label
    external :: closest_neighbor

    csf_parc(:, :, :) = 0
    do k=1,d
        do j=1,w
            do i=1,h
                if (csf_mask(i, j, k) > 0) then
                    call closest_neighbor(cortex, h, w, d, i, j, k, x_res, y_res, z_res, dist_cutoff, label)
                    csf_parc(i, j, k) = label
                end if
            end do
        end do
    end do
end subroutine parcellate


subroutine closest_neighbor(cortex, h, w, d, h0, w0, d0, x_res, y_res, z_res, dist_cutoff, label)
    implicit none
    integer, intent(in) :: h, w, d, h0, w0, d0, dist_cutoff
    integer, dimension(h, w, d), intent(in) :: cortex
    integer, intent(out) :: label
    real(8), intent(in) :: z_res, x_res, y_res

    integer :: i, j, k
    real(8) :: dist, best_dist
    integer :: best_label

    integer :: start_h, start_w, start_d
    integer :: end_h, end_w, end_d

    best_dist = dist_cutoff**2
    best_label = 0

    ! refine search space
    start_h = max(1, h0 - int(dist_cutoff/x_res, 4))
    start_w = max(1, w0 - int(dist_cutoff/y_res, 4))
    start_d = max(1, d0 - int(dist_cutoff/z_res, 4))
    end_h = min(h, h0 + int(dist_cutoff/x_res, 4))
    end_w = min(w, w0 + int(dist_cutoff/y_res, 4))
    end_d = min(d, d0 + int(dist_cutoff/z_res, 4))

    do k=start_d,end_d
        do j=start_w,end_w
            do i=start_h,end_h
                if (cortex(i, j, k) > 0) then
                    dist = (x_res**2)*real((i - h0)**2, 8) + (y_res**2)*real((j - w0)**2, 8) + (z_res**2)*real((k - d0)**2, 8)
                    if (dist < best_dist) then
                        best_dist = dist
                        best_label = cortex(i, j, k)
                    end if
                end if
            end do
        end do
    end do
    label = best_label
end subroutine closest_neighbor
