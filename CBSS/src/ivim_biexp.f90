! based on https://github.com/dipy/dipy/blob/master/dipy/reconst/ivim.py
! fit a biexponential model to IVIM data

! ! ====================Linear fitting====================
! subroutine lingress_full(x, n, y, slope, intercept)
!     implicit none
!     integer, intent(in) :: n
!     real(8), dimension(n), intent(in) :: x, y
!     real(8), intent(out) :: slope, intercept
!     real(8) :: mu_x, mu_y, var_x
! 
!     mu_x = sum(x)/n
!     mu_y = sum(y)/n
!     var_x = sum((X(:) - mu_x)**2)
!     slope = sum((X(:) - mu_x)*(Y(:) - mu_y))/var_x
!     intercept = mu_y - mu_x*slope
! end subroutine lingress_full
! 
! 
! subroutine initial_fit(bvals, t, signal, thres, signs, D, logS0)
!     implicit none
!     integer, intent(in) :: t
!     real(8), dimension(t), intent(in) :: bvals, signal
!     real(8), intent(in) :: thres
!     real(8), intent(in) :: signs
!     ! -1 for b<thres, +1 for b>thres
!     real(8), intent(out) :: D, logS0
! 
!     real(8), dimension(:), allocatable :: b_thres, sig_thres
!     integer :: len_b, index_b, i
!     external :: lingress_full
! 
!     len_b = 0
!     do i=1,t
!         if (signs*bvals(i) >= signs*thres) then
!             len_b = len_b + 1
!         end if
!     end do
! 
!     allocate(sig_thres(len_b))
!     allocate(b_thres(len_b))
! 
!     index_b = 0
!     do i=1,t
!         if (signs*bvals(i) >= signs*thres) then
!             index_b = index_b + 1
!             sig_thres(index_b) = signal(i)
!             b_thres(index_b) = bvals(i)
!         end if
!     end do
! 
!     sig_thres(:) = log(sig_thres(:))
!     call lingress_full(b_thres, len_b, sig_thres, D, logS0)
!     D = -D
!     deallocate(sig_thres)
!     deallocate(b_thres)
! end subroutine initial_fit
! 
! 
! ! ====================Initial bi-exponential fitting====================
! subroutine biexp_pred(b, t, logS0, D, D_star, f, S)
!     implicit none
!     integer, intent(in) :: t
!     real(8), dimension(t), intent(in) :: b
!     real(8), intent(in) :: logS0, D, D_star, f
!     real(8), dimension(t), intent(out) :: S
!     S = exp(logS0) * (f * exp(-b * D_star) + (1 - f) * exp(-b * D))
! end subroutine biexp_pred
! 
! 
! subroutine biexp_fit_initial(x_true, t, y_true, iter, out_vars)
!     ! estimate f and D*
!     use bobyqa_module, only: bobyqa
!     implicit none
!     integer, intent(in) :: t, iter
!     real(8), dimension(t), intent(in) :: x_true, y_true
!     real(8), dimension(4), intent(out) :: out_vars
! 
!     real(8), dimension(2) :: upper, lower, vars
!     real(8) :: rhobeg, rhoend
!     real(8) :: D, logS0, D_star, logS0_prime, frac
! 
!     external :: biexp_pred, initial_fit
! 
!     ! use linear regression to obtain initial estimates
!     call initial_fit(x_true, t, y_true, real(400., 8), real(1., 8), D, logS0)
!     D = max(D, 0.)
!     D = min(D, 0.01)
!     call initial_fit(x_true, t, y_true, real(200., 8), real(-1., 8), D_star, logS0_prime)
!     frac = 1 - exp(logS0)/exp(logS0_prime)
! 
!     ! set initial values
!     vars = (/D_star, frac/)
!     ! set bounds
! 
!     upper = (/0.15, 0.2/)
!     lower = (/1.5, 0./)
!     lower(1) = 1.5*D
!     rhobeg = 0.04
!     rhoend = 0.0004
! 
!     ! write(*, *) frac
!     ! write(*, *) D_star
!     ! write(*, *) D
!     ! write(*, *) S0
! 
!     call bobyqa(2, 5, vars, lower, upper, rhobeg, rhoend, 0, iter, calfun)
!     out_vars(1) = logS0_prime
!     out_vars(2) = D
!     out_vars(3) = D_star
!     out_vars(4) = frac
! 
!     contains
!         subroutine biexp_error(params, x, y, t, error)
!             implicit none
!             integer, intent(in) :: t
!             real(8), dimension(2), intent(in) :: params
!             real(8), dimension(t), intent(in) :: x, y
!             real(8), intent(out) :: error
!             real(8), dimension(t) :: y_pred
! 
!             call biexp_pred(x, t, logS0_prime, D, params(1), params(2), y_pred)
!             error = sum((y_pred - y)**2)/t
!         end subroutine biexp_error
! 
!         subroutine calfun(n, x, f)
!             integer,intent(in) :: n
!             real(8),dimension(:),intent(in) :: x
!             real(8),intent(out) :: f
!             call biexp_error(x, x_true, y_true, t, f)
!         end subroutine calfun
! end subroutine biexp_fit_initial
! 
! ! ====================Final bi-exponential fitting====================
! subroutine biexp_fit_final(x_true, t, y_true, iter, init_vars, final_vars)
!     ! estimate f and D*
!     use bobyqa_module, only: bobyqa
!     implicit none
!     integer, intent(in) :: t, iter
!     real(8), dimension(t), intent(in) :: x_true, y_true
!     real(8), dimension(4), intent(in) :: init_vars
!     real(8), dimension(4), intent(out) :: final_vars
! 
!     real(8), dimension(4) :: upper, lower, vars
!     real(8) :: rhobeg, rhoend
!     real(8) :: D, logS0_prime, D_star, frac
! 
!     external :: biexp_pred
! 
!     ! set initial values
!     vars(:) = init_vars(:)
!     logS0_prime = vars(1)
!     D = vars(2)
!     D_star = vars(3)
!     frac = vars(4)
! 
!     rhobeg = 0.04
!     rhoend = 0.0001
! 
!     ! set bounds
!     upper = (/0.5, 1500., 0.15, 0.2/)
!     upper(1) = 1.3*logS0_prime
!     upper(2) = 1500.*D
!     lower = (/1.5, 0.1, 1.5, 0./)
!     lower(1) = 0.7*logS0_prime
!     lower(3) = D*1.5
! 
!     vars(3) = max(vars(3), lower(3))
!     vars(2) = max(vars(2), lower(2))
! 
!     call bobyqa(4, 9, vars, lower, upper, rhobeg, rhoend, 0, iter, calfun)
!     final_vars(:) = vars(:)
!     final_vars(2) = vars(2)/1000.
!     contains
!         subroutine biexp_error(params, x, y, t, error)
!             implicit none
!             integer, intent(in) :: t
!             real(8), dimension(4), intent(in) :: params
!             real(8), dimension(t), intent(in) :: x, y
!             real(8), intent(out) :: error
!             real(8), dimension(t) :: y_pred
! 
!             call biexp_pred(x, t, params(1), params(2)/1000., params(3), params(4), y_pred)
!             error = sum((y_pred - y)**2)/t
!         end subroutine biexp_error
! 
!         subroutine calfun(n, x, f)
!             integer,intent(in) :: n
!             real(8),dimension(:),intent(in) :: x
!             real(8),intent(out) :: f
!             call biexp_error(x, x_true, y_true, t, f)
!         end subroutine calfun
! end subroutine biexp_fit_final
! 
! 
! subroutine biexp_fit_all(x_true, t, y_true, iter, final_vars)
!     implicit none
!     integer, intent(in) :: t, iter
!     real(8), dimension(t), intent(in) :: x_true, y_true
!     real(8), dimension(4), intent(out) :: final_vars
! 
!     real(8), dimension(4) :: init_vars
!     external :: biexp_fit_initial, biexp_fit_final
!     call biexp_fit_initial(x_true, t, y_true, iter, init_vars)
!     call biexp_fit_final(x_true, t, y_true, iter, init_vars, final_vars)
! end subroutine biexp_fit_all
! 
! 
! subroutine biexp_fit_all4d(img4d, t, x, y, z, bval, masks, maps)
!     ! Args:
!     ! `img4d`: motion corrected 4D IVIM data
!     ! `bval`: a vector of b values in the same order as the IVIM data
!     ! `masks`: brain mask. IVIM fitting is relatively computationally 
!     ! intensive, using a mask can reduce computation time
!     ! `maps`: H x W X D x 3, the maps are stored in the order of D, D*, f
!     !
!     ! Method:
!     ! 1. Use linear regression to estimate D, D* and f
!     ! 2. Use Powell method to refine the estimate of D* and f
!     integer, intent(in) :: t, x, y, z
!     real(8), dimension(t, x, y, z), intent(in) :: img4d
!     real(8), dimension(t), intent(in) :: bval
!     integer, dimension(x, y, z), intent(in) :: masks
!     real(8), dimension(x, y, z, 3), intent(out) :: maps
!     integer :: i, j, k
!     real(8), dimension(t) :: signal
!     real(8), dimension(4) :: one_vars
!     external :: biexp_fit_all
! 
!     maps(:, :, :, :) = 0
!     do k=1,z
!         do j=1,y
!             do i=1,x
!                 if (masks(i, j, k) > 0) then
!                     signal(:) = img4d(:, i, j, k)
!                     call biexp_fit_initial(bval, t, signal, 200, one_vars)
!                     maps(i, j, k, :) = one_vars(2:4)
!                 end if
!             end do
!         end do
!     end do
! end subroutine biexp_fit_all4d
