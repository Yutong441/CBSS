import numpy as np
import nibabel as nib
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.dti import design_matrix
from src import fort
from scipy.ndimage import gaussian_filter


def find_fw_sing_shell(fn_data, fn_bval, fn_bvec, save_dir, fn_mask=None):
    '''
    This function translates the original free water calculation script from
    python into fortran. It can achieve 10x speedup.
    3 files will be saved:
    FW.nii.gz, FA.nii.gz, MD.nii.gz
    '''
    nii = nib.load(fn_data)
    affine = nii.affine
    img_data = nii.get_fdata()
    del nii
    bvals, bvecs = read_bvals_bvecs(fn_bval, fn_bvec)

    if fn_mask is not None:
        niim = nib.load(fn_mask)
        brain_mask = niim.get_fdata()
    else:
        brain_mask = np.ones(img_data.shape[:3])

    # Construct the gradient table
    gtab = gradient_table(bvals, bvecs)
    W = design_matrix(gtab)

    # smooth the data
    fwhm = 1.25
    gauss_std = fwhm / np.sqrt(8 * np.log(2))
    # converting fwhm to Gaussian std
    data_smooth = np.zeros(img_data.shape)
    for v in range(img_data.shape[-1]):
        data_smooth[..., v] = gaussian_filter(img_data[..., v],
                                              sigma=gauss_std)

    S0 = np.mean(data_smooth[..., gtab.b0s_mask], axis=-1)
    img4d = np.asfortranarray(np.moveaxis(data_smooth, 3, 0))
    maps = fort.find_fw_sing_shell(W, img4d, brain_mask, S0)

    FW = nib.Nifti1Image(maps[..., 0], affine=affine)
    FW.to_filename(save_dir+"/FW.nii.gz")
    FA = nib.Nifti1Image(maps[..., 1], affine=affine)
    FA.to_filename(save_dir+"/FA.nii.gz")
    MD = nib.Nifti1Image(maps[..., 2], affine=affine)
    MD.to_filename(save_dir+"/MD.nii.gz")
