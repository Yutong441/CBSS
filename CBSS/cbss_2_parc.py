import os
import shutil
import numpy as np
from src import fort
import ants
import nibabel as nib
from cbss_1_reg import bash_in_python


def keep_GM(img):
    '''
    Only include gray matter
    exclude 2/41 (WM), 4/43 (lateral ventricle),
    14 (third ventricle), 15 (fourth ventricle)
    11/50 (caudate)
    12/51 (putamen) 13/52 (pallidum)
    10/49 (thalamus)
    28/60 (ventral DC)
    24 (CSF)
    '''
    gm = img.numpy()
    for i in [2, 41, 4, 5, 43, 44, 12, 13, 14, 15, 11, 50, 12, 51, 52, 24,
              10, 49, 28, 60]:
        gm *= (gm != i).astype(int)
    gm = ants.from_numpy(gm, origin=img.origin, direction=img.direction,
                         spacing=img.spacing)
    return gm


def remove_vent(img):
    img_np = img.numpy()
    csf = np.zeros(img.shape)
    for i in [4, 5, 43, 44, 13, 14]:
        csf += (img_np == i).astype(int)
    not_csf = ants.from_numpy(
        (csf == 0).astype(float), origin=img.origin, direction=img.direction,
        spacing=img.spacing)
    return not_csf


def keep_vent(img):
    img_np = img.numpy()
    csf = np.zeros(img.shape)
    for i in [4, 5, 43, 44, 14, 15]:
        csf += (img_np == i).astype(int)
    not_csf = ants.from_numpy(
        (csf >= 1).astype(float)*img_np,
        origin=img.origin, direction=img.direction, spacing=img.spacing)
    return not_csf


def parcellate_csf(root, free_water_thres=0.8, cutoff_dist=10):
    '''
    Parcellate peripheral CSF based on its adjacent gray matter
    Method:
    1. register freesurfer segmentation to native space (both are in native
    space, but the segmentation is of higher resolution)
    2. obtain peripheral CSF mask by:
        * remove skull
        * remove ventricle
        * remove lacunes
        * only free water-rich region is retained
    3. parcellate peripheral CSF
    '''
    pT2 = ants.image_read(root+"/pT2.nii.gz")
    MRI = ants.image_read(root+"/pT1_mri.nii.gz")
    seg = ants.image_read(root+"/pT2_seg.nii.gz")

    brain_mask = ants.threshold_image(seg, low_thresh=0.5, binary=True)
    MRI_brain = ants.mask_image(MRI, brain_mask)
    if not os.path.exists(root+"/tmp"):
        os.mkdir(root+"/tmp")
    tx = ants.registration(pT2, MRI_brain, type_of_transform="Rigid",
                           outprefix=root+"/tmp/ants")
    seg_out = ants.apply_transforms(pT2, seg, tx["fwdtransforms"],
                                    interpolator="genericLabel")
    ants.image_write(seg_out, root+"/seg.nii.gz")
    shutil.rmtree(root+"/tmp")

    seg_gray = keep_GM(seg_out)
    no_csf = remove_vent(seg_out)

    free_water = ants.image_read(root+"/FW.nii.gz")
    mask = ants.image_read(root+"/FW_synth_mask.nii.gz")
    # apply skull mask
    free_water = ants.mask_image(free_water, mask)
    # remove ventricle
    free_water = ants.mask_image(free_water, no_csf)
    # remove lacunes
    free_water = free_water.numpy()*(seg_out.numpy() == 24)
    # ensure only CSF is measured
    free_water *= (free_water >= free_water_thres)

    one_nib = nib.load(root+"/FW.nii.gz")
    header = one_nib.header["pixdim"]
    csf_label = fort.parcellate(free_water, seg_gray.numpy(),
                                header[1], header[2], header[3], cutoff_dist)

    csf_label = nib.Nifti1Image(csf_label, affine=one_nib.affine)
    csf_label.to_filename(root+"/CSF_parcel.nii.gz")


def parcellate_vent(root, free_water_thres=0.8, cutoff_dist=10):
    seg_out = ants.image_read(root+"/seg.nii.gz")
    vent = keep_vent(seg_out)
    vent_mask = ants.threshold_image(vent, low_thresh=1)

    free_water = ants.image_read(root+"/FW.nii.gz")
    mask = ants.image_read(root+"/FW_synth_mask.nii.gz")
    free_water = ants.mask_image(free_water, mask)
    free_water = ants.mask_image(free_water, vent_mask)
    csf_mask = ants.threshold_image(
        free_water, low_thresh=free_water_thres, binary=True)

    one_nib = nib.load(root+"/FW.nii.gz")
    header = one_nib.header["pixdim"]
    vent_label = fort.parcellate(csf_mask.numpy(), vent.numpy(),
                                 header[1], header[2], header[3], cutoff_dist)
    vent_label = nib.Nifti1Image(vent_label, affine=one_nib.affine)
    vent_label.to_filename(root+"/vent_parcel.nii.gz")


# --------------------Sulcus atlas map--------------------
def custom_mask(arr, labels):
    out = np.zeros(arr.shape)
    for i in labels:
        out += arr.astype(int) == i
    return out > 0


def fill_gray_for_sulcus(atlas_path, seg_path, save_path):
    atlas_nib = nib.load(atlas_path)
    atlas = atlas_nib.get_fdata()
    seg = nib.load(seg_path).get_fdata()
    all_labels = [
        3, 42, 7, 8, 46, 47, 16, 17, 53, 18, 54, 26, 28, 58, 60,
        1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010,
        1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020,
        1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030,
        1031, 1032, 1033, 1034, 1035, 2001, 2002, 2003, 2004, 2005,
        2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015,
        2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025,
        2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035
        ]
    gray = custom_mask(seg, all_labels)
    pix = atlas_nib.header["pixdim"]
    out_mask = fort.fill_gray(atlas, gray, pix[1], pix[2], pix[3], 3)
    out_mask = nib.Nifti1Image(out_mask, affine=atlas_nib.affine)
    out_mask.to_filename(save_path)


def create_sr(pT2_path, mask_path, tmp_dir):
    # create pseudo-T1 for registration
    if not os.path.exists(tmp_dir+"/pT1_reg.nii.gz"):
        print("creating high resolution pseudo-T1")
        if not os.path.exists(tmp_dir+"/pT1.nii.gz"):
            bash_in_python(
                "mri_synthsr --i {}".format(pT2_path) +
                " --o {}/pT1.nii.gz".format(tmp_dir) +
                " --cpu")

        pT2 = ants.image_read(pT2_path)
        pT1 = ants.image_read(tmp_dir+"/pT1.nii.gz")
        tx = ants.registration(pT2, pT1, type_of_transform="Rigid",
                               outprefix=tmp_dir+"/ants")
        mask = ants.image_read(mask_path)
        ants.image_write(ants.mask_image(tx["warpedmovout"], mask),
                         tmp_dir+"/pT1_reg.nii.gz")


def register_sulcus(pT2_path, seg_path, sulcus_atlas_path,
                    MNI_template_path, FW_path, FW_mask_path,
                    save_dir, free_water_thres=0.8,
                    cutoff_dist=10):
    '''
    Method:
    1. Register sulcus atlas map from MNI space to the pseudo-T2 space
    2. Use freesurfer to segment pseudo-T2 (already done)
    3. Extract GM from the freesurfer mask
    4. Fill the GM with the sulcus atlas map
    5. Use the GM labels to identify sulci
    '''
    tmp_dir = save_dir+"/tmp"
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    # create_seg(pT2_path, seg_path)
    # create_mask(FW_path, FW_mask_path)
    create_sr(pT2_path, FW_mask_path, tmp_dir)

    # register atlas to pT2
    MNI = ants.image_read(MNI_template_path)
    atlas = ants.image_read(sulcus_atlas_path)

    pT1 = ants.image_read(tmp_dir+"/pT1_reg.nii.gz")
    tx = ants.registration(pT1, MNI, type_of_transform="SyN",
                           outprefix=tmp_dir+"/ants")
    atlas_pt2 = ants.apply_transforms(pT1, atlas, tx["fwdtransforms"],
                                      interpolator="genericLabel")

    # fill the gray matter with the atlas
    ants.image_write(atlas_pt2, tmp_dir+"/sulcus_atlas_reg.nii.gz")
    fill_gray_for_sulcus(tmp_dir+"/sulcus_atlas_reg.nii.gz", seg_path,
                         tmp_dir+"/sulcus_atlas_filled.nii.gz")

    seg_out = ants.image_read(seg_path)
    seg_gray = ants.image_read(tmp_dir+"/sulcus_atlas_filled.nii.gz")
    no_csf = remove_vent(seg_out)

    free_water = ants.image_read(FW_path)
    mask = ants.image_read(FW_mask_path)
    # apply skull mask
    free_water = ants.mask_image(free_water, mask)
    # remove ventricle
    free_water = ants.mask_image(free_water, no_csf)
    # remove lacunes
    free_water = free_water.numpy()*(seg_out.numpy() == 24)
    # ensure only CSF is measured
    free_water *= (free_water >= free_water_thres)

    one_nib = nib.load(FW_path)
    header = one_nib.header["pixdim"]
    csf_label = fort.parcellate(free_water, seg_gray.numpy(),
                                header[1], header[2], header[3], cutoff_dist)

    csf_label = nib.Nifti1Image(csf_label, affine=one_nib.affine)
    csf_label.to_filename(save_dir+"/sulcus_cistern.nii.gz")


# ----------final pipeline----------
def parcellate(wdir, num_cpu=1, cutoff_dist=10, free_water_thres=0.8):
    '''
    1. perform cortical parcellation using mri_synthseg --parc
    2. synthesize MRI to do registration
    3. register cortical parcellation to pT1

    create 2 more files in `wdir`: `CSF_parcel.nii.gz`, `vent_parcel.nii.gz`
    '''
    if not os.path.exists(wdir+"/pT2_seg.nii.gz"):
        bash_in_python(
            "mri_synthseg --i {}".format(wdir+"/pT2.nii.gz") +
            " --o {}".format(wdir+"/pT2_seg.nii.gz") +
            " --resample {}".format(wdir+"/pT1_mri.nii.gz") +
            " --parc --robust --threads {}".format(num_cpu)
            )

    if not os.path.exists(wdir+"/CSF_parcel.nii.gz"):
        parcellate_csf(wdir, cutoff_dist=10, free_water_thres=0.8)

    if not os.path.exists(wdir+"/vent_parcel.nii.gz"):
        parcellate_vent(wdir, cutoff_dist=10, free_water_thres=0.8)

    if not os.path.exists(wdir+"/sulcus_cistern.nii.gz"):
        MNI_template_path = os.environ["FSLDIR"] + \
            "/data/standard/MNI152_T1_1mm_brain.nii.gz"
        sulcus_atlas_path = os.path.dirname(os.path.realpath(__file__)) + \
            "/../data/sulcus_corr.nii.gz"
        register_sulcus(
            wdir+"/pT2.nii.gz", wdir+"/seg.nii.gz",
            sulcus_atlas_path, MNI_template_path,
            wdir+"/FW.nii.gz", wdir+"/FW_synth_mask.nii.gz", wdir)
