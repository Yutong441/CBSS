'''
Analyse perfusion related data

* ASL
* IVIM: one-step biexponential fitting
* fMRI
1. slice timing correction: slicetimer
2. motion correction: mcflirt -in {} -out {}
3. Gaussian filtering
4. skullstrip
5. linear detrending
'''
import re
import os
import subprocess
import json
import shutil
import glob

import numpy as np
import pandas as pd
import nibabel as nib
import scipy
import ants
from src import fort


def bash_in_python(cmd):
    out, err = unix_cmd(cmd)
    show_error(err)
    return out, err


def show_error(err):
    if len(err) > 0:
        print(err)


def unix_cmd(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
    out, err = p.communicate()
    return out, err


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


def custom_mask(arr, labels):
    out = np.zeros(arr.shape)
    for i in labels:
        out += arr == i
    return out > 0


def get_tcustom(json_path, save_path):
    '''Obtain the tcustom file for slicetimer'''
    with open(json_path, "r") as f:
        jfile = json.load(f)

    TR = jfile["RepetitionTime"]
    slicetimes = np.array(jfile["SliceTiming"])
    np.savetxt(save_path, slicetimes/TR, fmt='%.4f')
    return TR


def get_design(img_path, save_path):
    img4d = nib.load(img_path)
    N = img4d.shape[-1]
    design_matrix = np.column_stack((np.arange(1, N + 1),
                                     np.ones(N)))
    np.savetxt(save_path, design_matrix, fmt="%.4f",
               delimiter="\t")


def filter_signal(data, tr, lowcut=0.02, highcut=0.04, order=9):
    nyq = (1 / tr) / 2
    low = lowcut / nyq
    high = highcut / nyq
    a, b = scipy.signal.butter(int(order), [low, high], btype="band")
    filt_data = scipy.signal.filtfilt(a, b, data, axis=-1)
    return filt_data


def filter_fit(fmri4d_path, mask_path, global_bold_mask_path,
               baseline_path, save_path,
               TR, lag_max, lowcut, highcut):
    '''
    Args:
        `fmri4d_path`: preprocessed fMRI data after motion correction, slice
        time correction, Gaussian filtering, linear detrending
        `mask_path`: mask to obtain the delay/CVR map
    '''
    img4d_nib = nib.load(fmri4d_path)
    img4d = img4d_nib.get_fdata()
    img4d_fil = filter_signal(img4d, TR, lowcut=lowcut, highcut=highcut)
    baseline = nib.load(baseline_path).get_fdata()

    mask = nib.load(mask_path).get_fdata()
    global_mask = nib.load(global_bold_mask_path).get_fdata()
    # reconstruct original BOLD signals: BOLD0 + Delta BOLD
    global_bold = (global_mask[..., None]*(img4d_fil + baseline[..., None])
                   ).sum(axis=(0, 1, 2))/global_mask.sum(axis=(0, 1, 2))

    max_delay = int(np.ceil(lag_max/TR))
    regress_y = img4d_fil/baseline[..., None]
    # if a voxel has an intensity of 0, make sure it is 0 in the y variable
    regress_y *= (baseline[..., None] > 0).astype(float)
    maps = fort.lag_regression4d(
        np.moveaxis(regress_y, 3, 0),
        mask, global_bold, max_delay)

    # find out global CVR
    nodelay = maps[..., 1][global_mask == 1]
    nodelay = np.mean(nodelay[nodelay > 0])
    delay = maps[..., 2][global_mask == 1]
    delay = np.mean(delay[delay > 0])

    # normalize by global CVR values
    maps[..., 1] /= nodelay
    maps[..., 2] /= delay

    maps = nib.Nifti1Image(maps, affine=img4d_nib.affine)
    maps.to_filename(save_path)


def fmri_preproc(img_path, wdir, json_path, save_dir, fwmh=8, lag_max=8):
    tmp_dir = save_dir+"/tmp"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    TR = get_tcustom(json_path, tmp_dir+"/slicetime.txt")
    sigma = fwmh/2.355
    if not os.path.exists(tmp_dir+"/gauss.nii.gz"):
        # motion correction, slice timing, Gaussian filtering, skullstripping
        bash_in_python(
            "mcflirt -in {}".format(img_path) +
            " -out {}/mc.nii.gz;".format(tmp_dir) +
            " slicetimer -i {}/mc.nii.gz".format(tmp_dir) +
            " -o {}/slicetime.nii.gz".format(tmp_dir) +
            " --tcustom={}/slicetime.txt;".format(tmp_dir) +
            " fslmaths {}/slicetime.nii.gz -s {}".format(tmp_dir, sigma) +
            " {}/gauss.nii.gz;".format(tmp_dir) +
            " fslselectvols -i {}/slicetime.nii.gz".format(tmp_dir) +
            " -o {}/f0.nii.gz --vols=0;".format(tmp_dir) +
            " bet {}/f0.nii.gz".format(tmp_dir) +
            " {}/f0_brain.nii.gz -m".format(tmp_dir)
        )

    if not os.path.exists(tmp_dir+"/detrend.nii.gz"):
        img4d_nib = nib.load(tmp_dir+"/gauss.nii.gz")
        img4d = img4d_nib.get_fdata()
        mask = nib.load(tmp_dir+"/f0_brain_mask.nii.gz").get_fdata()
        img4d_de = fort.detrend4d(np.moveaxis(img4d, 3, 0), mask)
        img4d_de = np.moveaxis(img4d_de, 0, 3)
        img4d_de_nib = nib.Nifti1Image(img4d_de[..., 1:],
                                       affine=img4d_nib.affine)
        img4d_de_nib.to_filename(tmp_dir+"/detrend.nii.gz")
        baseline = nib.Nifti1Image(img4d_de[..., 0],
                                   affine=img4d_nib.affine)
        baseline.to_filename(tmp_dir+"/baseline.nii.gz")

    f0 = ants.image_read(tmp_dir+"/f0_brain.nii.gz")
    warp_to_modality(f0, wdir, tmp_dir, tmp_dir+"/reg_fmri")
    shutil.copy2(tmp_dir+"/reg_fmri_gm.nii.gz", save_dir)
    shutil.copy2(tmp_dir+"/reg_fmri_sulcus.nii.gz", save_dir)

    if not os.path.exists(tmp_dir+"/maps_2_4.nii.gz"):
        # two sets of filtering
        # 0.01-0.1Hz: from the Amemiya 2014 study
        # 0.02-0.04Hz: from Liu 2017 study
        mask_path = tmp_dir+"/f0_brain_mask.nii.gz"
        filter_fit(tmp_dir+"/detrend.nii.gz", mask_path, mask_path,
                   tmp_dir+"/baseline.nii.gz",
                   tmp_dir+"/maps_1_1.nii.gz",
                   TR, lag_max, 0.01, 0.1)

        # obtain gray matter mask
        gm_mask_path = tmp_dir+"/reg_fmri_gm.nii.gz"
        filter_fit(tmp_dir+"/detrend.nii.gz", mask_path, gm_mask_path,
                   tmp_dir+"/baseline.nii.gz",
                   tmp_dir+"/maps_2_4.nii.gz",
                   TR, lag_max, 0.02, 0.04)

    img_nib = nib.load(tmp_dir+"/maps_1_1.nii.gz")
    img = img_nib.get_fdata()
    delay_cvr = nib.Nifti1Image(img[..., 0], affine=img_nib.affine)
    delay_cvr.to_filename(save_dir+"/fMRI_delay.nii.gz")

    img = nib.load(tmp_dir+"/maps_2_4.nii.gz").get_fdata()
    cvr = nib.Nifti1Image(img[..., 1], affine=img_nib.affine)
    cvr.to_filename(save_dir+"/fMRI_CVR.nii.gz")
    cvr = nib.Nifti1Image(img[..., 2], affine=img_nib.affine)
    cvr.to_filename(save_dir+"/fMRI_CVRlag.nii.gz")
    # shutil.copy2(tmp_dir+"/f0_brain.nii.gz", save_dir)


def ivim_map(dwi_fname, bval_fname, mask_fname, ofolder, bval_thres=2000):
    '''Use dipy to fit IVIM data'''
    mask = nib.load(mask_fname).get_fdata()
    bval = np.genfromtxt(bval_fname)
    ivim_nib = nib.load(dwi_fname)
    ivim = ivim_nib.get_fdata()
    maps = fort.biexp_fit_all4d(
        np.moveaxis(ivim, 3, 0)[bval <= bval_thres],
        bval[bval <= bval_thres], mask)

    Dstar = nib.Nifti1Image(maps[..., 1], affine=ivim_nib.affine)
    Dstar.to_filename(ofolder+"/Dstar_map.nii.gz")
    Fivim = nib.Nifti1Image(maps[..., 2], affine=ivim_nib.affine)
    Fivim.to_filename(ofolder+"/Fivim_map.nii.gz")

    # ivim = ivim[..., bval <= bval_thres]
    # bval = bval[bval <= bval_thres]
    # bvecs = np.zeros([3, len(bval)])
    # bvecs[0] = 1
    # gtab = gradient_table(bval, bvecs, b0_threshold=0)

    # ivimmodel = IvimModel(gtab, fit_method='trr')
    # ivimfit = ivimmodel.fit(ivim*mask[..., None])

    # Fivim = ivimfit.model_params[..., 1]
    # Fivim = nib.Nifti1Image(Fivim, affine=ivim_nib.affine)
    # Fivim.to_filename(ofolder+"/Fivim.nii.gz")

    # Dstar = ivimfit.model_params[..., 2]
    # Dstar = nib.Nifti1Image(Dstar, affine=ivim_nib.affine)
    # Dstar.to_filename(ofolder+"/Dstar.nii.gz")


def regional_val(one_map, atlas_mask, save_path, thres=20):
    out_dict = {}
    all_regions = np.unique(atlas_mask)
    all_regions = all_regions[all_regions != 0]

    for i in all_regions:
        diff_val = one_map[atlas_mask == i]
        # remove the values below 0
        diff_val = diff_val[diff_val > 0]
        diff_val = diff_val[np.isfinite(diff_val)]

        if len(diff_val) > thres:
            out_dict[i] = {
                "mean_val": diff_val.mean(),
                "median_val": np.median(diff_val),
                "max_val": diff_val.max(),
                "min_val": diff_val.min(),
                "sd": diff_val.std(),
                "vol": len(diff_val)
            }

    out_dict = pd.DataFrame.from_dict(out_dict, orient="index")
    out_dict.to_csv(save_path)


def register_gray(pT2_path, seg_path, sulcus_atlas_path,
                  MNI_template_path, FW_mask_path, tmp_dir):
    '''Register sulcus atlas to pseudo-T2'''
    # create_sr(pT2_path, FW_mask_path, tmp_dir)

    # register atlas to pT2
    MNI = ants.image_read(MNI_template_path)
    atlas = ants.image_read(sulcus_atlas_path)

    pT2 = ants.image_read(pT2_path)
    tx = ants.registration(pT2, MNI, type_of_transform="SyN",
                           outprefix=tmp_dir+"/ants")
    atlas_pt2 = ants.apply_transforms(pT2, atlas, tx["fwdtransforms"],
                                      interpolator="genericLabel")

    # fill the gray matter with the atlas
    ants.image_write(atlas_pt2, tmp_dir+"/sulcus_atlas_reg.nii.gz")
    fill_gray_for_sulcus(tmp_dir+"/sulcus_atlas_reg.nii.gz", seg_path,
                         tmp_dir+"/sulcus_atlas_filled.nii.gz")


def keep_cortex(img):
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
    all_labels = np.unique(gm)
    all_labels = all_labels[all_labels < 1000]
    for i in [3, 42]:
        all_labels = all_labels[all_labels != i]
    for i in all_labels:
        gm *= (gm != i).astype(int)
    gm = ants.from_numpy(gm, origin=img.origin, direction=img.direction,
                         spacing=img.spacing)
    return gm


def warp_to_modality(img, wdir, tmp_dir, prefix):
    dwi = ants.image_read(wdir+"/DTI/b0.nii.gz")
    tx = ants.registration(img, dwi, outprefix=prefix,
                           type_of_transform="Rigid")

    MNI_template_path = os.environ["FSLDIR"] + \
        "/data/standard/MNI152_T1_1mm_brain.nii.gz"
    sulcus_atlas_path = os.path.dirname(os.path.realpath(__file__)) + \
        "/../data/sulcus_corr.nii.gz"
    register_gray(
        wdir+"/DTI/pT2.nii.gz",
        wdir+"/DTI/seg.nii.gz",
        sulcus_atlas_path,
        MNI_template_path,
        wdir+"/DTI/FW_synth_mask.nii.gz",
        tmp_dir)

    sulcus_atlas = ants.image_read(tmp_dir+"/sulcus_atlas_filled.nii.gz")
    gm_atlas = ants.image_read(wdir+"/CBSS/DTI/seg.nii.gz")

    sulcus_atlas_reg = ants.apply_transforms(
        img, sulcus_atlas, tx["fwdtransforms"],
        interpolator="genericLabel")

    gm_atlas_reg = ants.apply_transforms(
        img, gm_atlas, tx["fwdtransforms"],
        interpolator="genericLabel")
    gm_atlas_reg = keep_cortex(gm_atlas_reg)
    gm = ants.threshold_image(gm_atlas_reg, low_thresh=0.5,
                              binary=True)

    ants.image_write(sulcus_atlas_reg, prefix+"_sulcus.nii.gz")
    ants.image_write(gm, prefix+"_gm.nii.gz")


def region_val_atlas(map_dict, prefix, tmp_dir, thres=20):
    sulcus_atlas_reg = ants.image_read(prefix+"_sulcus.nii.gz")
    gm = ants.image_read(prefix+"_gm.nii.gz")

    for key, val in map_dict.items():
        one_map = ants.image_read(val)
        regional_val(one_map.numpy(), gm.numpy(),
                     tmp_dir+"/stats_gm_{}.csv".format(key),
                     thres=thres)

        regional_val(one_map.numpy(), sulcus_atlas_reg.numpy(),
                     tmp_dir+"/stats_sulcus_{}.csv".format(key),
                     thres=thres)


def ivim_pipeline(wdir, bval_fname, tmp_dir):
    if not os.path.exists(tmp_dir+"/b0_mask.nii.gz"):
        bash_in_python(
            "mri_synthstrip -i {}/CBSS/b0_ivim.nii.gz".format(wdir) +
            " -o {}/b0_brain.nii.gz".format(tmp_dir) +
            " -m {}/b0_mask.nii.gz".format(tmp_dir)
            )
    if not os.path.exists(tmp_dir+"/Dstar_map.nii.gz"):
        ivim_map(wdir+"/DWI/corr_ivim.nii.gz", bval_fname,
                 tmp_dir+"/b0_mask.nii.gz", tmp_dir)

    # warp pT2 to IVIM
    ivim = ants.image_read(wdir+"/CBSS/b0_ivim.nii.gz")
    mask = ants.image_read(tmp_dir+"/b0_mask.nii.gz")
    ivim = ants.mask_image(ivim, mask)
    warp_to_modality(ivim, wdir, tmp_dir, tmp_dir+"/reg_ivim")
    map_dict = {"Dstar": tmp_dir+"/Dstar_map.nii.gz",
                "Fivim": tmp_dir+"/Fivim_map.nii.gz"}
    region_val_atlas(map_dict, tmp_dir+"/reg_ivim", tmp_dir)


def fMRI_pipeline(fmri_path, fmri_json_path, wdir, tmp_dir):
    fmri_preproc(fmri_path, wdir, fmri_json_path,
                 tmp_dir, fwmh=8, lag_max=8)
    map_dict = {"delay": tmp_dir+"/fMRI_delay.nii.gz",
                "CVR": tmp_dir+"/fMRI_CVR.nii.gz",
                "CVRlag": tmp_dir+"/fMRI_CVRlag.nii.gz",
                }
    region_val_atlas(map_dict, tmp_dir+"/reg_fmri", tmp_dir, thres=5)


def make_sr(brain_path, tmp_dir, mask_path=None):
    print("creating high resolution pseudo-T1")
    if mask_path is None:
        bash_in_python(
            "mri_synthstrip -i {}".format(brain_path) +
            " -o {}/brain.nii.gz".format(tmp_dir) +
            " -m {}/brain_mask.nii.gz".format(tmp_dir)
        )
        mask_path = tmp_dir+"/brain_mask.nii.gz"

    if not os.path.exists(tmp_dir+"/pT1.nii.gz"):
        bash_in_python(
            "mri_synthsr --i {}".format(brain_path) +
            " --o {}/pT1.nii.gz".format(tmp_dir) +
            " --cpu")

    pT2 = ants.image_read(brain_path)
    pT1 = ants.image_read(tmp_dir+"/pT1.nii.gz")
    tx = ants.registration(pT2, pT1, type_of_transform="Rigid",
                           outprefix=tmp_dir+"/ants")
    mask = ants.image_read(mask_path)
    out_img = ants.mask_image(tx["warpedmovout"], mask)
    ants.image_write(out_img, tmp_dir+"/pT1_reg.nii.gz")
    return out_img


def ASL_pipeline(asl_path, wdir, tmp_dir):
    CBF_nib = nib.load(asl_path)
    CBF = CBF_nib.get_fdata()
    if len(CBF.shape) == 4:
        CBF = CBF[..., 0]
    CBF = nib.Nifti1Image(CBF, affine=CBF_nib.affine)
    CBF.to_filename(tmp_dir+"/CBF.nii.gz")

    CBF = ants.image_read(tmp_dir+"/CBF.nii.gz")
    warp_to_modality(CBF, wdir, tmp_dir, tmp_dir+"/reg_asl")
    map_dict = {"CBF": tmp_dir+"/CBF.nii.gz"}
    region_val_atlas(map_dict, tmp_dir+"/reg_asl", tmp_dir, thres=10)


def clean_up(wdir):
    blood_dir = wdir+"/CBSS/blood"
    stats_dir = wdir+"/CBSS/blood_stats"
    if not os.path.exists(blood_dir):
        os.mkdir(blood_dir)
        os.mkdir(stats_dir)

    for i in ["Dstar_map", "Fivim_map", "CBF", "fMRI_delay",
              "fMRI_CVR", "fMRI_CVRlag"]:
        need_file = wdir+"/CBSS/tmp/{}.nii.gz".format(i)
        if os.path.exists(need_file):
            shutil.copy2(need_file, blood_dir)

    for i in glob.glob(wdir+"/CBSS/tmp/stats*.csv"):
        shutil.copy2(i, stats_dir)

    shutil.rmtree(wdir+"/CBSS/tmp")


def run_one_sample(ivim_path, ivim_bval_path, fmri_path,
                   fmri_json_path, asl_path, save_dir):
    '''
    Method:
    1. extract parameter maps for ASL, IVIM, fMRI
    2. register pseudo-T2 to ASL and fMRI
    3. warp the sulcus atlas mask to ASL and fMRI
    4. repeat the same procedure for IVIM
    '''
    tmp_dir = save_dir+"/CBSS/tmp"
    if os.path.exists(save_dir+"/CBSS"):
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)

    ivim_pipeline(save_dir, ivim_bval_path, tmp_dir)
    if os.path.exists(fmri_path):
        fMRI_pipeline(fmri_path, fmri_json_path, save_dir, tmp_dir)

    if os.path.exists(asl_path):
        ASL_pipeline(asl_path, save_dir, tmp_dir)

    clean_up(save_dir)


def ivim_config(config: dict, root: str, ID: str):
    cf = {}
    all_files = ["IVIM", "IVIM_bval", "fMRI", "fMRI_json", "ASL"]
    proceed = len(all_files)
    for i in all_files:
        cf[i] = _add_ID([config[i]], root, ID)[0]
        if os.path.exists(cf[i]):
            proceed -= 1

    cf["save"] = root+"/"+config["save"]+"/"+ID
    return cf, proceed


def modify_str(txt, ID):
    ''' translate regular expression in the data_regex.json file '''
    if "*8IDIDIDID" in txt:
        numer = re.sub("([a-z]|[A-Z])+", "", ID)
        ID_num = numer.zfill(8)
        return re.sub("8IDIDIDID", ID_num, txt)
    elif "IDIDIDID" in txt:
        return re.sub("IDIDIDID", ID, txt)


def _add_ID(img_list, root, ID):
    ilist = [root+"/"+modify_str(i, ID) for i in img_list]
    # account for regular expressions in json file
    for index, i in enumerate(ilist):
        if "*" in i:
            out = sorted(glob.glob(i))
            if len(out) > 0:
                ilist[index] = out[0]
    return ilist


def pipeline_all(root, data_name, num=0, arrayID=0):
    cf_name = os.path.dirname(os.path.realpath(__file__))+"/data_regex.json"
    with open(cf_name, "r") as f:
        cf = json.load(f)
    one_cf = cf[data_name]
    all_IDs = sorted(os.listdir(root+"/"+one_cf["ID"]))

    for index, i in enumerate(all_IDs):
        if index % num == arrayID:
            config, proceed = ivim_config(one_cf, root, i)
            if proceed <= 2:
                print(i)
                if not os.path.exists(config["save"]):
                    os.mkdir(config["save"])

                save_goal = config["save"]+"/CBSS/blood_stats/stats_gm_CVR.csv"
                if not os.path.exists(save_goal):
                    run_one_sample(
                        config["IVIM"], config["IVIM_bval"],
                        config["fMRI"], config["fMRI_json"],
                        config["ASL"], config["save"])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_name', type=str)
    parser.add_argument('--root', type=str)
    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()
    pipeline_all(args.root, args.data_name, num=args.num,
                 arrayID=args.arrayID)
