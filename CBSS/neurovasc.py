'''
1. Register sulcus atlas to the DTI space
2. Register ALFF and CBF maps to the DTI space
3. Neurovascular coupling
4. Coupling and ALFF in each region
'''
import os
import re
import glob
import json
import shutil
import subprocess
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


def get_tcustom(json_path, save_path):
    '''Obtain the tcustom file for slicetimer'''
    with open(json_path, "r") as f:
        jfile = json.load(f)

    TR = jfile["RepetitionTime"]
    slicetimes = np.array(jfile["SliceTiming"])
    np.savetxt(save_path, slicetimes/TR, fmt='%.4f')
    return TR


def custom_mask(arr, labels):
    out = np.zeros(arr.shape)
    for i in labels:
        out += arr == i
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
    '''
    Args:
        `img`: ants image object
        `wdir`: working directory for the patient
        `prefix`: ants image registration prefix
    '''
    dwi = ants.image_read(wdir+"/DTI/pT2.nii.gz")
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

    # register sulcus atlas to the image
    tx = ants.registration(img, dwi, outprefix=prefix,
                           type_of_transform="Rigid")
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


def detrend(inp_path, mask_path, out_path, method="quadratic"):
    img4d_nib = nib.load(inp_path)
    img4d = img4d_nib.get_fdata()
    mask = nib.load(mask_path).get_fdata()

    if method == "linear":
        img4d_de = fort.detrend4d(np.moveaxis(img4d, 3, 0), mask)
    else:
        N = img4d.shape[3]
        confound = np.stack(
            [np.arange(N), (np.arange(N))**2],
            axis=-1)
        img4d_de = fort.detrend_quad(np.moveaxis(img4d, 3, 0), mask,
                                     confound)

    img4d_de = np.moveaxis(img4d_de, 0, 3)
    img4d_de_nib = nib.Nifti1Image(img4d_de,
                                   affine=img4d_nib.affine)
    img4d_de_nib.to_filename(out_path)


def remove_outliers(img4d_nib, mask, confound_mat):
    img4d = img4d_nib.get_fdata()
    img4d = np.moveaxis(img4d, 3, 0)

    img4d_de = fort.detrend_quad(img4d, mask, confound_mat)
    img4d_de = np.moveaxis(img4d_de, 0, 3)
    img4d_de_nib = nib.Nifti1Image(img4d_de,
                                   affine=img4d_nib.affine)
    return img4d_de_nib


def filter_signal(data, tr, lowcut=0.02, highcut=0.04, order=9):
    nyq = (1 / tr) / 2
    low = lowcut / nyq
    high = highcut / nyq
    a, b = scipy.signal.butter(int(order), [low, high], btype="band")
    filt_data = scipy.signal.filtfilt(a, b, data, axis=-1)
    return filt_data


def filter_img4d(inp_path, out_path, TR, lowcut, highcut):
    img4d_nib = nib.load(inp_path)
    img4d = img4d_nib.get_fdata()
    baseline = img4d.mean(axis=3)[..., None]
    img4d_fil = filter_signal(img4d, TR, lowcut=lowcut, highcut=highcut)
    img4d_fil = nib.Nifti1Image(img4d_fil + baseline,
                                affine=img4d_nib.affine)
    img4d_fil.to_filename(out_path)
    del img4d_fil, img4d, img4d_nib


def remove_motion_outliers(img_list, mask_path, save_dir):
    bash_in_python(
        "fsl_motion_outliers -i {}".format(img_list[0]) +
        " -o {}/confound.txt -nomoco".format(save_dir) +
        " -m {}".format(mask_path)
            )

    for i in img_list:
        param = np.genfromtxt(save_dir+"/regress_param.txt")
        if os.path.exists(save_dir+"/confound.txt"):
            motion = np.genfromtxt(save_dir+"/confound.txt")
            param = np.concatenate([param, motion], axis=1)

        mask = nib.load(mask_path).get_fdata()
        img4d_nib = nib.load(i)
        img4d_de_nib = remove_outliers(img4d_nib, mask, param)
        img4d_de_nib.to_filename(save_dir+"/"+re.sub(
            ".nii.gz", "_demot.nii.gz", os.path.basename(i)))


def obtain_regress_param(mc_path, f0_path, seg_img_path, seg_path, save_dir):
    mc = ants.image_read(mc_path).numpy()
    f0 = ants.image_read(f0_path)
    seg_img = ants.image_read(seg_img_path)
    seg = ants.image_read(seg_path)

    tx = ants.registration(f0, seg_img, type_of_transform="Rigid",
                           outprefix=save_dir+"/f0")
    seg_reg = ants.apply_transforms(f0, seg, tx["fwdtransforms"],
                                    interpolator="genericLabel")
    CSF_mask = custom_mask(seg_reg.numpy(), [24])
    WM_mask = custom_mask(seg_reg.numpy(), [2, 41])
    WM_sig = np.mean(mc*WM_mask[..., None], axis=(0, 1, 2))
    CSF_sig = np.mean(mc*CSF_mask[..., None], axis=(0, 1, 2))

    motion = np.genfromtxt(save_dir+"/mc.nii.gz.par")
    motion = np.concatenate([motion, WM_sig[..., None], CSF_sig[..., None]],
                            axis=-1)
    np.savetxt(save_dir+"/regress_param.txt", motion)


def remove_initial_slices(img_path, save_path, n_initial=10):
    img_nib = nib.load(img_path)
    img = img_nib.get_fdata()
    img = img[..., n_initial:]
    nib.Nifti1Image(img, affine=img_nib.affine).to_filename(save_path)


def fmri_preproc(img_path, json_path, save_dir, fwmh=6, lag_max=8):
    tmp_dir = save_dir+"/tmp"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    TR = get_tcustom(json_path, tmp_dir+"/slicetime.txt")
    sigma = fwmh/2.355
    if not os.path.exists(tmp_dir+"/gauss.nii.gz"):
        # motion correction, slice timing, Gaussian filtering, skullstripping
        remove_initial_slices(img_path, tmp_dir+"/initial_f.nii.gz")
        bash_in_python(
            "mcflirt -in {}/initial_f.nii.gz".format(tmp_dir) +
            " -out {}/mc.nii.gz".format(tmp_dir) +
            " -plots;" +
            " slicetimer -i {}/mc.nii.gz".format(tmp_dir) +
            " -o {}/slicetime.nii.gz".format(tmp_dir) +
            " --tcustom={}/slicetime.txt;".format(tmp_dir) +
            " fslselectvols -i {}/slicetime.nii.gz".format(tmp_dir) +
            " -o {}/f0.nii.gz --vols=0;".format(tmp_dir) +
            " bet {}/f0.nii.gz".format(tmp_dir) +
            " {}/f0_brain.nii.gz -m;".format(tmp_dir) +
            # apply brain mask
            " fslmaths {}/mc.nii.gz".format(tmp_dir) +
            " -mul {}/f0_brain_mask.nii.gz".format(tmp_dir) +
            " {}/mc_ss.nii.gz;".format(tmp_dir) +
            # Gaussian smoothening
            " fslmaths {}/mc_ss.nii.gz -s {}".format(tmp_dir, sigma) +
            " {}/gauss.nii.gz".format(tmp_dir)
        )

    if not os.path.exists(tmp_dir+"/alff.nii.gz"):
        detrend(tmp_dir+"/gauss.nii.gz",
                tmp_dir+"/f0_brain_mask.nii.gz",
                tmp_dir+"/detrend.nii.gz", method="quadratic")

        obtain_regress_param(tmp_dir+"/mc_ss.nii.gz",
                             tmp_dir+"/f0_brain.nii.gz",
                             save_dir+"/DTI/pT2.nii.gz",
                             save_dir+"/DTI/seg.nii.gz",
                             tmp_dir)

        remove_motion_outliers(
            [tmp_dir+"/detrend.nii.gz"],
            tmp_dir+"/f0_brain_mask.nii.gz", tmp_dir)

        filter_img4d(
            tmp_dir+"/detrend_demot.nii.gz", tmp_dir+"/temp_fil.nii.gz",
            TR, 0.01, 0.1)

        cal_alff(tmp_dir+"/temp_fil.nii.gz",
                 tmp_dir+"/f0_brain_mask.nii.gz",
                 tmp_dir, TR)


def cal_alff(ori_path, mask_path, save_dir, TR):
    mask = nib.load(mask_path).get_fdata()
    img4d_nib = nib.load(ori_path)
    img4d = img4d_nib.get_fdata()
    # img4d -= img4d.mean(axis=3, keepdims=True)
    img4d_fft = np.fft.fft(img4d, axis=-1)
    freqs = np.fft.fftfreq(img4d.shape[-1], TR)
    band = np.abs(img4d_fft[..., (freqs >= 0.01)*(freqs <= 0.1)])
    alff = np.mean(band, axis=-1)

    # global_mean = np.sum(alff*mask)/np.sum(mask)
    alff_nib = nib.Nifti1Image(alff,
                               affine=img4d_nib.affine)
    alff_nib.to_filename(save_dir+"/alff.nii.gz")


def neurovasc_cor(img1, img2, atlas, thres=20):
    out_dict = {}
    all_regions = np.unique(atlas)
    all_regions = all_regions[all_regions != 0]

    for i in all_regions:
        val1 = img1[atlas == i]
        val2 = img2[atlas == i]

        val1 = val1[val2 > 0]
        val2 = val2[val2 > 0]
        val2 = val2[val1 > 0]
        val1 = val1[val1 > 0]

        val1 = val1[np.isfinite(val2)]
        val2 = val2[np.isfinite(val2)]
        val2 = val2[np.isfinite(val1)]
        val1 = val1[np.isfinite(val1)]

        if len(val1) > thres:
            out_dict[i] = {
                "corr": np.corrcoef(val1, val2)[0, 1],
                "ratio": np.median(val2/(val1 + 1)),
                "vol": len(val1)
            }

    out_dict = pd.DataFrame.from_dict(out_dict, orient="index")
    return out_dict


def get_neurovasc(cbf_path, f0_path, b0_path, sulcus_atlas_path,
                  tmp_dir):
    # register fMRI and ASL to DTI
    CBF = ants.image_read(cbf_path)
    ALFF = ants.image_read(tmp_dir+"/alff.nii.gz")
    f0 = ants.image_read(f0_path)
    b0 = ants.image_read(b0_path)

    tx = ants.registration(b0, f0, type_of_transform="Rigid",
                           outprefix=tmp_dir+"/f0")
    alff_reg = ants.apply_transforms(b0, ALFF, tx["fwdtransforms"])
    cbf_reg = ants.registration(b0, CBF, type_of_transform="Rigid",
                                outprefix=tmp_dir+"/cbf")

    # for each sulcus, get neurovascular coupling
    sulcus_atlas = ants.image_read(sulcus_atlas_path)
    NVC = neurovasc_cor(alff_reg.numpy(), cbf_reg["warpedmovout"].numpy(),
                        sulcus_atlas.numpy())
    NVC.to_csv(tmp_dir+"/NVC.csv")

    affine = nib.load(b0_path).affine
    NVC = cbf_reg["warpedmovout"].numpy()/(alff_reg.numpy() + 1)
    NVC[alff_reg.numpy() == 0] = 0
    nib.Nifti1Image(NVC, affine=affine).to_filename(tmp_dir+"/b0_NVC.nii.gz")


def neurovasc_fMRI(fmri_path, fmri_json_path, cbf_path, save_dir):
    tmp_dir = save_dir+"/tmp"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    # obtain f0 image
    fmri_preproc(fmri_path, fmri_json_path, save_dir)
    f0 = ants.image_read(tmp_dir+"/f0_brain.nii.gz")
    # register sulcus to DTI
    # register f0 to DTI
    warp_to_modality(f0, save_dir, tmp_dir, tmp_dir+"/reg_fmri")

    # summarize ALFF statistics
    map_dict = {"ALFF": tmp_dir+"/alff.nii.gz"}
    region_val_atlas(map_dict, tmp_dir+"/reg_fmri", tmp_dir, thres=5)

    # neurovascular coupling
    get_neurovasc(cbf_path, tmp_dir+"/f0_brain.nii.gz",
                  save_dir+"/DTI/pT2.nii.gz",
                  tmp_dir+"/sulcus_atlas_filled.nii.gz",
                  tmp_dir)
    cleanup(tmp_dir, save_dir)


def cleanup(tmp_dir, wdir):
    # save the NVC and ALFF results
    shutil.copyfile(tmp_dir+"/NVC.csv",
                    wdir+"/CBSS/blood_stats/stats_sulcus_NVC.csv")
    shutil.copyfile(tmp_dir+"/stats_gm_ALFF.csv",
                    wdir+"/CBSS/blood_stats/stats_gm_ALFF.csv")
    shutil.copyfile(tmp_dir+"/stats_sulcus_ALFF.csv",
                    wdir+"/CBSS/blood_stats/stats_sulcus_ALFF.csv")
    shutil.copyfile(tmp_dir+"/b0_NVC.nii.gz",
                    wdir+"/fMRI/figs/b0_NVC.nii.gz")
    shutil.rmtree(tmp_dir)


def ivim_config(config: dict, root: str, ID: str):
    cf = {}
    all_files = ["fMRI", "fMRI_json", "ASL"]
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
            if proceed == 0:
                print(i)
                if not os.path.exists(config["save"]):
                    os.mkdir(config["save"])

                save_goal = config["save"]+"/CBSS/blood_stats/stats_gm_ALFF.csv"
                if not os.path.exists(save_goal):
                    if os.path.exists(config["save"]+"/CBSS/blood_stats"):
                        neurovasc_fMRI(
                            config["fMRI"], config["fMRI_json"],
                            config["save"]+"/CBSS/blood/CBF.nii.gz",
                            config["save"])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_name', type=str, default="CBSS")
    parser.add_argument('--root', type=str)
    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()
    pipeline_all(args.root, args.data_name, num=args.num,
                 arrayID=args.arrayID)
