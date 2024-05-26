# obtain the values in IVIM
import re
import os
import numpy as np
import ants
import pandas as pd
from src import fort
from get_fw import extract_b0


# register CSF mask to IVIM
def reg_csf_img(ivim_b0, dwi_b0, CSF_mask, wdir):
    tx = ants.registration(ivim_b0, dwi_b0, type_of_transforms="Rigid",
                           outprefix=wdir+"/tmp/ants")
    reg_csf = ants.apply_transforms(ivim_b0, CSF_mask, tx["fwdtransforms"],
                                    interpolator="genericLabel")
    return reg_csf


def warp_all_maps(wdir, ivim_b0, dwi_b0):
    for i in ["CSF_parcel.nii.gz", "vent_parcel.nii.gz",
              "sulcus_cistern.nii.gz"]:
        CSF_mask = ants.image_read(wdir+"/"+i)
        reg_csf = reg_csf_img(ivim_b0, dwi_b0, CSF_mask, wdir)
        out_name = re.sub(".nii.gz", "_reg.nii.gz", i)
        ants.image_write(reg_csf, wdir+"/"+out_name)


def reg_csf_wrap(ivim_path, ivim_b0_path, dwi_path, dwi_b0_path, wdir,
                 ivim_b_thres=200):
    if not os.path.exists(wdir+"/tmp"):
        os.mkdir(wdir+"/tmp")

    extract_b0(ivim_path, ivim_b0_path, wdir+"/tmp/b0_ivim.nii.gz",
               threshold=2)
    extract_b0(dwi_path, dwi_b0_path, wdir+"/tmp/b0_dwi.nii.gz",
               threshold=10)

    ivim = ants.image_read(ivim_path)
    ivim_b0 = ants.image_read(wdir+"/tmp/b0_ivim.nii.gz")
    dwi_b0 = ants.image_read(wdir+"/tmp/b0_dwi.nii.gz")
    warp_all_maps(wdir, ivim_b0, dwi_b0)

    ivim = np.moveaxis(ivim.numpy(), 3, 0)
    bvals = np.genfromtxt(ivim_b0_path).reshape([-1])
    ivim = ivim[bvals <= ivim_b_thres]
    bvals = bvals[bvals <= ivim_b_thres]

    out_dict = {"CSF": "CSF_parcel_reg.nii.gz",
                "vent": "vent_parcel_reg.nii.gz",
                "sulcus": "sulcus_cistern_reg.nii.gz"}

    for key, val in out_dict.items():
        reg_csf = ants.image_read(wdir+"/"+val)
        csf_mask = reg_csf.numpy() != 0

        fast_diffu = fort.fit_ivim(ivim, bvals, ivim_b0.numpy(), csf_mask)
        fast_diffu = ants.from_numpy(fast_diffu, origin=ivim_b0.origin,
                                     direction=ivim_b0.direction,
                                     spacing=ivim_b0.spacing)
        ants.image_write(fast_diffu, wdir+"/{}_diff.nii.gz".format(key))


def extract_ivim_val(wdir, mode="CSF"):
    fast_diffu = ants.image_read(wdir+"/{}_diff.nii.gz".format(mode)).numpy()
    if mode != "sulcus":
        mask_mode = mode+"_parcel_reg.nii.gz"
    else:
        mask_mode = "sulcus_cistern_reg.nii.gz"
    CSF_mask = ants.image_read(wdir+"/{}".format(mask_mode)
                               ).numpy()
    all_regions = np.unique(CSF_mask)
    all_regions = all_regions[all_regions != 0]

    out_dict = {}
    for i in all_regions:
        if np.sum(CSF_mask == i) > 20:
            diff_val = fast_diffu[CSF_mask == i]
            # remove the values below 0
            diff_val = diff_val[diff_val > 0]
            diff_val = diff_val[np.isfinite(diff_val)]
            out_dict[i] = {
                "mean_val": diff_val.mean(),
                "median_val": np.median(diff_val),
                "max_val": diff_val.max(),
                "min_val": diff_val.min(),
                "sd": diff_val.std(),
                "vol": len(diff_val)
            }

    out_dict = pd.DataFrame.from_dict(out_dict, orient="index")
    out_dict.to_csv(wdir+"/sum_stat_{}.csv".format(mode))


def extract_global_val(wdir):
    out_dict = {}
    all_vals = []
    out_dict = {"CSF": "CSF_parcel_reg.nii.gz",
                "vent": "vent_parcel_reg.nii.gz",
                "sulcus": "sulcus_cistern_reg.nii.gz"}
    for mode, filename in out_dict.items():
        mode2 = "CSF" if mode == "sulcus" else mode
        fast_diffu = ants.image_read(
            wdir+"/{}_diff.nii.gz".format(mode2)).numpy()
        CSF_mask = ants.image_read(
            wdir+"/{}".format(filename)).numpy()
        diff_val = fast_diffu[CSF_mask > 0]
        # remove the values below 0
        diff_val = diff_val[diff_val > 0]
        diff_val = diff_val[np.isfinite(diff_val)]
        all_vals.append(diff_val)
        out_dict[mode] = {
            "mean_val": diff_val.mean(),
            "median_val": np.median(diff_val),
            "max_val": diff_val.max(),
            "min_val": diff_val.min(),
            "sd": diff_val.std(),
            "vol": len(diff_val)
        }

    all_vals = np.concatenate(all_vals, axis=0)
    out_dict["all"] = {
        "mean_val": all_vals.mean(),
        "median_val": np.median(all_vals),
        "max_val": all_vals.max(),
        "min_val": all_vals.min(),
        "sd": all_vals.std(),
        "vol": len(all_vals)
    }
    out_dict = pd.DataFrame.from_dict(out_dict, orient="index")
    out_dict.to_csv(wdir+"/sum_stat_global.csv")


def get_ivim(ivim_bval_path, dwi_bval_path, save_dir):
    reg_csf_wrap(save_dir+"/DWI/corr_ivim.nii.gz",
                 ivim_bval_path,
                 save_dir+"/DWI/corr.nii.gz",
                 dwi_bval_path,
                 save_dir)

    for i in ["CSF", "vent", "sulcus"]:
        extract_ivim_val(save_dir, mode=i)
    extract_global_val(save_dir)
