# create pseudo-T2 contrast
import os
import re
import subprocess
import shutil
import numpy as np
import ants
import nibabel as nib
from get_fw import find_fw


def WM_prob(FA):
    mask = ants.get_mask(FA)
    out = ants.atropos(a=FA, i='kmeans[2]', x=mask, m='[0.2,1x1x1]')
    cor0 = np.corrcoef(out["probabilityimages"][0].flatten(), FA.flatten())
    cor1 = np.corrcoef(out["probabilityimages"][1].flatten(), FA.flatten())
    if cor0[0, 1] > cor1[0, 1]:
        return out["probabilityimages"][0]
    else:
        return out["probabilityimages"][1]
    return out


def pseudo_T2(fw_path, fa_path, mask_path, save_path):
    mask = ants.image_read(mask_path)
    fw_nib = nib.load(fw_path)
    affine = fw_nib.affine
    del fw_nib

    FW = ants.image_read(fw_path).numpy()
    FA = ants.image_read(fa_path)
    FA = ants.mask_image(FA, mask)
    WM = WM_prob(FA)
    WM = WM.numpy()
    GM = 1 - (FW + WM)
    T2 = GM + 2*FW
    T2 = nib.Nifti1Image(T2*mask.numpy(), affine=affine)
    T2.to_filename(save_path)


def pseudo_T1(fw_path, fa_path, mask_path, save_path):
    mask = ants.image_read(mask_path)
    fw_nib = nib.load(fw_path)
    affine = fw_nib.affine
    del fw_nib

    FW = ants.image_read(fw_path).numpy()
    FA = ants.image_read(fa_path)
    FA = ants.mask_image(FA, mask)
    WM = WM_prob(FA)
    WM = WM.numpy()
    GM = 1 - (FW + WM)
    T1 = GM + 2*WM + 0.2*FW
    # consider adding non-zero weight to T1
    T1 = nib.Nifti1Image(T1*mask.numpy(), affine=affine)
    T1.to_filename(save_path)


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


def skullstrip(skull_path, save_path):
    if not os.path.exists(save_path):
        mask_path = re.sub(".nii.gz", "_mask.nii.gz", save_path)
        cmd = "mri_synthstrip -i {} -o {} -m {}".format(
            skull_path, save_path, mask_path)
        bash_in_python(cmd)


def cbss_1_reg(fn_data, fn_bval, fn_bvec, save_dir):
    if not os.path.exists(save_dir+"/pT2.nii.gz"):
        tmp_dir = save_dir+"/tmp"
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)

        find_fw(fn_data, fn_bval, fn_bvec, tmp_dir)
        skullstrip(tmp_dir+"/FW.nii.gz", tmp_dir+"/FW_synth.nii.gz")
        pseudo_T2(tmp_dir+"/FW.nii.gz", tmp_dir+"/FA.nii.gz",
                  tmp_dir+"/FW_synth_mask.nii.gz",
                  tmp_dir+"/pT2.nii.gz")

        save_files = ["FW.nii.gz", "FA.nii.gz", "pT2.nii.gz",
                      "FW_synth_mask.nii.gz"]
        for i in save_files:
            shutil.copyfile(tmp_dir+"/"+i, save_dir+"/"+i)
        shutil.rmtree(tmp_dir)


def make_template(img_path_list, save_dir):
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    template = os.path.dirname(os.path.realpath(__file__)) + \
        "/../data/MNI_pT2.nii.gz"

    tmp_dir = save_dir+"/tmp"
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    # for key, val in img_path_list.items():
    #     shutil.copyfile(val, tmp_dir+"/"+key+".nii.gz")

    # cmd = "cd {};".format(tmp_dir) + \
    #     "buildtemplateparallel.sh -d 3 -j 1 -o {}/pT2_".format(save_dir) + \
    #     " -n 0 -s MI -i 4 -m 30x50x20 -t GR -z {}".format(template) + \
    #     " *.nii.gz; cd -"
    # bash_in_python(cmd)

    template = ants.image_read(template)
    img_list = [ants.image_read(i) for i in img_path_list]

    built = ants.build_template(
        initial_template=template,
        image_list=img_list,
        iterations=4,
        gradient_step=0.2,
        blending_weight=0.75,
        weights=None,
        useNoRigid=True,
        outprefix=tmp_dir+"/ants"
    )
    ants.image_write(built, save_dir+"/group_template.nii.gz")
    shutil.rmtree(tmp_dir)


def warp_to_template(img_dir, template_path, free_water_thres=0.8):
    pT2 = ants.image_read(img_dir+"/pT2.nii.gz")
    template = ants.image_read(template_path)
    tmp_dir = img_dir+"/tmp"
    free_water = ants.image_read(img_dir+"/FW.nii.gz")
    mask = ants.image_read(img_dir+"/FW_synth_mask.nii.gz")
    free_water = ants.mask_image(free_water, mask)
    csf_mask = ants.threshold_image(free_water, low_thresh=free_water_thres,
                                    binary=True)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    tx = ants.registration(template, pT2, type_of_transform="SyN",
                           outprefix=tmp_dir+"/ants")
    mask_templ = ants.apply_transforms(template, csf_mask, tx["fwdtransforms"])
    ants.image_write(mask_templ, img_dir+"/csf_templ.nii.gz")
    shutil.rmtree(tmp_dir)
