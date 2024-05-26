'''
Investigate regions in close promixity to the third ventricle

Lateral ventricle:
    1. dilate lateral and third ventricle masks 3 times
    2. find the overlap between the two masks, which is assumed to be the
    foramen of Monro
    3. calculate distance of any voxels away from the overlap region
    4. within lateral ventricle, define 3 zones, which are 5mm, 5-10mm and
    10-15mm away from the foramen of Monro
    5. calculate median and mean IVIm value

Third ventricle:
    1. start from foramen of Mnro
    2. calculate distance of any voxels away from the overlap region

Fourth ventricle:
    1. calculate distance of any voxels away from the fourth ventricle
    2. within the fourth ventricle, define zones as a function of distance away
    from the fourth ventricle
    3. the first zone is defined as the closest region to the fourth ventricle
    that has at least 10 voxels

Define the zones in DTI, which has higher resolution then IVIM
'''
import os
import numpy as np
import pandas as pd
import subprocess
import nibabel as nib
import skimage
import ants


def remove_extreme(val):
    val = val[val > 0]
    val = val[np.isfinite(val)]
    return val


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


def custom_mask(arr, labels):
    out = np.zeros(arr.shape)
    for i in labels:
        out += arr.astype(int) == i
    return (out > 0).astype(float)


def dilation(img, times=3):
    '''times: number of times to dilate the mask'''
    img_mask = img.copy()
    for i in range(times):
        img_mask = skimage.morphology.binary_dilation(img_mask)
    return img_mask


def find_foramen_Monro(vent_mask, times=1):
    '''
    Args:
        `vent_mask`: ventricular parcellation mask (must be in sufficient
        resolution, i.e., in DTI)
    '''
    lateral_ventricle = custom_mask(vent_mask, [4, 43])
    third_ventricle = (vent_mask == 14).astype(float)
    third_ori = third_ventricle.copy()

    for i in range(times):
        lateral_ventricle = dilation(lateral_ventricle)
        third_ventricle = dilation(third_ventricle)

    # find overlap, which must be located in the third ventricle
    overlap = third_ventricle * lateral_ventricle * third_ori
    return overlap.astype(float)


def regions_from_foramen(vent_path, tmp_dir, niimath_path,
                         region_size=5, region_num=4):
    vent_nib = nib.load(vent_path)
    vent_seg = vent_nib.get_fdata()
    # pixdim_z = vent_nib.header["pixdim"][3]
    overlap = find_foramen_Monro(vent_seg, times=1)
    overlap = nib.Nifti1Image(overlap, affine=vent_nib.affine)
    overlap.to_filename(tmp_dir+"/Monro.nii.gz")

    bash_in_python(
        "{}/niimath {}/Monro.nii.gz -binv -edt {}/dist.nii.gz".format(
            niimath_path, tmp_dir, tmp_dir))

    dist = nib.load(tmp_dir+"/dist.nii.gz").get_fdata()
    lateral_ventricle = custom_mask(vent_seg, [4, 43])
    label_map = np.zeros(dist.shape)

    smallest_dist = lateral_ventricle*dist
    smallest_dist = np.min(smallest_dist[smallest_dist > 0])
    largest_dist = region_num*region_size + 1
    regions = (np.arange(0, largest_dist, region_size)).astype(float)
    regions += smallest_dist

    for index, i in enumerate(regions):
        prev_dist = regions[index-1] if index > 0 else 0
        sel_reg = (dist < i)*(dist >= prev_dist)*lateral_ventricle
        sel_reg = sel_reg.astype(float)
        label_map[sel_reg == 1] = index

    label_map = nib.Nifti1Image(label_map, affine=vent_nib.affine)
    label_map.to_filename(tmp_dir+"/regions_from_Monro.nii.gz")


def third_from_foramen(vent_path, tmp_dir, region_size=5, region_num=4):
    vent_nib = nib.load(vent_path)
    vent_seg = vent_nib.get_fdata()

    dist = nib.load(tmp_dir+"/dist.nii.gz").get_fdata()
    third_ventricle = custom_mask(vent_seg, [14])
    label_map = np.zeros(dist.shape)

    smallest_dist = third_ventricle*dist
    smallest_dist = np.min(smallest_dist[smallest_dist > 0])
    largest_dist = region_num*region_size + 1
    regions = (np.arange(0, largest_dist, region_size)).astype(float)
    regions += smallest_dist

    for index, i in enumerate(regions):
        prev_dist = regions[index-1] if index > 0 else 0
        sel_reg = (dist < i)*(dist >= prev_dist)*third_ventricle
        sel_reg = sel_reg.astype(float)
        label_map[sel_reg == 1] = index

    label_map = nib.Nifti1Image(label_map, affine=vent_nib.affine)
    label_map.to_filename(tmp_dir+"/third_from_Monro.nii.gz")


def dist_from_third_vent(vent_path, tmp_dir, niimath_path,
                         region_size=5, region_num=4, thres=10):
    '''thres: number of voxels in the fourth ventricle to be considered as
    being in the fourth ventricle'''
    vent_nib = nib.load(vent_path)
    vent_seg = vent_nib.get_fdata()
    third_ventricle = (vent_seg == 14).astype(float)
    third_ventricle = nib.Nifti1Image(third_ventricle, affine=vent_nib.affine)
    third_ventricle.to_filename(tmp_dir+"/third_ventricle.nii.gz")

    bash_in_python(
        "{}/niimath {}/third_ventricle.nii.gz".format(niimath_path, tmp_dir) +
        " -binv -edt {}/dist_third.nii.gz".format(tmp_dir))

    dist = nib.load(tmp_dir+"/dist_third.nii.gz").get_fdata()
    fourth_ventricle = (vent_seg == 15).astype(float)
    label_map = np.zeros(dist.shape)

    smallest_dist = fourth_ventricle*dist
    smallest_dist = np.min(smallest_dist[smallest_dist > 0])
    largest_dist = region_num*region_size + 1
    regions = (np.arange(0, largest_dist, region_size)).astype(float)
    regions += smallest_dist

    for index, i in enumerate(regions):
        prev_dist = regions[index-1] if index > 0 else 0
        sel_reg = (dist < i)*(dist >= prev_dist)*fourth_ventricle
        sel_reg = sel_reg.astype(float)

        if sel_reg.sum() > thres:
            label_map[sel_reg == 1] = index

    label_map = nib.Nifti1Image(label_map, affine=vent_nib.affine)
    label_map.to_filename(tmp_dir+"/regions_from_third.nii.gz")


def zone_stats(diff_mat, zone_mask, key="Monro"):
    all_regions = np.unique(zone_mask)
    all_regions = all_regions[all_regions != 0]

    out_dict = {}
    for index, i in enumerate(all_regions):
        vals = remove_extreme(diff_mat[zone_mask == i])
        if len(vals) > 0:
            out_dict[key+"_"+"zone"+str(int(i))] = {
                "mean_val": np.mean(vals),
                "median_val": np.median(vals),
                "max_val": np.max(vals),
                "min_val": np.min(vals),
                "sd": np.std(vals),
                "vol": len(vals)
            }
        else:
            out_dict[key+"_"+"zone"+str(int(i))] = {
                "mean_val": np.nan,
                "median_val": np.nan,
                "max_val": np.nan,
                "min_val": np.nan,
                "sd": np.nan,
                "vol": 0
            }
    return pd.DataFrame.from_dict(out_dict, orient="index")


def warp_to_IVIM(vent_diff_path, b0_ivim_path, b0_dti_path, vent_path,
                 niimath_path, tmp_dir, save_path=None, registration=True):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    regions_from_foramen(vent_path, tmp_dir, niimath_path)
    third_from_foramen(vent_path, tmp_dir)
    dist_from_third_vent(vent_path, tmp_dir, niimath_path)

    lateral = ants.image_read(tmp_dir+"/regions_from_Monro.nii.gz")
    third = ants.image_read(tmp_dir+"/third_from_Monro.nii.gz")
    fourth = ants.image_read(tmp_dir+"/regions_from_third.nii.gz")

    if registration:
        b0_ivim = ants.image_read(b0_ivim_path)
        b0_dti = ants.image_read(b0_dti_path)
        tx = ants.registration(b0_ivim, b0_dti, type_of_transform="Rigid",
                               outprefix=tmp_dir+"/ants")
        lateral = ants.apply_transforms(
                b0_ivim, lateral, tx["fwdtransforms"],
                interpolator="genericLabel")
        third = ants.apply_transforms(b0_ivim, third, tx["fwdtransforms"],
                                      interpolator="genericLabel")
        fourth = ants.apply_transforms(b0_ivim, fourth, tx["fwdtransforms"],
                                       interpolator="genericLabel")

    vent_diff = nib.load(vent_diff_path).get_fdata()
    lateral_dict = zone_stats(vent_diff, lateral.numpy())
    third_dict = zone_stats(vent_diff, third.numpy(), key="third")
    fourth_dict = zone_stats(vent_diff, fourth.numpy(), key="fourth")
    out_dict = pd.concat([lateral_dict, third_dict, fourth_dict], axis=0)

    if save_path is None:
        return out_dict
    else:
        # out_df = pd.DataFrame.from_dict(out_dict, orient="index")
        out_dict.to_csv(save_path)

# if __name__ == "__main__":
#     import argparse
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--root', type=str)
#     parser.add_argument('--niimath_path', type=str)
#     parser.add_argument('--save_path', type=str)
#     parser.add_argument('--registration', type=int, default=1)
#     args = parser.parse_args()
#
#     all_df = []
#     for i in sorted(os.listdir(args.root)):
#         root_dir = args.root+"/"+i
#         if os.path.exists(root_dir+"/CBSS/vent_diff.nii.gz"):
#             one_df = warp_to_IVIM(
#                 root_dir+"/CBSS/vent_diff.nii.gz",
#                 root_dir+"/CBSS/b0_ivim.nii.gz",
#                 root_dir+"/DTI/b0.nii.gz",
#                 root_dir+"/CBSS/vent_parcel.nii.gz",
#                 args.niimath_path,
#                 root_dir+"/CBSS/tmp/",
#                 registration=args.registration == 1)
#             one_df.columns = [i]
#             all_df.append(one_df)
#
#     all_df = pd.concat(all_df, axis=1)
#     all_df.T.to_csv(args.save_path)
