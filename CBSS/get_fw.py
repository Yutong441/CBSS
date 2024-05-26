import os
import re
import subprocess
import numpy as np
import nibabel as nib


def read_bval(filename):
    '''
    Read a bval file as numpy array
    '''
    with open(filename) as f:
        bval = f.readlines()

    try:
        bval = "".join(bval).split(" ")
        bval = [int(re.sub("\\\n", "", i)) for i in bval]
    except Exception:
        bval = "".join(bval).split("\t")
        bval = [int(re.sub("\\\n", "", i)) for i in bval]
    return np.array(bval)


def bash_in_python(bashcommand):
    ID = ''.join([str(i) for i in np.random.choice(9, 10)])
    sh_script = 'tmp_test'+ID+'.sh'
    with open(sh_script, 'w') as f:
        f.write('#!/bin/bash \n')
        f.write(bashcommand)
    subprocess.run(['chmod', '+x', sh_script], stderr=subprocess.PIPE)
    subprocess.call('./'+sh_script)
    os.remove(sh_script)


def b0_indices(bval_path, threshold=10):
    with open(bval_path, "r") as f:
        bvals = f.readline()

    bvals = re.sub("\n", "", bvals)
    bvals = re.sub("\t", " ", bvals)
    bvals = [i for i in bvals.split(" ") if i != ""]
    ind_list = []
    for index, i in enumerate(bvals):
        if int(i) <= threshold:  # sometimes b0 volumes are acquired under b=5
            ind_list.append(index)

    return ind_list


def extract_b0(DTI_path, bval_path, save_path, threshold=10):
    # img_nib = nib.load(DTI_path)
    # affine = img_nib.affine
    # img = img_nib.get_fdata()
    # del img_nib

    b0_ind = b0_indices(bval_path, threshold=threshold)
    b0_indices_list = ",".join([str(i) for i in b0_ind])
    print("extracting volume {}".format(b0_indices_list))
    # b0_all = np.stack([img[..., i] for i in b0_ind], axis=-1)

    b0_all_path = os.path.dirname(save_path)+"/b0_all.nii.gz"
    cmd = "fslselectvols -i {} -o {} --vols={};".format(
        DTI_path, b0_all_path, b0_indices_list
    )
    cmd += "fslmaths {} -Tmean {}".format(
        b0_all_path, save_path
    )

    bash_in_python(cmd)
    # b0_avg = b0_all.mean(-1)
    # b0_out = nib.Nifti1Image(b0_avg, affine=affine)
    # b0_out.to_filename(save_path)


def extract_single_shell(data_path, bval_path, bvec_path, save_dir,
                         single_shell_thres=1300):
    bvecs = np.genfromtxt(bvec_path)
    bval = np.genfromtxt(bval_path)
    single_shell = bval < single_shell_thres

    # dwi_nib = nib.load(data_path)
    # affine = dwi_nib.affine
    # dwi = dwi_nib.get_fdata()
    # del dwi_nib
    # dwi = dwi[..., single_shell]
    # dwi_sing = nib.Nifti1Image(dwi, affine=affine)
    # dwi_sing.to_filename(save_dir+"/corr.nii.gz")

    sing_indices_list = np.arange(len(bval))[single_shell]
    sing_indices_list = ",".join([str(i) for i in sing_indices_list])
    cmd = "fslselectvols -i {} -o {} --vols={}".format(
        data_path, save_dir+"/corr.nii.gz", sing_indices_list
    )
    bash_in_python(cmd)

    np.savetxt(save_dir+"/corr_bvecs",
               bvecs[:, bval < single_shell_thres], fmt="%.6f")
    np.savetxt(save_dir+"/corr_bvals",
               bval[bval < single_shell_thres].reshape([1, -1]),
               fmt="%.0f")


def find_fw(data, bval, bvec, save_dir, mask=None, skullstrip=False,
            method="fortran"):
    if not os.path.exists(save_dir+"/FW.nii.gz"):
        extract_b0(data, bval, save_dir+"/b0.nii.gz")
        if mask is None:
            if skullstrip:
                bash_in_python(
                    "bet {}/b0.nii.gz {}/b0_ss.nii.gz -f 0.3 -R -m".format(
                        save_dir, save_dir))
            else:
                img = nib.load(save_dir+"/b0.nii.gz")
                blank = nib.Nifti1Image(np.ones(img.shape), affine=img.affine,
                                        dtype=np.int16)
                blank.to_filename(save_dir+"/b0_ss_mask.nii.gz")

        extract_single_shell(data, bval, bvec, save_dir)
        if method == "fortran":
            from fw_fort import find_fw_sing_shell
            find_fw_sing_shell(
                save_dir+"/corr.nii.gz",
                save_dir+"/corr_bvals",
                save_dir+"/corr_bvecs",
                save_dir,
                save_dir+"/b0_ss_mask.nii.gz"
            )
        # else:
        #     from fw_mrn import find_fw_sing_shell
        #     find_fw_sing_shell(save_dir+"/corr.nii.gz",
        #                        save_dir+"/b0_ss_mask.nii.gz",
        #                        save_dir+"/corr_bvals",
        #                        save_dir+"/corr_bvecs",
        #                        save_dir)

        os.remove(save_dir+"/corr.nii.gz")
        os.remove(save_dir+"/corr_bvals")
        os.remove(save_dir+"/corr_bvecs")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=str)
    parser.add_argument('--bval', type=str)
    parser.add_argument('--bvec', type=str)
    parser.add_argument('--save_dir', type=str)
    args = parser.parse_args()
    find_fw(args.data, args.bval, args.bvec, args.save_dir)
