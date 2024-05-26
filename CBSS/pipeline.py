import os
import shutil
import json
import glob
import re

import cbss_1_reg as cb1
import cbss_2_parc as cb2
import cbss_3_ivim as cb3
import cbss_4_vent as cb4


def makedir(wdir):
    if not os.path.exists(wdir):
        os.mkdir(wdir)


def clean_up(wdir):
    # save CBSS file
    makedir(wdir+"/CBSS")
    makedir(wdir+"/CBSS/stats")
    makedir(wdir+"/CBSS/parcel")
    makedir(wdir+"/CBSS/diff")
    makedir(wdir+"/CBSS/DTI")

    os.rename(wdir+"/tmp/b0_ivim.nii.gz", wdir+"/CBSS/b0_ivim.nii.gz")
    # shutil.move(wdir+"/sulcus_reg", wdir+"/CBSS/sulcus_reg")
    for i in ["FW.nii.gz", "FW_synth_mask.nii.gz", "pT2.nii.gz",
              "seg.nii.gz"]:
        shutil.copy2(wdir+"/"+i, wdir+"/CBSS/DTI/")
        os.remove(wdir+"/"+i)

    for i in glob.glob(wdir+"/*parcel*.nii.gz"):
        shutil.copy2(i, wdir+"/CBSS/parcel/")
        os.remove(i)

    for i in glob.glob(wdir+"/sulcus_cistern*.nii.gz"):
        shutil.copy2(i, wdir+"/CBSS/parcel/")
        os.remove(i)

    for i in glob.glob(wdir+"/*diff.nii.gz"):
        shutil.copy2(i, wdir+"/CBSS/diff/")
        os.remove(i)

    for i in glob.glob(wdir+"/*.csv"):
        shutil.copy2(i, wdir+"/CBSS/stats")
        os.remove(i)

    # save DTI and IVIM files
    if not os.path.exists(wdir+"/DWI"):
        os.mkdir(wdir+"/DWI")
    for i in ["corr.nii.gz", "corr_ivim.nii.gz"]:
        if os.path.exists(wdir+"/"+i):
            shutil.copy2(wdir+"/"+i, wdir+"/DWI")

    shutil.rmtree(wdir+"/tmp")
    for i in glob.glob(wdir+"/*.nii.gz"):
        if os.path.exists(i):
            os.remove(i)


def preproc(inp_path, out_path):
    if not os.path.exists(out_path):
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        out_name = os.path.basename(out_path)
        inp_name = os.path.basename(inp_path)
        ori_dir = os.getcwd()

        shutil.copy2(inp_path, out_dir)
        os.chdir(out_dir)
        cb1.bash_in_python("fslreorient2std {} reori.nii.gz".format(inp_name))
        cb1.bash_in_python("eddy_correct reori.nii.gz {} 0".format(out_name))

        os.remove("reori.nii.gz")
        os.remove(inp_name)
        os.remove(re.sub(".nii.gz", ".ecclog", out_name))
        os.chdir(ori_dir)


def cbss(dwi_data_path, dwi_bval_path, dwi_bvec_path,
         ivim_data_path, ivim_bval_path, niimath_path,
         save_dir):
    '''
    Output: in the `save_dir`:
    DWI/corr.nii.gz
    DWI/corr_bvals
    DWI/corr_bvecs
    DWI/corr_ivim.nii.gz

    # image files
    CBSS/CSF_diff.nii.gz: pseudodiffusion ADC in the peripheral CSF region
    CBSS/vent_diff.nii.gz: pseudodiffusion ADC in the ventricle region
    CBSS/CSF_parcel_reg.nii.gz: peripheral CSF parcellation
    CBSS/vent_parcel_reg.nii.gz: ventricle parcellation
    CBSS/b0_ivim.nii.gz
    CBSS/FW.nii.gz: free water map (without skullstripping)
    CBSS/pT2.nii.gz: pseudo-T2 contrast image

    # summary statistics
    CBSS/sum_stat_global.csv
    CBSS/sum_stat_csf.csv
    CBSS/sum_stat_vent.csv
    '''
    # do motion correction
    preproc(dwi_data_path, save_dir+"/DWI/corr.nii.gz")
    preproc(ivim_data_path, save_dir+"/DWI/corr_ivim.nii.gz")

    # CBSS
    cb1.cbss_1_reg(save_dir+"/DWI/corr.nii.gz", dwi_bval_path, dwi_bvec_path,
                   save_dir)
    cb2.parcellate(save_dir)
    cb3.get_ivim(ivim_bval_path, dwi_bval_path, save_dir)
    cb4.warp_to_IVIM(
        save_dir+"/vent_diff.nii.gz",
        save_dir+"/tmp/b0_ivim.nii.gz",
        save_dir+"/tmp/b0_dwi.nii.gz",
        save_dir+"/vent_parcel.nii.gz",
        niimath_path, save_dir+"/tmp/",
        save_dir+"/sum_vent_region.csv",
        registration=True)

    # clean up folders
    clean_up(save_dir)
    shutil.copyfile(dwi_bval_path, save_dir+"/DWI/corr_bvals")
    shutil.copyfile(dwi_bvec_path, save_dir+"/DWI/corr_bvecs")


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


def ivim_config(config: dict, root: str, ID: str):
    cf = {}
    all_files = ["DWI", "bvec", "bval", "IVIM", "IVIM_bval"]
    proceed = len(all_files)
    for i in all_files:
        cf[i] = _add_ID([config[i]], root, ID)[0]
        if os.path.exists(cf[i]):
            proceed -= 1

    cf["save"] = root+"/"+config["save"]+"/"+ID
    return cf, proceed


def pipeline(root, data_name, niimath_path, num=0, arrayID=0):
    cf_name = os.path.dirname(os.path.realpath(__file__))+"/data_regex.json"
    with open(cf_name, "r") as f:
        cf = json.load(f)
    one_cf = cf[data_name]
    all_IDs = sorted(os.listdir(root+"/"+one_cf["ID"]))

    for index, i in enumerate(all_IDs):
        if index % num == arrayID:
            config, proceed = ivim_config(one_cf, root, i)
            if not os.path.exists(config["save"]):
                os.mkdir(config["save"])

            if proceed == 0:
                if not os.path.exists(config["save"]+"/CBSS/stats/sum_stat_CSF.csv"):
                    cbss(config["DWI"], config["bval"],
                         config["bvec"], config["IVIM"],
                         config["IVIM_bval"], niimath_path, config["save"])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_name', type=str, default="CBSS")
    parser.add_argument('--root', type=str)
    parser.add_argument("--niimath_path", type=str, default="None")
    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()
    if args.niimath_path == "None":
        niimath_p = os.path.realpath(__file__) + \
            "/../software/niimath/build/bin/"
    else:
        niimath_p = args.niimath_path

    pipeline(args.root, args.data_name, num=args.num,
             arrayID=args.arrayID, niimath_path=niimath_p)
