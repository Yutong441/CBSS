import numpy as np
import re
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import nibabel as nib


def disp_overlay(img, overlay, mask=None, save_path=None, width=6,
                 height=3, magni=6, exclude_slices=30, reverse_color=False):
    N = width*height
    D = img.shape[2]
    seq = np.linspace(exclude_slices, D - exclude_slices, num=N)
    seq = [int(i) for i in seq]

    if not reverse_color:
        cmap = mpl.cm.get_cmap("jet").copy()
    else:
        cmap = mpl.cm.get_cmap("jet_r").copy()
    cmap.set_under(color='black')
    mpl.rcParams.update({'font.size': 22})

    if mask is None:
        mask = np.ones(img.shape)
        # mask = (np.abs(overlay) > 0.001).astype(float)

    fig, axes = plt.subplots(height, width, figsize=(magni*width, magni*height))
    maxval = (mask*overlay).max() + 0.001
    minval = (mask*overlay).min() + 0.001

    overlay += (minval - 0.001)*(1 - mask)
    for index, i in enumerate(seq):
        axes[index // width, index % width].imshow(
            np.rot90(img[..., i]), cmap="gray")
        im = axes[index // width, index % width].imshow(
            np.rot90(overlay[..., i]), cmap=cmap, alpha=0.5,
            vmin=minval, vmax=maxval)
        axes[index // width, index % width].axis("off")
        axes[index // width, index % width].set_title("z={}".format(int(i)))

    cbar = plt.colorbar(im, ax=axes.ravel().tolist())
    cbar.ax.set_title('Pseudodiffusivity (mm$^2$/s)')
    if save_path is not None:
        plt.savefig(save_path, dpi=400, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def generate_overlay(parcel_map, labels, values):
    '''
    Args:
        `parcel_map`: parcellation atlas
        `value_dict`: a dictionary of labels (1, 2, 3 etc) and the value of
        diffusivitiy, correlation coeffficients etc
    '''
    overlay = np.zeros(parcel_map.shape)
    for key, val in zip(labels, values):
        if not np.isfinite(val):
            val = 0
        overlay[parcel_map == key] = val
    return overlay


def rename_labels(ori_names, conversion_df_path=None):
    if conversion_df_path is None:
        conversion_df_path = os.path.dirname(os.path.realpath(__file__)) + \
            "/../data/CBSS_atlas.csv"

    conversion = pd.read_csv(conversion_df_path)
    new_names = ori_names.copy()
    for index, i in enumerate(ori_names):
        new_names[index] = conversion[conversion["atlas"] == i]["label"].values
        new_names[index] = int(re.sub("C", "", list(new_names[index])[0]))
    return new_names


def plot_values(df_path, save_dir, parcel_path=None, template_path=None,
                **kwargs):
    df = pd.read_csv(df_path, index_col=[0])
    relabel = rename_labels(list(df.index))
    if parcel_path is None:
        parcel_path = os.path.dirname(os.path.realpath(__file__)) + \
            "/../data/sulcus_vent.nii.gz"

    if template_path is None:
        template_path = os.environ["FSLDIR"] + \
            "/data/standard/MNI152_T1_1mm_brain.nii.gz"

    MNI = nib.load(template_path).get_fdata()
    parcel_nib = nib.load(parcel_path)
    parcel_map = parcel_nib.get_fdata()
    for feature in df.columns:
        overlay = generate_overlay(parcel_map, relabel, df[feature].values)
        save_path = save_dir+"/Map_"+feature+".jpg"
        disp_overlay(MNI, overlay, mask=parcel_map > 0,
                     save_path=save_path, **kwargs)
        overlay = nib.Nifti1Image(overlay, affine=parcel_nib.affine)
        overlay.to_filename(save_dir+"/Map_"+feature+".nii.gz")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--df_path', type=str,
                        default="results/sum_stats.csv")
    parser.add_argument('--save_dir', type=str,
                        default="reports/figures/")
    args = parser.parse_args()
    plot_values(args.df_path, args.save_dir)
