import pandas as pd
import h5py
import allel
import numpy as np


#FUNCTIONS
def load_hdf5_data(hdf5_fn, chrom, s1, s2, gdistkey=None):

    samples = h5py.File(hdf5_fn)[chrom]["samples"][:]
    sample_name = [sid.decode() for sid in samples.tolist()]

    idx1 = np.array([sample_name.index(sid) for sid in s1])
    idx2 = np.array([sample_name.index(sid) for sid in s2])

    h5fh = h5py.File(hdf5_fn, mode="r")[chrom]
    g = allel.GenotypeCArray.from_hdf5(h5fh["calldata"]["genotype"])

    pos = allel.SortedIndex(h5fh["variants"]["POS"][:])
    if gdistkey is not None:
        gdist = h5fh["variants"][gdistkey][:]
    else:
        gdist = None

    return g.take(idx1, axis=1), g.take(idx2, axis=1), pos, gdist


def load_text_format_data(mapfn, pop_a_fn, pop_b_fn):

    tbl = pd.read_csv(mapfn, sep=" ",
                      names=["ID", "CHROM", "GDist", "POS", "REF", "ALT"])

    vartbl = allel.VariantChunkedTable(tbl.to_records(), index="POS")

    d1 = np.loadtxt(pop_a_fn, dtype="int8")
    geno1 = allel.GenotypeCArray(d1.reshape((d1.shape[0], -1, 2)))

    d2 = np.loadtxt(pop_b_fn, dtype="int8")
    geno2 = allel.GenotypeCArray(d2.reshape((d2.shape[0], -1, 2)))

    return geno1, geno2, allel.SortedIndex(vartbl.POS[:]), vartbl.GDist[:]


def tabulate_results(chrom, model_li, null_li, selectionc,
                     counts, count_avail, windows, edges):

    lidf = pd.DataFrame(np.vstack((model_li, null_li, selectionc, counts, count_avail)).T,
                        columns=["modelL", "nullL", "sel_coef", "nSNPs", "nSNPs_avail"])

    # these are the nominal windows
    winf = pd.DataFrame(windows, columns=["start", "stop"])

    # these are the "real" windows. Gives a guide to how close we are.
    realf = pd.DataFrame(edges, columns=["pos_start", "pos_stop"])

    out = pd.concat([winf, realf, lidf], axis=1)

    out["xpclr"] = 2 * (out.modelL - out.nullL)
    out["xpclr_norm"] = (out.xpclr - np.nanmean(out.xpclr))/np.nanstd(out.xpclr)

    out.insert(0, "chrom", np.repeat(chrom, len(out)))

    string_id = ["{0}_{1:08d}_{2:08d}".format(r.chrom, r.start, r.stop)
                 for i, r in out.iterrows()]
    out.insert(0, "id", string_id)

    return out

