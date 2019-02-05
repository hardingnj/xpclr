import pandas as pd
import allel
import numpy as np
import logging
logger = logging.getLogger(__name__)


# FUNCTIONS
def load_hdf5_data(hdf5_fn, chrom, s1, s2, gdistkey=None):

    import hdf5

    samples1 = get_sample_ids(s1)
    samples2 = get_sample_ids(s2)

    samples_x = h5py.File(hdf5_fn)[chrom]["samples"][:]
    sample_name = [sid.decode() for sid in samples_x.tolist()]

    idx1 = np.array([sample_name.index(sid) for sid in samples1])
    idx2 = np.array([sample_name.index(sid) for sid in samples2])

    h5fh = h5py.File(hdf5_fn, mode="r")[chrom]
    g = allel.GenotypeChunkedArray.from_hdf5(h5fh["calldata"]["genotype"])

    pos = allel.SortedIndex(h5fh["variants"]["POS"][:])
    if gdistkey is not None:
        gdist = h5fh["variants"][gdistkey][:]
    else:
        gdist = None

    return g.take(idx1, axis=1), g.take(idx2, axis=1), pos, gdist


def load_zarr_data(zarr_fn, chrom, s1, s2, gdistkey=None):

    import zarr

    samples1 = get_sample_ids(s1)
    samples2 = get_sample_ids(s2)

    zfh = zarr.open_group(zarr_fn, mode="r")[chrom]

    samples_x = zfh["samples"][:]
    sample_name = [sid.decode() for sid in samples_x.tolist()]

    idx1 = np.array([sample_name.index(sid) for sid in samples1])
    idx2 = np.array([sample_name.index(sid) for sid in samples2])

    g = allel.GenotypeChunkedArray(zfh["calldata"]["genotype"])

    pos = allel.SortedIndex(zfh["variants"]["POS"][:])
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
    geno1 = allel.GenotypeChunkedArray(d1.reshape((d1.shape[0], -1, 2)))

    d2 = np.loadtxt(pop_b_fn, dtype="int8")
    geno2 = allel.GenotypeChunkedArray(d2.reshape((d2.shape[0], -1, 2)))

    return geno1, geno2, allel.SortedIndex(vartbl.POS[:]), vartbl.GDist[:]


# function that either splits a string or reads a file
def get_sample_ids(sample_input):

    if "," in sample_input:
        # assume split and return
        logger.debug("Assuming sample IDs given as comma-separated strings.")
        samples = sample_input.split(",")

    else:
        logger.debug("Assuming sample IDs provided in a file.")
        with open(sample_input, "r") as reader:
            samples = [x.strip() for x in reader.readlines()]

    return samples


def load_vcf_wrapper(path, seqid, samples, samples_path):

    callset = allel.read_vcf(
        path,
        region=seqid,
        fields=['variants/POS', 'calldata/GT', 'samples'],
        tabix="tabix",
        samples=samples)

    assert "samples" in callset.keys(), "None of the samples provided in {0!r} are found in {1!r}".format(
        samples_path, path)

    p = allel.SortedIndex(callset["variants/POS"])
    g = allel.GenotypeArray(callset['calldata/GT'])

    return p, g


def load_vcf_format_data(vcf_fn, chrom, s1, s2, gdistkey=None):

    #    geno1, geno2, pos = q, q, q
    samples1 = get_sample_ids(s1)
    samples2 = get_sample_ids(s2)
    pos1, geno1 = load_vcf_wrapper(vcf_fn, chrom, samples1, s1)
    pos2, geno2 = load_vcf_wrapper(vcf_fn, chrom, samples2, s2)

    assert np.array_equal(pos1, pos2), "POS fields not the same"
    assert geno1.shape[0] == pos1.shape[0], "For samples 1, genotypes do not match positions"
    assert geno2.shape[0] == pos2.shape[0], "For samples 2, genotypes do not match positions"
    assert geno1.shape[1] == len(samples1)
    assert geno2.shape[1] == len(samples2)

    return geno1, geno2, pos1, None


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

