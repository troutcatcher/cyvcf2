from ._version import __version__
from .cyvcf2 import (VCF, Variant, Writer, r_ as r_unphased, par_relatedness,
                     par_het)
Reader = VCFReader = VCF


def read_gt012(path, transpose=False, gts012=True, strict_gt=False,
               parse_threads=None):
    """Read a whole VCF's genotypes as one 0/1/2 numpy int8 matrix.

    Returns an array of shape (n_variants, n_samples) with HOM_REF=0,
    HET=1, HOM_ALT=2 and UNKNOWN=3 (or, with gts012=False, cyvcf2's legacy
    coding where 2 and 3 are swapped). With transpose=True the returned
    array is the (n_samples, n_variants) transpose - a zero-copy view, so
    it is F-contiguous; wrap it in numpy.ascontiguousarray() if C order is
    required.

    This is the fastest way to get a genotype matrix: the file is opened
    parsing only the GT FORMAT field, the matrix is pre-allocated in full
    using the record count from the tabix/CSI index when one is present,
    and rows are written by a bulk C loop without creating any Variant
    objects. parse_threads adds htslib worker threads that decompress and
    parse lines in the background (e.g. parse_threads=4).
    """
    kwargs = {"gts012": gts012, "strict_gt": strict_gt}
    if not str(path).endswith((".bcf", ".bcf.gz")):
        kwargs["format_fields"] = ["GT"]
        if parse_threads:
            kwargs["parse_threads"] = parse_threads
    elif parse_threads:
        # BCF records are binary; only bgzf decompression can be threaded
        kwargs["threads"] = parse_threads
    vcf = VCF(path, **kwargs)
    try:
        return vcf.gt_types_matrix(transpose=transpose)
    finally:
        vcf.close()
