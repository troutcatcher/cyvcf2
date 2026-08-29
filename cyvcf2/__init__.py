from ._version import __version__
from .cyvcf2 import (VCF, Variant, Writer, r_ as r_unphased, par_relatedness,
                     par_het)
Reader = VCFReader = VCF


def read_gt012(path, transpose=False, gts012=True, strict_gt=False,
               parse_threads=None, shards=None):
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

    shards=N instead splits a tabix-indexed .vcf.gz into N ranges at
    record boundaries taken from the index (hts_idx_split) and reads them
    on N threads, each with its own decompressor and parser - this scales
    with cores better than the parse_threads pipeline. It falls back to
    the parse_threads path when the file cannot be sharded (no index, too
    small, BCF, or an htslib without the API).
    """
    if shards and int(shards) > 1 and not str(path).endswith((".bcf", ".bcf.gz")):
        parts = _read_gt012_sharded(path, int(shards), gts012, strict_gt)
        if parts is not None:
            import numpy as np
            out = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=0)
            return out.T if transpose else out
        if parse_threads is None:
            parse_threads = int(shards)
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


def _read_gt012_sharded(path, n, gts012, strict_gt):
    """Split `path` into byte ranges at record boundaries from its tabix
    index and read them in parallel, one VCF reader per thread with the
    GIL released; returns the per-shard matrices in file order, or None
    when the file cannot be sharded."""
    import threading

    probe = VCF(path)
    try:
        splits = probe._shard_offsets(n)
    finally:
        probe.close()
    if not splits:
        return None

    bounds = [-1] + list(splits) + [-1]
    readers = []
    try:
        for _ in range(len(bounds) - 1):
            readers.append(VCF(path, gts012=gts012, strict_gt=strict_gt,
                               format_fields=["GT"]))
        if any(r.lazy_format_mode != "native" for r in readers):
            return None  # line-stripping fallback cannot seek mid-file
        results = [None] * len(readers)
        errors = []

        def run(i):
            try:
                results[i] = readers[i]._gt012_shard(bounds[i], bounds[i + 1])
            except BaseException as e:
                errors.append(e)

        threads = [threading.Thread(target=run, args=(i,))
                   for i in range(len(readers))]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        if errors:
            raise errors[0]
        return results
    finally:
        for r in readers:
            r.close()
