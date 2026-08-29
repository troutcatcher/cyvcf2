#include <helpers.h>
#include <string.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// map one raw int8 GT value to the int32 value bcf_get_genotypes would produce
static inline int32_t gt8_to_i32(int8_t v) {
    if (v == bcf_int8_missing) return bcf_int32_missing;
    if (v == bcf_int8_vector_end) return bcf_int32_vector_end;
    return v;
}

// genotype classes before HOM_ALT/UNKNOWN are mapped for gts012
enum { GT012_HOM_REF = 0, GT012_HET = 1, GT012_HOM_ALT = 2, GT012_UNKNOWN = 3 };

// classify one genotype of `ploidy` raw int8 GT values into a GT012_* class.
// This mirrors as_gts exactly: the int8 -> int32 sentinel mapping preserves
// every comparison as_gts makes on the values bcf_get_genotypes produces.
// The single source of truth for the 012 semantics: the diploid lookup table
// and both generic-ploidy writers below are all derived from it.
static int gt012_classify8(const int8_t *g, int ploidy, int strict_gt) {
    int k, missing = 0;
    for (k = 0; k < ploidy; k++) {
        if (bcf_gt_is_missing(gt8_to_i32(g[k]))) missing += 1;
    }
    if (missing == ploidy || (missing != 0 && strict_gt)) return GT012_UNKNOWN;
    if (ploidy == 1 || g[1] == bcf_int8_vector_end) {
        int a = bcf_gt_allele(g[0]);
        if (a == 0) return GT012_HOM_REF;
        if (a == 1) return GT012_HOM_ALT;
        return GT012_UNKNOWN;
    }
    {
        int a = bcf_gt_allele(g[0]);
        int b = bcf_gt_allele(g[1]);
        // 0/. counts as hom ref: a single missing allele has no alts
        if ((a == 0 && b == 0) || (missing > 0 && (a == 0 || b == 0)))
            return GT012_HOM_REF;
        if (a == 1 && b == 1) return GT012_HOM_ALT;
        if (a != b) return GT012_HET;
        return GT012_HOM_ALT;
    }
}

// gt012_pair_table[(r0 << 8) | r1] holds the genotype class of the diploid
// genotype (r0, r1): the non-strict class in the low nibble and the strict_gt
// class in the high nibble. gt_idx_table[r] holds the allele index of one raw
// int8 GT value. Branchy per-genotype classification is unpredictable on real
// data, so precomputing every (r0, r1) pair once is a large win.
static uint8_t gt012_pair_table[1 << 16];
static int32_t gt_idx_table[1 << 8];
static int gt012_tables_ready = 0;

static void gt012_tables_init(void) {
    int i0, i1;
    for (i0 = 0; i0 < 256; i0++) {
        int8_t r0 = (int8_t)i0;
        int32_t v0 = gt8_to_i32(r0);
        gt_idx_table[i0] = (v0 >= 0) ? bcf_gt_allele(v0) : v0;
        for (i1 = 0; i1 < 256; i1++) {
            int8_t g[2] = { r0, (int8_t)i1 };
            gt012_pair_table[(i0 << 8) | i1] =
                (uint8_t)(gt012_classify8(g, 2, 0) |
                          (gt012_classify8(g, 2, 1) << 4));
        }
    }
    gt012_tables_ready = 1;
}

// Row-writer used by the bulk gt_types reader: identical classification to
// gt_types_012_from_int8 below, but writes only the int8 0/1/2/3 codes (no
// allele indexes or phasing) straight into one row of the caller's matrix.
// Callers hold the GIL, which serializes the lazy table initialization.
int gt_types_012_row_from_int8(const int8_t *data, int num_samples, int ploidy,
                               int strict_gt, int HOM_ALT, int UNKNOWN,
                               int8_t *out) {
    int j = 0, i;
    const int ngts = num_samples * ploidy;
    int8_t class_map[4];
    class_map[GT012_HOM_REF] = 0;
    class_map[GT012_HET] = 1;
    class_map[GT012_HOM_ALT] = (int8_t)HOM_ALT;
    class_map[GT012_UNKNOWN] = (int8_t)UNKNOWN;

    if (!gt012_tables_ready) gt012_tables_init();

    if (ploidy == 2) {
        const int shift = strict_gt ? 4 : 0;
        for (i = 0, j = 0; j < num_samples; i += 2, j++) {
            unsigned r0 = (uint8_t)data[i], r1 = (uint8_t)data[i + 1];
            out[j] = class_map[(gt012_pair_table[(r0 << 8) | r1] >> shift) & 0xF];
        }
        return j;
    }

    for (i = 0; i < ngts; i += ploidy) {
        out[j++] = class_map[gt012_classify8(data + i, ploidy, strict_gt)];
    }
    return j;
}

// Single-pass version of the gt_types pipeline for the common case where
// htslib stores GT with per-value type BCF_BT_INT8 (fewer than ~64 alleles).
// It reads the packed int8 FORMAT data directly (no intermediate int32 copy)
// and fills, per sample:
//   gt_types  (num_samples)          the 0/1/2/3 genotype codes (as as_gts)
//   gt_idxs   (num_samples * ploidy) the allele index of each chromosome,
//                                    -1 for '.', negative sentinels preserved
//                                    as their int32 equivalents
//   gt_phased (num_samples)          1 when the sample's genotype is phased
// The results are identical to running bcf_get_genotypes + the allele-index
// extraction + as_gts, because the int8 -> int32 sentinel mapping preserves
// every comparison these steps make.
// Callers hold the GIL, which serializes the lazy table initialization.
int gt_types_012_from_int8(const int8_t *data, int num_samples, int ploidy,
                           int strict_gt, int HOM_ALT, int UNKNOWN,
                           int32_t *gt_types, int32_t *gt_idxs,
                           int32_t *gt_phased) {
    int j = 0, i, k;
    const int ngts = num_samples * ploidy;
    int32_t class_map[4];
    class_map[GT012_HOM_REF] = 0;
    class_map[GT012_HET] = 1;
    class_map[GT012_HOM_ALT] = HOM_ALT;
    class_map[GT012_UNKNOWN] = UNKNOWN;

    if (!gt012_tables_ready) gt012_tables_init();

    if (ploidy == 2) {
        const int shift = strict_gt ? 4 : 0;
        for (i = 0, j = 0; j < num_samples; i += 2, j++) {
            unsigned r0 = (uint8_t)data[i], r1 = (uint8_t)data[i + 1];
            gt_idxs[i] = gt_idx_table[r0];
            gt_idxs[i + 1] = gt_idx_table[r1];
            gt_phased[j] = (data[i] > 0) && (r1 & 1);
            gt_types[j] = class_map[(gt012_pair_table[(r0 << 8) | r1] >> shift) & 0xF];
        }
        return j;
    }

    for (i = 0; i < ngts; i += ploidy) {
        for (k = 0; k < ploidy; k++) {
            int32_t v = gt8_to_i32(data[i + k]);
            gt_idxs[i + k] = (v >= 0) ? bcf_gt_allele(v) : v;
        }
        gt_phased[j] = (data[i] > 0) && (i + 1 < ngts) && bcf_gt_is_phased(data[i + 1]);
        gt_types[j++] = class_map[gt012_classify8(data + i, ploidy, strict_gt)];
    }
    return j;
}

// use htslib's native lazy FORMAT parsing when this htslib provides it
// (bcf_hdr_set_parse_formats); otherwise report -2 so the caller falls back
// to stripping the line itself.
int cyvcf2_hdr_set_parse_formats(bcf_hdr_t *hdr, const char *fmts) {
#ifdef CYVCF2_HAVE_PARSE_FORMATS
    return bcf_hdr_set_parse_formats(hdr, fmts);
#else
    (void)hdr;
    (void)fmts;
    return -2;
#endif
}

// use htslib's threaded VCF text parsing when this htslib provides it
// (bcf_set_parse_threads); otherwise report -2 so the caller can warn.
#ifdef CYVCF2_HAVE_PARSE_CB
// Post-parse callback registered by VCF.gt_types_matrix: runs on htslib's
// parse worker threads and writes each record's 012 row straight into the
// preallocated matrix, so the conversion overlaps with parsing instead of
// serializing on the reading thread. Rows are addressed by the record's file
// ordinal; failures are flagged (not raised) and checked after the read.
static int gt012_matrix_cb(const bcf_hdr_t *h, bcf1_t *v, int64_t ordinal,
                           void *data) {
    gt012_sink_t *s = (gt012_sink_t *)data;
    int64_t row = ordinal - s->base;
    bcf_fmt_t *fmt;
    int8_t *out;
    if (row < 0) return 0; // delivered before the matrix read began
    if (row >= s->cap) {   // the index under-counted the file's records
        __atomic_fetch_add(&s->overflow, 1, __ATOMIC_RELAXED);
        return 0;
    }
    if (bcf_unpack(v, BCF_UN_FMT) < 0) {
        __atomic_fetch_add(&s->badrec, 1, __ATOMIC_RELAXED);
        return 0;
    }
    out = s->out + (size_t)row * s->n_samples;
    fmt = s->gt_fmt_id >= 0 ? bcf_get_fmt_id(v, s->gt_fmt_id) : NULL;
    if (fmt && fmt->type == BCF_BT_INT8 && fmt->n > 0 && fmt->p) {
        gt_types_012_row_from_int8((int8_t *)fmt->p, s->n_samples, fmt->n,
                                   s->strict_gt, s->hom_alt, s->unknown, out);
    } else {
        // GT stored wider than int8 (very many alleles), or absent
        int32_t *tmp = NULL;
        int ntmp = 0, i;
        int n = bcf_get_genotypes(h, v, &tmp, &ntmp);
        if (n <= 0 || n / s->n_samples == 0) {
            free(tmp);
            __atomic_fetch_add(&s->badrec, 1, __ATOMIC_RELAXED);
            return 0;
        }
        as_gts(tmp, s->n_samples, n / s->n_samples, s->strict_gt, s->hom_alt,
               s->unknown);
        for (i = 0; i < s->n_samples; i++) out[i] = (int8_t)tmp[i];
        free(tmp);
    }
    if (s->pos) s->pos[row] = v->pos + 1;
    return 0;
}
#endif

// virtual-offset split points for sharded reading (hts_idx_split);
// -2 when this htslib has no such API
int cyvcf2_idx_split(const hts_idx_t *idx, int nranges, uint64_t **starts,
                     int *nout) {
#ifdef CYVCF2_HAVE_IDX_SPLIT
    return hts_idx_split(idx, nranges, starts, nout);
#else
    (void)idx;
    (void)nranges;
    *starts = NULL;
    *nout = 0;
    return -2;
#endif
}

// register (or clear, sink == NULL) the matrix-fill callback on this
// reader's threaded parser; -2 when this htslib has no such API
int cyvcf2_set_gt012_sink(htsFile *fp, gt012_sink_t *sink) {
#ifdef CYVCF2_HAVE_PARSE_CB
    return bcf_set_parse_callback(fp, sink ? gt012_matrix_cb : NULL, sink);
#else
    (void)fp;
    (void)sink;
    return -2;
#endif
}

int cyvcf2_set_parse_threads(htsFile *fp, int n) {
#ifdef CYVCF2_HAVE_PARSE_THREADS
    return bcf_set_parse_threads(fp, n);
#else
    (void)fp;
    (void)n;
    return -2;
#endif
}

// does the comma-separated list `keep` contain the name at name[0..len)?
static int keep_has(const char *keep, const char *name, int len) {
    const char *p = keep;
    while (*p) {
        const char *q = p;
        while (*q && *q != ',') q++;
        if ((int)(q - p) == len && memcmp(p, name, len) == 0) return 1;
        p = (*q == ',') ? q + 1 : q;
    }
    return 0;
}

// matches htslib's MAX_N_FMT so the strip fallback and native lazy
// parsing handle the same envelope of records
#define CYVCF2_MAX_FMT 255

// Rewrite one VCF data line in place so that only the FORMAT fields named in
// `keep` (a comma-separated list, e.g. "GT" or "GT,DS") remain; the values of
// every other FORMAT field are removed from each sample column before the
// line reaches vcf_parse, so htslib never spends time converting them. Kept
// fields stay in their original order. Returns the new line length, or -1
// when the line is left unmodified: fewer than 10 columns, no sample columns,
// more than CYVCF2_MAX_FMT FORMAT fields, or none/all of the kept fields
// present in this line's FORMAT.
int vcf_line_strip_format(char *line, int len, const char *keep) {
    char *end = line + len;
    char *p = line, *q, *w, *r;
    char *fs[CYVCF2_MAX_FMT];
    int fl[CYVCF2_MAX_FMT];
    int kept[CYVCF2_MAX_FMT];
    int t, i, n_fmt = 0, n_kept = 0;

    for (t = 0; t < 8; t++) { // skip CHROM..INFO; FORMAT is column 9
        p = memchr(p, '\t', end - p);
        if (p == NULL) return -1;
        p++;
    }
    char *fmt_start = p;
    char *fmt_end = memchr(p, '\t', end - p);
    if (fmt_end == NULL) return -1; // no sample columns

    q = fmt_start;
    while (q <= fmt_end) {
        if (n_fmt == CYVCF2_MAX_FMT) return -1;
        char *e = memchr(q, ':', fmt_end - q);
        if (e == NULL) e = fmt_end;
        fs[n_fmt] = q;
        fl[n_fmt] = (int)(e - q);
        n_fmt++;
        if (e == fmt_end) break;
        q = e + 1;
    }
    for (i = 0; i < n_fmt; i++) {
        if (keep_has(keep, fs[i], fl[i])) kept[n_kept++] = i;
    }
    if (n_kept == n_fmt) return -1;

    w = fmt_start;
    if (n_kept == 0) {
        // none of the kept fields is present: drop the FORMAT data entirely
        // ("." columns), matching bcf_hdr_set_parse_formats() semantics
        *w++ = '.';
        r = fmt_end;
        while (r < end) {
            *w++ = '\t';
            r++;
            char *cell_end = memchr(r, '\t', end - r);
            if (cell_end == NULL) cell_end = end;
            *w++ = '.';
            r = cell_end;
        }
        *w = '\0';
        return (int)(w - line);
    }
    for (i = 0; i < n_kept; i++) {
        if (i) *w++ = ':';
        memmove(w, fs[kept[i]], fl[kept[i]]);
        w += fl[kept[i]];
    }
    r = fmt_end;
    while (r < end) {
        // r is at the '\t' introducing the next sample column
        *w++ = '\t';
        r++;
        char *cell_end = memchr(r, '\t', end - r);
        if (cell_end == NULL) cell_end = end;
        char *ts[CYVCF2_MAX_FMT];
        int tl[CYVCF2_MAX_FMT];
        int nt = 0;
        q = r;
        while (q <= cell_end && nt < n_fmt) {
            // any surplus colons are folded into the last field; such a line
            // is malformed and fails in vcf_parse either way
            char *e = (nt == n_fmt - 1) ? cell_end : memchr(q, ':', cell_end - q);
            if (e == NULL) e = cell_end;
            ts[nt] = q;
            tl[nt] = (int)(e - q);
            nt++;
            if (e == cell_end) break;
            q = e + 1;
        }
        int emitted = 0;
        for (i = 0; i < n_kept; i++) {
            if (kept[i] >= nt) break; // cell dropped its trailing fields
            if (emitted) *w++ = ':';
            memmove(w, ts[kept[i]], tl[kept[i]]);
            w += tl[kept[i]];
            emitted = 1;
        }
        if (!emitted) *w++ = '.';
        r = cell_end;
    }
    *w = '\0';
    return (int)(w - line);
}

int as_gts(int32_t *gts, int num_samples, int ploidy, int strict_gt, int HOM_ALT, int UNKNOWN) {
    int j = 0, i, k;
    int missing= 0;
    for (i = 0; i < ploidy * num_samples; i += ploidy){
        missing = 0;
        for (k = 0; k < ploidy; k++) {
            if bcf_gt_is_missing(gts[i+k])  {
                missing += 1;
            }
        }
        if (missing == ploidy) {
            gts[j++] = UNKNOWN; // unknown
            continue;
        } else if ( (missing != 0) && (strict_gt == 1) ) {
            gts[j++] = UNKNOWN; // unknown
            continue;
        }

        if(ploidy == 1 || gts[i+1] == bcf_int32_vector_end) {
            int a = bcf_gt_allele(gts[i]);
            if (a == 0) {
                   gts[j++] = 0;
            } else if (a == 1) {
                gts[j++] = HOM_ALT;
            } else {
                gts[j++] = UNKNOWN;
            }
            continue;
        }

        int a = bcf_gt_allele(gts[i]);
        int b = bcf_gt_allele(gts[i+1]);

        if((a == 0) && (b == 0)) {
            gts[j++] = 0; //  HOM_REF
            continue;
        }
        //fprintf(stderr, "i: %d\tmissing:%d\ta:%d\tb:%d\n", i/ploidy, missing, a, b);
        if ((missing > 0) && ((a == 0) || (b == 0))) {
            // if a single allele is missing e.g 0/. it's still encoded as hom ref because it has no alts
            gts[j++] = 0; // HOM_REF
            continue;
        }
        else if((a == 1) && (b == 1)) {
            gts[j] = HOM_ALT; //  HOM_ALT
        }
        else if((a != b)) {
            gts[j] = 1; //  HET
        }
        else if((a == b)) {
            gts[j] = HOM_ALT; //  HOM_ALT
        } else {
            gts[j] = UNKNOWN; // unknown
        }
        j++;
    }
    return j;
}

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

// this is taken directly from atks/vt
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq)
{
    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    int32_t *lens = (int32_t*) malloc(m*sizeof(int32_t));
    khint_t k;
    int found = 0;

    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        lens[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if (lens[tid] > 0 && sscanf(kh_val(d, k).hrec[0]->vals[lens[tid]],"%d",&j) )
            lens[tid] = j;
    if(lens[tid] > 0){
      found++;
    }
    }
    *nseq = m;
    // found is used to check that we actually got the lengths.
    if(found == 0){
      *nseq = -1;
    }
    return lens;
}


