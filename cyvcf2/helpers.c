#include <helpers.h>
#include <string.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

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

// matches htslib's MAX_N_FMT
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
        // ("." columns)
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


