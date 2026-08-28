#include <helpers.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

// map one raw int8 GT value to the int32 value bcf_get_genotypes would produce
static inline int32_t gt8_to_i32(int8_t v) {
    if (v == bcf_int8_missing) return bcf_int32_missing;
    if (v == bcf_int8_vector_end) return bcf_int32_vector_end;
    return v;
}

// genotype classes before HOM_ALT/UNKNOWN are mapped for gts012
enum { GT012_HOM_REF = 0, GT012_HET = 1, GT012_HOM_ALT = 2, GT012_UNKNOWN = 3 };

// classify one diploid genotype (the two raw int8 GT values); this mirrors
// as_gts exactly, evaluated on the int32 values bcf_get_genotypes produces
static int gt012_classify_pair(int8_t r0, int8_t r1, int strict_gt) {
    int32_t v0 = gt8_to_i32(r0), v1 = gt8_to_i32(r1);
    int missing = bcf_gt_is_missing(v0) + bcf_gt_is_missing(v1);
    if (missing == 2) return GT012_UNKNOWN;
    if (missing != 0 && strict_gt) return GT012_UNKNOWN;
    if (v1 == bcf_int32_vector_end) {
        int a = bcf_gt_allele(v0);
        if (a == 0) return GT012_HOM_REF;
        if (a == 1) return GT012_HOM_ALT;
        return GT012_UNKNOWN;
    }
    int a = bcf_gt_allele(v0);
    int b = bcf_gt_allele(v1);
    if (a == 0 && b == 0) return GT012_HOM_REF;
    // if a single allele is missing e.g 0/. it's still encoded as hom ref
    // because it has no alts
    if (missing > 0 && (a == 0 || b == 0)) return GT012_HOM_REF;
    if (a == 1 && b == 1) return GT012_HOM_ALT;
    if (a != b) return GT012_HET;
    return GT012_HOM_ALT;
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
            int8_t r1 = (int8_t)i1;
            gt012_pair_table[(i0 << 8) | i1] =
                (uint8_t)(gt012_classify_pair(r0, r1, 0) |
                          (gt012_classify_pair(r0, r1, 1) << 4));
        }
    }
    gt012_tables_ready = 1;
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

    if (!gt012_tables_ready) gt012_tables_init();

    if (ploidy == 2) {
        const int shift = strict_gt ? 4 : 0;
        int32_t class_map[4];
        class_map[GT012_HOM_REF] = 0;
        class_map[GT012_HET] = 1;
        class_map[GT012_HOM_ALT] = HOM_ALT;
        class_map[GT012_UNKNOWN] = UNKNOWN;
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
        int missing = 0;
        for (k = 0; k < ploidy; k++) {
            int32_t v = gt8_to_i32(data[i + k]);
            gt_idxs[i + k] = (v >= 0) ? bcf_gt_allele(v) : v;
            if (bcf_gt_is_missing(v)) {
                missing += 1;
            }
        }
        gt_phased[j] = (data[i] > 0) && (i + 1 < ngts) && bcf_gt_is_phased(data[i + 1]);

        if (missing == ploidy || (missing != 0 && strict_gt == 1)) {
            gt_types[j++] = UNKNOWN;
            continue;
        }

        if (ploidy == 1 || data[i + 1] == bcf_int8_vector_end) {
            int a = bcf_gt_allele(data[i]);
            if (a == 0) {
                gt_types[j++] = 0;
            } else if (a == 1) {
                gt_types[j++] = HOM_ALT;
            } else {
                gt_types[j++] = UNKNOWN;
            }
            continue;
        }

        int a = bcf_gt_allele(data[i]);
        int b = bcf_gt_allele(data[i + 1]);

        if ((a == 0) && (b == 0)) {
            gt_types[j++] = 0; // HOM_REF
        } else if ((missing > 0) && ((a == 0) || (b == 0))) {
            // if a single allele is missing e.g 0/. it's still encoded as hom
            // ref because it has no alts
            gt_types[j++] = 0; // HOM_REF
        } else if ((a == 1) && (b == 1)) {
            gt_types[j++] = HOM_ALT;
        } else if (a != b) {
            gt_types[j++] = 1; // HET
        } else {
            gt_types[j++] = HOM_ALT;
        }
    }
    return j;
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


