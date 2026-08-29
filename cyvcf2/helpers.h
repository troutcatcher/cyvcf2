#include <htslib/vcf.h>
#include <htslib/hts_log.h>
#include <htslib/khash.h>

int as_gts(int *gts, int num_samples, int ploidy, int strict_gt, int HOM_ALT, int UNKNOWN);
int gt_types_012_from_int8(const int8_t *data, int num_samples, int ploidy,
                           int strict_gt, int HOM_ALT, int UNKNOWN,
                           int32_t *gt_types, int32_t *gt_idxs,
                           int32_t *gt_phased);
int gt_types_012_row_from_int8(const int8_t *data, int num_samples, int ploidy,
                               int strict_gt, int HOM_ALT, int UNKNOWN,
                               int8_t *out);
int vcf_line_strip_format(char *line, int len, const char *keep);
int cyvcf2_hdr_set_parse_formats(bcf_hdr_t *hdr, const char *fmts);
int cyvcf2_set_parse_threads(htsFile *fp, int n);

// destination and settings for the worker-side 012 matrix fill
typedef struct {
    int8_t *out;        // (cap x n_samples) row-major matrix
    int64_t *pos;       // optional POS output (cap entries), or NULL
    int64_t base;       // file ordinal of the matrix's first row
    long cap;           // rows allocated
    int n_samples;
    int strict_gt, hom_alt, unknown;
    int gt_fmt_id;      // header id of FORMAT/GT, or -1
    long overflow;      // records beyond cap (index under-counted)
    long badrec;        // records whose genotypes could not be read
} gt012_sink_t;

int cyvcf2_set_gt012_sink(htsFile *fp, gt012_sink_t *sink);
int cyvcf2_idx_split(const hts_idx_t *idx, int nranges, uint64_t **starts,
                     int *nout);
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq);
