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
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq);
