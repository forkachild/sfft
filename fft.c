#include "fft.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

/**
 * @brief O(log n) Ultra fast log base-2 of only 2^n numbers. For others, the
 * result/behavior is invalid/undefined.
 *
 * Since 2^n numbers will have only one '1' bit, we just need to shift
 * right until we find it, and that's log-base-2.
 *
 * @param n The input. *MUST BE A POWER OF 2*
 * @return The log-base-2 result
 */
static inline unsigned int fast_log_2(unsigned int n) {
    unsigned int shift_count, i;

    // Keep shifting right until we find a '1' at the LSB
    // and that's when we know we hit the jackpot!
    for (shift_count = 0, i = n; (i & 0b1) == 0; i >>= 1, shift_count++)
        ;

    return shift_count;
}

/**
 * @brief O(n) order reverse bits.
 *
 * @param n Value to be bit-reversed
 * @param bit_depth Number of bits to be reversed
 * @return `n` bit-reversed
 */
static inline unsigned int reverse_bits(unsigned int n,
                                        unsigned int bit_depth) {
    unsigned int output, i;

    // Simple, left shift one, right shift the other, bleh!
    for (output = 0, i = 0; i < bit_depth; i++, n >>= 1)
        output = (output << 1) | (n & 0b1);

    return output;
}

static inline void fill_reversed_indices(unsigned int *reversed_indices,
                                         unsigned int size) {
    unsigned int bit_depth, i;

    // Number of bits required
    bit_depth = fast_log_2(size);

    for (i = 0; i < size; i++)
        reversed_indices[i] = reverse_bits(i, bit_depth);
}

static inline void fill_twiddles(float complex *twiddles, unsigned int size) {
    float angle_per_sample;
    unsigned int half_size, i;

    // Mr. Clean
    half_size = size / 2;

    // -2pi/N, constant, reducing the calculations
    angle_per_sample = -(float)M_PI / half_size;

    // Cache the twiddle factors
    // Compromise some space for HUGE performance gain
    // Cache locality baby!
    for (i = 0; i < half_size; i++)
        twiddles[i] = cexp(angle_per_sample * i * I);
}

static inline void fill_twiddles_d(double complex *twiddles,
                                   unsigned int size) {
    double angle_per_sample;
    unsigned int half_size, i;

    // Mr. Clean
    half_size = size / 2;

    // -2pi/N, constant, reducing the calculations
    angle_per_sample = -M_PI / half_size;

    // Cache the twiddle factors
    // Compromise some space for HUGE performance gain
    // Cache locality baby!
    for (i = 0; i < half_size; i++)
        twiddles[i] = cexp(angle_per_sample * i * I);
}

void sfft_init(sfft_ctx_t *ctx, unsigned int size) {
    unsigned int *reversed_indices;
    float complex *twiddles;

    reversed_indices = malloc(size * sizeof(unsigned int));
    if (reversed_indices == NULL)
        return;

    twiddles = malloc((size / 2) * sizeof(float complex));
    if (twiddles == NULL) {
        free(reversed_indices);
        return;
    }

    fill_reversed_indices(reversed_indices, size);
    fill_twiddles(twiddles, size);

    ctx->size = size;
    ctx->reversed_indices = reversed_indices;
    ctx->twiddles = twiddles;
}

void sfft_init_d(sfft_ctx_d_t *ctx, unsigned int size) {
    unsigned int *reversed_indices;
    double complex *twiddles;

    if (ctx->reversed_indices != NULL || ctx->twiddles != NULL)
        return;

    reversed_indices = malloc(size * sizeof(unsigned int));
    if (reversed_indices == NULL)
        return;

    twiddles = malloc((size / 2) * sizeof(double complex));
    if (twiddles == NULL) {
        free(reversed_indices);
        return;
    }

    fill_reversed_indices(reversed_indices, size);
    fill_twiddles_d(twiddles, size);

    ctx->size = size;
    ctx->reversed_indices = reversed_indices;
    ctx->twiddles = twiddles;
}

void sfft_rad2_dit(const sfft_ctx_t *ctx, float complex *x) {
    unsigned int half_size, set_count, ops_per_set, set, start, butterfly,
        butterfly_top_idx, butterfly_bottom_idx;
    float complex twiddle, butterfly_top, butterfly_bottom;

    // Mr. Clean
    half_size = ctx->size / 2;

    // Perform the stages
    // i is the number of sets to perform
    // eg For N = 8
    // i = 4,2,1
    for (set_count = half_size; set_count >= 1; set_count >>= 1) {
        // No of operations per set in this stage
        // ops_per_set = 1,2,4
        ops_per_set = half_size / set_count;

        // Loop over sets
        // j is the set #
        for (set = 0; set < set_count; set++) {
            // Start the butterflies
            start = set * ops_per_set * 2;

            // Loop over butterflies
            for (butterfly = 0; butterfly < ops_per_set; butterfly++) {
                butterfly_top_idx = ctx->reversed_indices[start + butterfly];
                butterfly_bottom_idx =
                    ctx->reversed_indices[start + butterfly + ops_per_set];

                // Determine the twiddle to pre-multiply the lower
                // half of the butterfly
                twiddle =
                    ctx->twiddles[butterfly * set_count]; // Cache hit baby!
                butterfly_top = x[butterfly_top_idx];
                butterfly_bottom = twiddle * x[butterfly_bottom_idx];

                // Finally the butterfly
                x[butterfly_top_idx] = butterfly_top + butterfly_bottom;
                x[butterfly_bottom_idx] = butterfly_top - butterfly_bottom;
            }
        }
    }
}

void sfft_rad2_dit_d(const sfft_ctx_d_t *ctx, double complex *x) {
    unsigned int half_size, set_count, ops_per_set, set, start, butterfly,
        butterfly_top_idx, butterfly_bottom_idx;
    double complex twiddle, butterfly_top, butterfly_bottom;

    // Mr. Clean
    half_size = ctx->size / 2;

    // Perform the stages
    // i is the number of sets to perform
    // eg For N = 8
    // i = 4,2,1
    for (set_count = half_size; set_count >= 1; set_count >>= 1) {
        // No of operations per set in this stage
        // ops_per_set = 1,2,4
        ops_per_set = half_size / set_count;

        // Loop over sets
        // j is the set #
        for (set = 0; set < set_count; set++) {
            // Start the butterflies
            start = set * ops_per_set * 2;

            // Loop over butterflies
            for (butterfly = 0; butterfly < ops_per_set; butterfly++) {
                butterfly_top_idx = ctx->reversed_indices[start + butterfly];
                butterfly_bottom_idx =
                    ctx->reversed_indices[start + butterfly + ops_per_set];

                // Determine the twiddle to pre-multiply the lower
                // half of the butterfly
                twiddle =
                    ctx->twiddles[butterfly * set_count]; // Cache hit baby!
                butterfly_top = x[butterfly_top_idx];
                butterfly_bottom = twiddle * x[butterfly_bottom_idx];

                // Finally the butterfly
                x[butterfly_top_idx] = butterfly_top + butterfly_bottom;
                x[butterfly_bottom_idx] = butterfly_top - butterfly_bottom;
            }
        }
    }
}

void sfft_rad2_dif(const sfft_ctx_t *ctx, float complex *x) {
    unsigned int half_size, set_count, ops_per_set, set, start, butterfly,
        butterfly_top_idx, butterfly_bottom_idx;
    float complex twiddle, butterfly_top, butterfly_bottom;

    // Mr. Clean!
    half_size = ctx->size / 2;

    // Perform the stages
    // i is the number of sets to perform
    // eg For N = 8
    // i = 1,2,4
    for (set_count = 1; set_count <= half_size; set_count <<= 1) {
        // No of operations per set in this stage
        // ops_per_set = 4,2,1
        ops_per_set = half_size / set_count;

        // Loop over sets
        for (set = 0; set < set_count; set++) {
            // Start the butterflies
            start = set * ops_per_set * 2;

            // Loop over butterflies
            for (butterfly = 0; butterfly < ops_per_set; butterfly++) {
                butterfly_top_idx = start + butterfly;
                butterfly_bottom_idx = start + butterfly + ops_per_set;

                // Determine the twiddle to pre-multiply the lower
                // half of the butterfly
                // Cache hit baby!
                twiddle = ctx->twiddles[butterfly * set_count];
                butterfly_top = x[butterfly_top_idx];
                butterfly_bottom = x[butterfly_bottom_idx];

                // Finally the butterfly
                x[butterfly_top_idx] = butterfly_top + butterfly_bottom;
                x[butterfly_bottom_idx] =
                    twiddle * (butterfly_top - butterfly_bottom);
            }
        }
    }
}

void sfft_rad2_dif_d(const sfft_ctx_d_t *ctx, double complex *x) {
    unsigned int half_size, set_count, ops_per_set, set, start, butterfly,
        butterfly_top_idx, butterfly_bottom_idx;
    double complex twiddle, butterfly_top, butterfly_bottom;

    // Mr. Clean
    half_size = ctx->size / 2;

    // Perform the stages
    // i is the number of sets to perform
    // eg For N = 8
    // i = 1,2,4
    for (set_count = 1; set_count <= half_size; set_count <<= 1) {
        // No of operations per set in this stage
        // ops_per_set = 4,2,1
        ops_per_set = half_size / set_count;

        // Loop over sets
        for (set = 0; set < set_count; set++) {
            // Start the butterflies
            start = set * ops_per_set * 2;

            // Loop over butterflies
            for (butterfly = 0; butterfly < ops_per_set; butterfly++) {
                butterfly_top_idx = start + butterfly;
                butterfly_bottom_idx = start + butterfly + ops_per_set;

                // Determine the twiddle to pre-multiply the lower
                // half of the butterfly
                // Cache hit baby!
                twiddle = ctx->twiddles[butterfly * set_count];
                butterfly_top = x[butterfly_top_idx];
                butterfly_bottom = x[butterfly_bottom_idx];

                // Finally the butterfly
                x[butterfly_top_idx] = butterfly_top + butterfly_bottom;
                x[butterfly_bottom_idx] =
                    twiddle * (butterfly_top - butterfly_bottom);
            }
        }
    }
}

unsigned int sfft_reverse_idx(const sfft_ctx_t *ctx, unsigned int idx) {
    return ctx->reversed_indices[idx];
}

unsigned int sfft_reverse_idx_d(const sfft_ctx_d_t *ctx, unsigned int idx) {
    return ctx->reversed_indices[idx];
}

void sfft_deinit(sfft_ctx_t *ctx) {
    free((void *)ctx->reversed_indices);
    free((void *)ctx->twiddles);

    ctx->size = 0;
    ctx->reversed_indices = NULL;
    ctx->twiddles = NULL;
}

void sfft_deinit_d(sfft_ctx_d_t *ctx) {
    free((void *)ctx->reversed_indices);
    free((void *)ctx->twiddles);

    ctx->size = 0;
    ctx->reversed_indices = NULL;
    ctx->twiddles = NULL;
}