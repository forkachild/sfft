#ifndef _FFT_H_
#define _FFT_H_

#include <complex.h>

/**
 * @brief Single-precision context which needs to be initiated for a constant
 * input-size.
 *
 * @param size Number of signals in the input (Must be a power of 2)
 * @param reversed_indices `size` sized cache with each element representing
 * its bit-reversed value
 * @param twiddles Complex constants used to calculate the FFT
 */
typedef struct sfft_ctx {
    unsigned int size;
    unsigned int *reversed_indices;
    float complex *twiddles;
} sfft_ctx_t;

/**
 * @brief Double-precision context which needs to be initiated for a constant
 * input-size.
 *
 * @param size Number of signals in the input (Must be a power of 2)
 * @param reversed_indices `size` sized cache with each element representing
 * its bit-reversed value
 * @param twiddles Complex constants used to calculate the FFT
 */
typedef struct sfft_ctx_d {
    unsigned int size;
    unsigned int *reversed_indices;
    double complex *twiddles;
} sfft_ctx_d_t;

/**
 * @brief Initialize a single-precision context.
 *
 * @param ctx An un-initialized single-precision context
 * @param size Number of input signals to initialize for
 */
void sfft_init(sfft_ctx_t *ctx, unsigned int size);

/**
 * @brief Initialize a double-precision context.
 *
 * @param ctx An un-initialized single-precision context
 * @param size Number of input signals to initialize for
 */
void sfft_init_d(sfft_ctx_d_t *ctx, unsigned int size);

/**
 * @brief Single-precision decimation-in-time in-place FFT. Use this when the
 * access-locality of the result is of more priority than the computation speed.
 *
 * @param ctx An initialized single-precision context
 * @param x Array of input signals. Must be the same size as the one which was
 * used to initialize the context. After execution, this array will contain the
 * result, accessible in-order.
 */
void sfft_rad2_dit(const sfft_ctx_t *ctx, float complex *x);

/**
 * @brief Double-precision decimation-in-time in-place FFT. Use this when the
 * access-locality of the result is of more priority than the computation speed.
 *
 * @param ctx An initialized double-precision context
 * @param x Array of input signals. Must be the same size as the one which was
 * used to initialize the context. After execution, this array will contain the
 * result, accessible in-order.
 */
void sfft_rad2_dit_d(const sfft_ctx_d_t *ctx, double complex *x);

/**
 * @brief Single-precision decimation-in-frequency in-place FFT. Use this when
 * the access-locality of the result is of less priority than the computation
 * speed.
 *
 * @param ctx An initialized single-precision context
 * @param x Array of input signals. Must be the same size as the one which was
 * used to initialize the context. After execution, this array will contain the
 * result, accessible in bit-reversed order through `sfft_reverse_idx()`
 * function.
 */
void sfft_rad2_dif(const sfft_ctx_t *ctx, float complex *x);

/**
 * @brief Double-precision decimation-in-frequency in-place FFT. Use this when
 * the access-locality of the result is of less priority than the computation
 * speed.
 *
 * @param ctx An initialized double-precision context
 * @param x Array of input signals. Must be the same size as the one which was
 * used to initialize the context. After execution, this array will contain the
 * result, accessible in bit-reversed order through `sfft_reverse_idx()`
 * function.
 */
void sfft_rad2_dif_d(const sfft_ctx_d_t *ctx, double complex *x);

/**
 * @brief Generate the bit-reversed index from a given index relevant
 * to the initialized context
 *
 * @param ctx An initialized single-precision context
 * @param idx The index to reverse
 * @returns Bit-reversed index
 */
unsigned int sfft_reverse_idx(const sfft_ctx_t *ctx, unsigned int idx);

/**
 * @brief Generate the bit-reversed index from a given index relevant
 * to the initialized context
 *
 * @param ctx An initialized double-precision context
 * @param idx The index to reverse
 * @returns Bit-reversed index
 */
unsigned int sfft_reverse_idx_d(const sfft_ctx_d_t *ctx, unsigned int idx);

/**
 * @brief De-initialize a single-precision context
 *
 * @param ctx An initialized single-precision context
 */
void sfft_deinit(sfft_ctx_t *ctx);

/**
 * @brief De-initialize a double-precision context
 *
 * @param ctx An initialized double-precision context
 */
void sfft_deinit_d(sfft_ctx_d_t *ctx);

#endif