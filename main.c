#include "fft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SAMPLE_COUNT 64

static float complex samples[SAMPLE_COUNT];

void reset_signal() {
    int i;

    for (i = 0; i < SAMPLE_COUNT; i++)
        samples[i] = 0.f;
}

void add_signal(float freq, float amplitude) {
    int i;
    float d_angle;

    d_angle = 2.f * freq * (float)M_PI / SAMPLE_COUNT;

    for (i = 0; i < SAMPLE_COUNT; i++)
        samples[i] += amplitude * cosf(i * d_angle);
}

void envelope() {
    int i;
    float d_angle;

    d_angle = (float)M_PI / SAMPLE_COUNT;

    for (i = 0; i < SAMPLE_COUNT; i++) {
        float value = sinf(i * d_angle);
        samples[i] *= value * value;
    }
}

int main() {
    sfft_ctx_t fft;
    int i;

    reset_signal();
    add_signal(12.1f, 1.f);
    envelope();

    sfft_init(&fft, SAMPLE_COUNT);
    sfft_rad2_dif(&fft, (float complex *)&samples);

    for (i = 0; i < SAMPLE_COUNT / 2; i++)
        printf("%.1f ", cabsf(samples[sfft_reverse_idx(&fft, i)]));

    printf("\n");

    sfft_deinit(&fft);
    return EXIT_SUCCESS;
}