#ifndef _FFT_FACTORS_8_32B_INT_H_
#define _FFT_FACTORS_8_32B_INT_H_

#include <stdint.h>
#include "fxp.h"
#include "defines.h"

/*
 * Fixed-point format (X-HEEP CGRA complex):
 * 32-bit signed, Q16.15:
 * - 1 sign bit
 * - 16 integer bits
 * - 15 fractional bits
 *
 * Twiddles for W8^k = exp(-j*2*pi*k/8), k=0..3
 */
const fxp f_real[FFT_SIZE/2] = {
   32768,   //  1.0
   23170,   //  cos(pi/4)
       0,   //  cos(pi/2)
  -23170    //  cos(3pi/4)
};

const fxp f_imag[FFT_SIZE/2] = {
       0,   // -sin(0)
  -23170,   // -sin(pi/4)
  -32768,   // -sin(pi/2)
  -23170    // -sin(3pi/4)
};

#endif // _FFT_FACTORS_8_32B_INT_H_