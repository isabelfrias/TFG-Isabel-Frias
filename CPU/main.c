/*
 * CPU version of the FFT algorithm
 */

#include <stdio.h>
#include <stdlib.h>

#include "csr.h"
#include "hart.h"
#include "handler.h"
#include "core_v_mini_mcu.h"
#include "rv_plic.h"
#include "rv_plic_regs.h"
#include "heepsilon.h"
#include "fxp.h"
#include "defines.h"
#include "fft_data.h"

#if FFT_SIZE == 8
  #include "fft_factors_8_32b_int.h"
#endif
#if FFT_SIZE==256
  #include "fft_factors_256_32b_int.h"
#endif
#if FFT_SIZE==512
  #include "fft_factors_512_32b_int.h"
#endif
#if FFT_SIZE==1024
  #include "fft_factors_1024_32b_int.h"
#endif
#if FFT_SIZE==2048
  #include "fft_factors_2048_32b_int.h"
#endif

#ifndef CPU_FREQ
#define CPU_FREQ 100000000 // 100000000 Hz (100000 kHz)
#endif

/* --------------------------------------------------------------------------
 *                     Global variables
 * --------------------------------------------------------------------------*/

// FFT radix-2 variables
fxp in_X_r[FFT_SIZE] __attribute__ ((aligned (4))) = { 0 };
fxp in_X_i[FFT_SIZE] __attribute__ ((aligned (4))) = { 0 };

/* --------------------------------------------------------------------------
 *                     Auxiliary functions
 * --------------------------------------------------------------------------*/

 // Function to reverse the bits of an index (used for bit-reversal in FFT)
 uint16_t ReverseBits(uint16_t index, uint16_t numBits) {
  uint16_t i, rev;

  for (i=rev=0; i<numBits; i++) {
    rev = (rev << 1) | (index & 1);
    index >>= 1;
  }

  return rev;
}

// Returns the number of bits needed to represent a power of two (e.g., 8 -> 3, 16 -> 4, etc.)
uint16_t NumberOfBitsNeeded(uint16_t powerOfTwo) {
  uint16_t i;

  if (powerOfTwo < 2) {
   return 0; // should not happen
  }

  for (i=0;; i++) {
    if (powerOfTwo & (1 << i))
      return i;
  }
}

/* --------------------------------------------------------------------------
 *                     Main function
 * --------------------------------------------------------------------------*/

int main(void) {
  printf("=== FFT RISC-V===\n");

  uint32_t prog_cycles_start = 0;
  uint32_t prog_cycles_end = 0;
  uint32_t bitrev_cycles_start = 0;
  uint32_t bitrev_cycles_end = 0;
  uint32_t fft_cycles_start = 0;
  uint32_t fft_cycles_end = 0;

  CSR_READ(CSR_REG_MCYCLE, &prog_cycles_start);

  //Bit reverse------------------------------------------
  uint16_t numBits = NumberOfBitsNeeded ( FFT_SIZE );

  // Load input data into the FFT input buffers (Real and Imaginary parts in separate arrays)
  for (int i = 0; i < FFT_SIZE; i++) {
    in_X_r[i] = input_signal[2*i + 0];
    in_X_i[i] = input_signal[2*i + 1];
  }

  CSR_READ(CSR_REG_MCYCLE, &bitrev_cycles_start);
  // Bit-reversal in-place for DIT FFT input ordering
  for (int i = 0; i < FFT_SIZE; i++) {
    uint16_t j = ReverseBits((uint16_t)i, numBits);
    if (j > i) {
      fxp tmp;
      tmp = in_X_r[i];
      in_X_r[i] = in_X_r[j];
      in_X_r[j] = tmp;

      tmp = in_X_i[i];
      in_X_i[i] = in_X_i[j];
      in_X_i[j] = tmp;
    }
  }
  
  CSR_READ(CSR_REG_MCYCLE, &bitrev_cycles_end);
  printf("Bitreverse done\n");
  printf("Input loaded for FFT_SIZE=%d\n", FFT_SIZE);

  //FFT Radix-2 Computation-------------------------------------
  CSR_READ(CSR_REG_MCYCLE, &fft_cycles_start);
  for(int i = 0; i < numBits; i++){
    uint16_t L = 1 << i; // Number of blocks in this stage
    
    for(int j = 0; j < FFT_SIZE; j += 2*L) { // For each block
      for(int k = 0; k < L; k++) { // For each butterfly in the block
        uint16_t index1 = j + k; // top index
        uint16_t index2 = j + k + L; //bottom index
        uint16_t tw_idx = FFT_SIZE * k / 2 / L; // Twiddle factor index
        
        fxp product_r = fxp_mult(in_X_r[index2], f_real[tw_idx]) - fxp_mult(in_X_i[index2], f_imag[tw_idx]);
        fxp product_i = fxp_mult(in_X_r[index2], f_imag[tw_idx]) + fxp_mult(in_X_i[index2], f_real[tw_idx]);

        fxp top_r, bottom_r, top_i, bottom_i;
        top_r = product_r + in_X_r[index1];
        top_i = product_i + in_X_i[index1];

        bottom_r = in_X_r[index1] - product_r;
        bottom_i = in_X_i[index1] - product_i;

        in_X_r[index1] = top_r;
        in_X_i[index1] = top_i;

        in_X_r[index2] = bottom_r;
        in_X_i[index2] = bottom_i;
      }
    }
  }
  
  CSR_READ(CSR_REG_MCYCLE, &fft_cycles_end);
  printf("FFT done\n");
  CSR_READ(CSR_REG_MCYCLE, &prog_cycles_end);


  /* -------------------------------------------------------------------------------
  *
  * CHECK RESULTS
  * 
  ------------------------------------------------------------------------------- */

  int32_t errors=0;
  for (int i=0; i<FFT_SIZE; i++) {
    if(in_X_r[i] != exp_output_real[i] ||
        in_X_i[i] != exp_output_imag[i]) {
          errors++;
          if (errors <= 10){
            printf("Error at index %d: got (%d, %d), expected (%d, %d)\n", i, in_X_r[i], in_X_i[i], exp_output_real[i], exp_output_imag[i]);
          }
    }
  }

  /* -------------------------------------------------------------------------------
  *
  * Performance Counters and Metrics
  * 
  ------------------------------------------------------------------------------- */

  uint32_t fft_cycles_total = fft_cycles_end - fft_cycles_start;
  uint32_t bitrev_cycles_total = bitrev_cycles_end - bitrev_cycles_start;
  uint32_t prog_cycles_total = prog_cycles_end - prog_cycles_start;

  unsigned long bitrev_time_us = (unsigned long)(((uint64_t)bitrev_cycles_total * 1000000ULL) / (uint64_t)CPU_FREQ);
  unsigned long fft_cycles_time_cpu_us = (unsigned long)(((uint64_t)fft_cycles_total * 1000000ULL) / (uint64_t)CPU_FREQ);
  unsigned long prog_time_us = (unsigned long)(((uint64_t)prog_cycles_total * 1000000ULL) / (uint64_t)CPU_FREQ);

  printf("Bit-reversal total cycles (mcycle, CPU only): %lu (%lu us)\n",
    (unsigned long)bitrev_cycles_total,
    bitrev_time_us);
  printf("FFT section total cycles (mcycle, CPU only): %lu (%lu us)\n",
    (unsigned long)fft_cycles_total,
    fft_cycles_time_cpu_us);
  printf("Program total cycles (mcycle, CPU only): %lu (%lu us)\n",
    (unsigned long)prog_cycles_total,
    prog_time_us);

  double power = 2.5e-3; //Paper reference: 2.5 mW (2.5e-3 W)
  double fft_time = (double)fft_cycles_total / (double)CPU_FREQ;
  double fft_energy = power * fft_time;
  uint64_t fft_energy_nj = (uint64_t)(fft_energy * 1e9);

  printf("FFT energy: %lu nJ\n", (unsigned long)fft_energy_nj);

  double prog_time = (double)prog_cycles_total / (double)CPU_FREQ;
  double prog_energy = power * prog_time;
  uint64_t prog_energy_nj = (uint64_t)(prog_energy * 1e9);
  
  printf("Program energy: %lu nJ\n", (unsigned long)prog_energy_nj);

  return 0;
}
 
