/*
 * CGRA_O1
 * FFT implementation (First optimization)
 * In this version, the bus is configured in ntoM format, which allows concurrent memory accesses from the CGRA, reducing stalls and improving performance compared to the onetoM configuration.
 * In addition, the CGRA computes all the butterflies in one stage increasing granularity, which avoids multiple kernel launches to the CGRA and continuous reconfiguration of write and read pointers in each iteration. 
 * This means a great reduction of the overhead due to latency of launching kernels and waiting for interrupts for each butterfly block, in consequence, overall performance is improved.
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
#include "cgra.h"
#include "cgra_bitstream.h"
#include "fxp.h"
#include "defines.h"
#include "fft_data.h"
#include "timer_sdk.h"

// Check CGRA size
#if CGRA_N_COLS != 4 || CGRA_N_ROWS != 4
  #error This example requires a 4x4 CGRA
#endif

// make clock-show
#ifndef CGRA_FREQ
#define CGRA_FREQ 100000000 // 100000000 Hz (100000 kHz)
#endif

#ifndef CPU_FREQ
#define CPU_FREQ 100000000 // 100000000 Hz (100000 kHz)
#endif

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

/* --------------------------------------------------------------------------
 *                     Global variables
 * --------------------------------------------------------------------------*/

// FFT radix-2 variables
fxp in_X_r[FFT_SIZE] __attribute__ ((aligned (4)));
fxp in_X_i[FFT_SIZE] __attribute__ ((aligned (4)));

/* --------------------------------------------------------------------------
 *                     Interrupt Handler
 * --------------------------------------------------------------------------*/

volatile int8_t cgra_intr_flag;

 void handler_irq_cgra(uint32_t id) {
    cgra_intr_flag = 1;
}

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
  printf("=== CGRA_O1 ===\n");

  uint32_t prog_cycles_start = 0;
  uint32_t prog_cycles_end = 0;
  uint32_t bitrev_cycles_start = 0;
  uint32_t bitrev_cycles_end = 0;
  uint32_t fft_cycles_start = 0;
  uint32_t fft_cycles_end = 0;
  CSR_READ(CSR_REG_MCYCLE, &prog_cycles_start);

  // Address Injection for SWI

  // Init Context Memory
  cgra_cmem_init(cgra_cmem_bitstream, cgra_kmem_bitstream);

  // Init Interrupts
  plic_Init();
  plic_irq_set_priority(CGRA_INTR, 1);
  plic_irq_set_enabled(CGRA_INTR, kPlicToggleEnabled);
  plic_assign_external_irq_handler(CGRA_INTR, (void*)&handler_irq_cgra);

  // Enable Global Interrupts
  CSR_SET_BITS(CSR_REG_MSTATUS, 0x8);
  CSR_SET_BITS(CSR_REG_MIE, 1 << 11);
  cgra_intr_flag = 0;

  cgra_t cgra;
  cgra.base_addr = mmio_region_from_addr((uintptr_t)CGRA_PERIPH_START_ADDRESS);

  cgra_wait_ready(&cgra);
  cgra_perf_cnt_enable(&cgra, 1);

  /* -------------------------------------------------------------------------------
  *
  * BIT REVERSE
  * 
  ------------------------------------------------------------------------------- */

  uint16_t stages = NumberOfBitsNeeded ( FFT_SIZE );

  // Load input data into the FFT input buffers (Real and Imaginary parts in separate arrays)
  for (int i = 0; i < FFT_SIZE; i++) {
    in_X_r[i] = input_signal[2*i + 0];
    in_X_i[i] = input_signal[2*i + 1];
  }

  CSR_READ(CSR_REG_MCYCLE, &bitrev_cycles_start);
  // Bit-reversal in-place for DIT FFT input ordering
  for (int i = 0; i < FFT_SIZE; i++) {
    uint16_t j = ReverseBits((uint16_t)i, stages);
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

  /* -------------------------------------------------------------------------------
  *
  * FFT RADIX-2 COMPUTATION
  * 
  ------------------------------------------------------------------------------- */

  printf("Run a complex FFT of %d points on CGRA...\n", FFT_SIZE);

  // AO wall-time timer: counts real elapsed time including WFI sleep.
  // End-to-end time aproximation
  timer_cycles_init();
  timer_start();

  cgra_perf_cnt_reset(&cgra);
  CSR_READ(CSR_REG_MCYCLE, &fft_cycles_start);

  // Buffers for input and output of the CGRA kernel, organized in columns as expected by the CGRA.
  static int32_t in_col0[FFT_SIZE] __attribute__((aligned(4)));
  static int32_t in_col1[FFT_SIZE/2 + 1] __attribute__((aligned(4))); //+1 because we need to store FFT_SIZE/2, the limit of the inner loop
  static int32_t in_col2[FFT_SIZE] __attribute__((aligned(4)));
  static int32_t in_col3[FFT_SIZE/2] __attribute__((aligned(4)));

  static int32_t out_col0[FFT_SIZE/2] __attribute__((aligned(4)));
  static int32_t out_col1[FFT_SIZE/2] __attribute__((aligned(4)));
  static int32_t out_col2[FFT_SIZE/2] __attribute__((aligned(4)));
  static int32_t out_col3[FFT_SIZE/2] __attribute__((aligned(4)));

  // Program stream base pointers
  cgra_wait_ready(&cgra);
  cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col0, 0);
  cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col1, 1);
  cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col2, 2);
  cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col3, 3);

  cgra_set_write_ptr(&cgra, (uint32_t)(uintptr_t)out_col0, 0);
  cgra_set_write_ptr(&cgra, (uint32_t)(uintptr_t)out_col1, 1);
  cgra_set_write_ptr(&cgra, (uint32_t)(uintptr_t)out_col2, 2);
  cgra_set_write_ptr(&cgra, (uint32_t)(uintptr_t)out_col3, 3);

  for (int i = 0; i < stages; i ++) { // For each FFT stage batch
    int L = 1 << i;
    int tw_step = FFT_SIZE / (2 * L);
    
    // Prepare pointers to input buffers for this stage
    int32_t *c0 = in_col0;
    int32_t *c1 = in_col1;
    int32_t *c2 = in_col2;
    int32_t *c3 = in_col3;
    *c1++ = (FFT_SIZE / 2); // Number of butterflies in this kernel launch
 
    // Prepare input data for stage
    for (int j = 0; j < FFT_SIZE; j += 2 * L){ // For each block of butterflies in this stage
      int tw_idx = 0;
      for(int k = 0; k < L; k++){ //For each butterfly in the block
        int index1 = j + k; // Index of the top element in the butterfly
        int index2 = index1 + L; // Index of the bottom element in the butterfly
        
        *c0++ = in_X_r[index2]; // bottom_r
        *c0++ = in_X_r[index1]; // top_r

        *c1++ = in_X_i[index2]; // bottom_i

        *c2++ = f_imag[tw_idx]; // wi
        *c2++ = in_X_i[index1]; // top_i

        *c3++ = f_real[tw_idx]; // wr

        tw_idx += tw_step;
      }
    }

    // Launch CGRA kernel for this stage
    cgra_intr_flag = 0;
    cgra_set_kernel(&cgra, CGRA_KERNEL);
    while (cgra_intr_flag == 0) {
      wait_for_interrupt();
    }

    // Prepare pointers to output buffers for this stage
    int32_t *o0 = out_col0;
    int32_t *o1 = out_col1;
    int32_t *o2 = out_col2;
    int32_t *o3 = out_col3;

    for (int j = 0; j < FFT_SIZE; j += 2 * L){ // For each block of butterflies in this stage
      for(int k = 0; k < L; k++){ // For each butterfly in the block
        int index1 = j + k; // Top index
        int index2 = index1 + L; // Bottom index

        in_X_r[index1] = (fxp)(*o1++); // top_r
        in_X_i[index1] = (fxp)(*o3++); // top_i

        in_X_r[index2] = (fxp)(*o0++); // bottom_r
        in_X_i[index2] = (fxp)(*o2++); // bottom_i
      }
    }
  }

  uint32_t fft_wall_cycles = timer_stop(); // Latency end-to-end including WFI sleep
  CSR_READ(CSR_REG_MCYCLE, &fft_cycles_end);
  printf("FFT done\n");

  /* -------------------------------------------------------------------------------
  *
  * CHECK RESULTS
  * 
  ------------------------------------------------------------------------------- */
  int32_t errors=0;
  printf("Results:\n");
  for (int i=0; i<FFT_SIZE; i++) {
    if(in_X_r[i] != exp_output_real[i] ||
        in_X_i[i] != exp_output_imag[i]) {
          printf("ERROR at index %d: ", i);
          errors++;
    } else  if (FFT_SIZE <= 8) { // Print all results for small FFT sizes
        printf("[%d] OUT Re=0x%08x Im=0x%08x | EXP Re=0x%08x Im=0x%08x\n",
        i,
        (uint32_t)in_X_r[i], (uint32_t)in_X_i[i],
        (uint32_t)exp_output_real[i], (uint32_t)exp_output_imag[i]);
    }
  }

  printf("CGRA FFT computation finished with %d errors\n", errors);
  CSR_READ(CSR_REG_MCYCLE, &prog_cycles_end);

  /* -------------------------------------------------------------------------------
  *
  * Performance Counters and Metrics
  * 
  ------------------------------------------------------------------------------- */
  // CGRA performance counters are per column. Active = configuration + execution.
  uint32_t cgra_cycles_active_max = 0;
  uint32_t cgra_cycles_stall_max = 0;
  uint32_t cgra_cycles_productive_max = 0;
  
  for (uint8_t col = 0; col < CGRA_N_COLS; col++) {
    uint32_t col_cycles = cgra_perf_cnt_get_col_active(&cgra, col);
    uint32_t col_stalls = cgra_perf_cnt_get_col_stall(&cgra, col);
    uint32_t col_productive = col_cycles - col_stalls;

    if (col_cycles > cgra_cycles_active_max) {
      cgra_cycles_active_max = col_cycles;
    }
    if (col_stalls > cgra_cycles_stall_max) {
      cgra_cycles_stall_max = col_stalls;
    }
    if (col_productive > cgra_cycles_productive_max) {
      cgra_cycles_productive_max = col_productive;
    }

    printf("CGRA column %u: active=%lu stall=%lu productive=%lu cycles\n",
      (unsigned)col,
      (unsigned long)col_cycles,
      (unsigned long)col_stalls,
      (unsigned long)col_productive);
  }

  unsigned long cgra_active_time_us = (unsigned long)(((uint64_t)cgra_cycles_active_max * 1000000ULL) / (uint64_t)CGRA_FREQ);
  unsigned long cgra_productive_time_us = (unsigned long)(((uint64_t)cgra_cycles_productive_max * 1000000ULL) / (uint64_t)CGRA_FREQ);

  printf("CGRA active time (configuration + execution): %lu CGRA cycles (%lu us)\n",
    (unsigned long)cgra_cycles_active_max,
    cgra_active_time_us);
  printf("CGRA stall cycles (max col): %lu CGRA cycles\n",
    (unsigned long)cgra_cycles_stall_max);
  printf("CGRA productive cycles (max per-column active-stall): %lu CGRA cycles (%lu us)\n",
      (unsigned long)cgra_cycles_productive_max,
      cgra_productive_time_us);

  unsigned long fft_wall_time_us = (unsigned long)(((uint64_t)fft_wall_cycles * 1000000ULL) / (uint64_t)CGRA_FREQ);

  printf("FFT offload wall-time (AO timer, includes WFI): %lu wall cycles (%lu us)\n",
    (unsigned long)fft_wall_cycles,
    fft_wall_time_us);
  printf("CGRA productive time (from perf counters): %lu us\n",
    cgra_productive_time_us);
  
  // Total CGRA kernels executed
  uint32_t cgra_kernels = cgra_perf_cnt_get_kernel(&cgra);

  // Total cycles spent by the CPU in the FFT offload window.
  // This is not the CGRA execution time; CGRA time is measured by perf counters above.
  uint32_t fft_cycles_total = fft_cycles_end - fft_cycles_start;

  // Total cycles for the bit-reversal section on the CPU.
  uint32_t bitrev_cycles_total = bitrev_cycles_end - bitrev_cycles_start;

  // Total cycles for the whole program.
  uint32_t prog_cycles_total = prog_cycles_end - prog_cycles_start;

  printf("CGRA kernels executed: %lu\n", (unsigned long)cgra_kernels);

  unsigned long bitrev_time_us = (unsigned long)(((uint64_t)bitrev_cycles_total * 1000000ULL) / (uint64_t)CPU_FREQ);
  unsigned long fft_cycles_time_cpu_us = (unsigned long)(((uint64_t)fft_cycles_total * 1000000ULL) / (uint64_t)CPU_FREQ);
  unsigned long prog_time_us = (unsigned long)(((uint64_t)prog_cycles_total * 1000000ULL) / (uint64_t)CPU_FREQ);

  long overhead_time_us = (long)fft_wall_time_us - (long)cgra_active_time_us;

  printf("FFT wall-time (AO timer): %lu us\n", fft_wall_time_us);
  printf("FFT section total cycles (mcycle, CPU only): %lu (%lu us)\n",
    (unsigned long)fft_cycles_total,
    fft_cycles_time_cpu_us);
  printf("Overhead time (wall - CGRA active): %ld us\n", overhead_time_us);
  printf("Bit-reversal total cycles (mcycle, CPU only): %lu (%lu us)\n",
    (unsigned long)bitrev_cycles_total,
    bitrev_time_us);
  printf("Program total cycles (mcycle, CPU only): %lu (%lu us)\n",
    (unsigned long)prog_cycles_total,
    prog_time_us);

  return errors ? EXIT_FAILURE : EXIT_SUCCESS;

}
 