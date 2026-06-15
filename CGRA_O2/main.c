/*
 * CGRA_O2
 * FFT implementation (Second optimization)
 * Loop unrolling of tje FFT stages in the CGRA kernel, which allows to compute multiple stages in a single kernel launch, reducing the overhead of kernel launches and reconfiguration of read/write pointers.
 * The unrolling factor is defined by the CGRA_STAGES macro, which can be set in this file. The maximum value of CGRA_STAGES is the number of stages in the FFT, which is log2(FFT_SIZE).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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

#ifndef CGRA_STAGES
#define CGRA_STAGES 11 // Number of stages computed in each CGRA kernel launch
#endif

#if CGRA_STAGES < 1
  #error CGRA_STAGES must be >= 1
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
  printf("=== FFT CGRA (opti 2) ===\n");

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

  //SAFETY CHECK: CGRA_STAGES validation
  // Validate that CGRA_STAGES doesn't exceed the number of FFT stages
  if (CGRA_STAGES > stages) {
    printf("ERROR: CGRA_STAGES (%d) cannot be greater than FFT stages (%d)\n", CGRA_STAGES, stages);
    return EXIT_FAILURE;
  }

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
  CSR_READ(CSR_REG_MCYCLE, &fft_cycles_start); // CPU-side window only; CGRA time is measured by perf counters below.
  
   // Buffers for input and output of the CGRA kernel, organized in columns as expected by the CGRA.
  static int32_t in_col0[2 + CGRA_STAGES] __attribute__((aligned(4)));
  static int32_t in_col1[1] __attribute__((aligned(4)));
  static int32_t in_col2[3] __attribute__((aligned(4)));
  static int32_t in_col3[2 + CGRA_STAGES] __attribute__((aligned(4)));

  // Beginning FFT loop
  int stages_this_launch = CGRA_STAGES;
  for (int i = 0; i < stages; i += stages_this_launch) { // For each FFT stage 
    int ind0 = 0, ind1 = 0, ind2 = 0, ind3 = 0;

    in_col0[ind0++] = (int32_t)(uintptr_t)&in_X_r[0];
    in_col0[ind0++] = (int32_t)(uintptr_t)&in_X_r[0];

    in_col1[ind1++] = (int32_t)(uintptr_t)&in_X_i[0];

    in_col2[ind2++] = (int32_t)(uintptr_t)&f_imag[0];
    in_col2[ind2++] = (int32_t)(FFT_SIZE * (int)sizeof(fxp));
    in_col2[ind2++] = (int32_t)(uintptr_t)&in_X_i[0];

    in_col3[ind3++] = (int32_t)(uintptr_t)&f_real[0];
    if (i + stages_this_launch > stages) { // If remaining stages are less than CGRA_STAGES, adjust
      stages_this_launch = stages - i;
    }
    in_col3[ind3++] = (int32_t) stages_this_launch;

    for (int s = 0; s < stages_this_launch; s++) { // For each stage in this CGRA kernel launch
      int L = 1 << (i + s); // Current block size
      int tw_step = FFT_SIZE / (2 * L); // Current twiddle factor step size

      in_col0[ind0++] = (int32_t)(L * (int)sizeof(fxp));
      in_col3[ind3++] = (int32_t)(tw_step * (int)sizeof(fxp));
    }

    // Program stream base pointers
    cgra_wait_ready(&cgra);
    cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col0, 0);
    cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col1, 1);
    cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col2, 2);
    cgra_set_read_ptr(&cgra, (uint32_t)(uintptr_t)in_col3, 3);

    // Launch CGRA kernel for this stage
    cgra_intr_flag = 0;
    cgra_set_kernel(&cgra, CGRA_KERNEL);
    while (cgra_intr_flag == 0) {
      wait_for_interrupt();
    }
  }
  uint32_t fft_wall_cycles = timer_stop(); // Latency end-to-end including WFI sleep
  CSR_READ(CSR_REG_MCYCLE, &fft_cycles_end);
  printf("FFT done\n");
  printf("%d stages computed in each CGRA kernel launch\n", CGRA_STAGES);

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
  printf("CGRA stall cycles (max col): %lu CGRA cycles\n", (unsigned long)cgra_cycles_stall_max);
  printf("CGRA productive cycles (max per-column active-stall): %lu CGRA cycles (%lu us)\n",
      (unsigned long)cgra_cycles_productive_max,
      cgra_productive_time_us);

    unsigned long fft_wall_time_us = (unsigned long)(((uint64_t)fft_wall_cycles * 1000000ULL) / (uint64_t)CGRA_FREQ);

    printf("FFT offload wall-time (AO timer, includes WFI): %lu wall cycles (%lu us)\n",
      (unsigned long)fft_wall_cycles,
      fft_wall_time_us);
    printf("CGRA productive time (from perf counters): %lu us\n", cgra_productive_time_us);
  
  // Total CGRA kernels executed                         
  uint32_t cgra_kernels = cgra_perf_cnt_get_kernel(&cgra);

  // Total cycles spent by the CPU in the FFT offload window.
  // This is not the CGRA execution time; CGRA time is measured by perf counters above.
  uint32_t fft_cycles_total = fft_cycles_end - fft_cycles_start;

  // Total cycles for the bit-reversal section on the CPU.
  uint32_t bitrev_cycles_total = bitrev_cycles_end - bitrev_cycles_start;

  // Total cycles for the whole program.
  uint32_t prog_cycles_total = prog_cycles_end - prog_cycles_start;
  
  // Print summary of performance metrics
  printf("CGRA kernels executed: %lu\n", (unsigned long)cgra_kernels);
  printf("CGRA active cycles total (max col): %lu\n", (unsigned long)cgra_cycles_active_max);

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
