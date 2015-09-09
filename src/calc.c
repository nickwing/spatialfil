/* C code for convolution kernel application */
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

/* kernel function for normalized results */
void applyKernel(double *dataInput, double *Kernel, int *extralines, int *kindex,
                 int *Nrow, int *Ncol, int *Nslices, double *dataOutput) {
  int ksize = (*extralines * 2 + 1) * (*extralines * 2 + 1); // size of convolution kernel matrix
  int m, n, p, q; // counters
  long index;     // vector index for calculating value
  double accData, accKernel; // accumulators for data and kernel
  accKernel = 0;
  for (m=0; m < ksize; m++) accKernel = accKernel + Kernel[m];  // calculate accumulated kernel
  for (m=0; m < *Nslices; m++) {                            // first loop over slices
    for (n = *extralines; n < *Ncol - *extralines; n++) {   // loop among columns
      for (p = *extralines; p < *Nrow - *extralines; p++) { // loop among rows
        //printf("\nSlice=%d, Column=%d, Line=%d, pixel no=%d", m,n,p, m * (*Nrow * *Ncol) + n * *Ncol + p);
        accData = 0;
        index = m * (*Nrow * *Ncol) + n * *Ncol + p;  // index of the voxel to be processed
        for (q=0; q < ksize; q++) {  // loop among kernel matrix
          accData = accData + dataInput[index + kindex[q]] * Kernel[q];
        }
        dataOutput[index] = accData / accKernel;
      }
    }
  }
}

/* kernel function without normalized results */
void applyKernelWithoutNorm(double *dataInput, double *Kernel, int *extralines, int *kindex,
                 int *Nrow, int *Ncol, int *Nslices, double *dataOutput) {
  int ksize = (*extralines * 2 + 1) * (*extralines * 2 + 1); // size of convolution kernel matrix
  int m, n, p, q; // counters
  long index;     // vector index for calculating value
  double accData; // accumulators for data and kernel
  for (m=0; m < *Nslices; m++) {                            // first loop over slices
    for (n = *extralines; n < *Ncol - *extralines; n++) {   // loop among columns
      for (p = *extralines; p < *Nrow - *extralines; p++) { // loop among rows
        //printf("\nSlice=%d, Column=%d, Line=%d, pixel no=%d", m,n,p, m * (*Nrow * *Ncol) + n * *Ncol + p);
        accData = 0;
        index = m * (*Nrow * *Ncol) + n * *Ncol + p;  // index of the voxel to be processed
        for (q=0; q < ksize; q++) {  // loop among kernel matrix
          accData = accData + dataInput[index + kindex[q]] * Kernel[q];
        }
        dataOutput[index] = accData;
      }
    }
  }
}