/* C code for convolution kernel application */
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

void applyKernel(double *dataInput, double *Kernel, int *extralines,
                 int *Nrow, int *Ncol, int *Nslices, double *dataOutput) {
  int ksize = *extralines * 2 + 1;            // size of convolution kernel matrix
  int m, n, p; // counters
  int index;   // vector index for calculating value
  for (m=0; m < *Nslices; m++) {              // first loop over slices
    for (n = *extralines; n < *Ncol - *extralines; n++) {   // loop among columns
      for (p = *extralines; p < *Nrow - *extralines; p++) { // loop among rows
        //printf("\nSlice=%d, Column=%d, Line=%d, pixel no=%d", m,n,p, m * (*Nrow * *Ncol) + n * *Ncol + p);
        index = m * (*Nrow * *Ncol) + n * *Ncol + p;

      }
    }
  }
  printf("\n");
}
