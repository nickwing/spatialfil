#' Function for creating convolution kernel for Gaussian and Laplacian of Gaussian (LoG) filter
#' @description This function creates the convolution kernel for applying a filter to an array/matrix
#' @param sigma The \code{numeric} value of standard deviation for the Gaussian or LoG filter
#' @param k \code{character} value: \code{gaussian} for Gaussian kernel, \code{LoG} for Laplacian of Gaussian kernel, \code{sharpen} for
#' 3x3 convolution matrix for sharping the edges and \code{edge} for a 3x3 convolution matrix to trace the edges
#' @details The convolution kernel is a matrix that is used by \code{spacialfil} function over a matrix, or array, for filtering
#' the data. \emph{Gaussian}  kernel is calculated starting from the 2 dimension, isotropic, Gaussian distribution:
#' \deqn{G(x)=\frac{1}{2\pi\sigma^{2}}e^{-\frac{x^{2}+y^{2}}{2\sigma^{2}}}} \emph{Laplacian of Gaussian} kernel applies
#' a second derivative to enhance regions of rapid intensity changes:
#' \deqn{LoG\left ( x,y \right )=\frac{-1}{\pi\sigma^{4}}\left ( 1-\frac{x^{2}+y^{2}}{2\sigma^{2}}\right ) e^{-\frac{x^{2}+y^{2}}{2\sigma^{2}}}} the use of the underlying Gaussian kernel (so the name
#' Laplacian of Gaussian or \emph{LoG}) is needed to reduce the effect of high frequency noise that can affect the signal
#' distribution. \emph{Edge} kernel is a 3x3 cnvolution kernel used to enhance the edges of the matrix. \emph{Sharpen} enhance the detail
#' (but also the noise) in original dataset
#' @return A matrix with convolution kernel with size varying according the value of \code{sigma}
#' @export
#' @examples # creates a convolution kernel with Gaussian function and sigma = 1.4
#'  K <- convKernel(sigma = 1.4, k = 'gaussian')
convKernel <- function(sigma = 1.4, k = c('gaussian','LoG','sharpen','edge')) {
  k <- match.arg(k)
  l <- sigma * 7
  # check if odd and if not increas by one
  if (l%%2==0) l <- l + 1
  # dynamic adaptation of kernel size according the value of sigma
  x <- c(-floor(l/2):floor(l/2))
  y <- c(-floor(l/2):floor(l/2))
  if (k=='gaussian')
    M <- outer(X = x, Y = y, FUN = function(X, Y) return(1/(2*pi*sigma^2)*exp(-(X^2+Y^2)/(2*sigma^2))))

  if (k=='LoG')
    M <- outer(X = x, Y = y, FUN = function(X, Y) return(-1/(pi*sigma^4)*(1-(X^2+Y^2)/(2*sigma^2))*exp(-(X^2+Y^2)/(2*sigma^2))))

  if (k=='sharpen') M <- matrix(data = c(0,-1,0,-1,5,-1,0,-1,0), nrow = 3)
  if (k=='edge') M <- matrix(data = c(0,1,0,1,-4,1,0,1,0), nrow = 3)
  return(M)
}
