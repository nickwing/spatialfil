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
#' @return An object of class \code{convKern} with the \code{matrix} of convolution kernel whose size varies according the value of \code{sigma} (in case of
#' \code{gaussian} or \code{LoG} option selected), and \code(k) being the convolution kernel type label
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
  if (k=='laplacian') M <- matrix(data = c(0,1,0,1,-4,1,0,1,0), nrow = 3)
  # create S3 class
  output <- list('matrix' = M, 'kernel' = k)
  class(output) <- 'convKern'
  return(output)
}


#' Function for applying convolution kernel to a matrix or array
#' @description This function applies the a convolution kernel based filter to a \code{matrix} or \code{array} object type.
#' @param x An object of class \code{matrix} or \code{array}
#' @param kernel A \code{matrix} containing the values chosen as convolution kernel
#' @details The application of a convolution kernel over a 2D matrix dataset allows to apply functions as smoothing or edge detection.
#' The aim of this function is to filter 2D matrices in order to help signal finding across (images-derived) data. It is also possible to
#' filter 3D arrays considering them as slices of a series of images to be processed. Higher dimensions arrays are not allowed.
#' The \code{kernel} parameter is a simple \bold{square matrix} with an odd number of rows/columns, that can be pre-calculated by using
#' the function \code{\link{convKernel}}. Not square matrices or matrices with even number of rows/columns will exit an error.
#' @return An object with the same size of \code{x} containing data processed by convolution kernel
#' @import abind
#' @useDynLib spacialfil
#' @export
applyFilter <- function(x, kernel) {
  # check the k entry
  if (class(kernel)!='convKern') stop('kernel MUST be a convKern class object')
  # prepare matrix if x is matrix by adding required extra lines
  # number of extralines to be added to matrix and arrays (extraplanes)
  extralines <- dim(kernel$matrix)[1] %/% 2
  # now resizes array or matrix adding as new extralines (or planes) as the length
  # of half side of kernel - 1 (the center pixel)
  if (class(x)=='matrix') {
    for (n in 1:extralines) x <- cbind(x[,1], x, x[,ncol(x)])
    for (n in 1:extralines) x <- rbind(x[1,], x, x[nrow(x),])
  }
  if (class(x)=='array') {
    for (n in 1:extralines) x <- abind(x[,1,], x, x[,dim(x)[2],], along = 2)
    for (n in 1:extralines) x <- abind(x[1,,], x, x[dim(x)[1],,], along = 1)
  }
  Nrow <- dim(x)[1]
  Ncol <- dim(x)[2]
  if (class(x)=='matrix') Nslices <- 1 else Nslices <- dim(x)[3]
  dataOutput <- x
  # building the index matrix for transporting the kernel on the array or matrix
  kindex <- c()
  for (n in -extralines:extralines)  # index for columns
    for (m in -extralines:extralines)  # index for rows
      kindex <- c(kindex, n*Nrow + m)

  result <- .C('applyKernel', as.double(x), as.double(kernel$matrix), as.integer(extralines), as.integer(kindex),
               as.integer(Nrow), as.integer(Ncol), as.integer(Nslices), as.double(dataOutput))
  if (class(x)=='matrix') {
    output <- matrix(data = result[[8]], nrow = Nrow)
    # resize matrix to original size
    output <- output[(extralines+1):(nrow(output) - extralines), (extralines+1):(ncol(output) - extralines)]
  }
  if (class(x)=='array')  {
    output <- array(data = result[[8]], dim = c(Nrow, Ncol, Nslices))
    output <- output[(extralines+1):(nrow(output) - extralines), (extralines+1):(ncol(output) - extralines),]
  }
  return(output)
}
