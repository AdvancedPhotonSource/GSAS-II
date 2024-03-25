/* Find "bad" regions on an image and create a pixel mask to remove them.
   This works by masking pixels that are m*sigma outside the range of the
   median at that radial distance. This is based on the AIRXD C++ code, 
   https://github.com/AdvancedPhotonSource/AIRXD-ML-PUB, developed by 
   Howard Yanxon, Wenqian Xu and James Weng.

    This is much faster than GSASIIimage.AutoPixelMask, which does almost 
    exactly the same computation, but uses pure Python/numpy code. 

    Version 2 (in progress, with speedup)
*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

// Report a Python TypeError Exception, message is reported as error
void * makeTypeException(const char *message) {
    PyErr_SetString(PyExc_TypeError, message);
    return NULL;
}
// Report a Python RuntimeError Exception, message is reported as error
void * makeRunException(const char *message) {
    PyErr_SetString(PyExc_RuntimeError, message);
    return NULL;
}

// Find median from an array
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }
typedef int elem_type;
elem_type quick_select_median(elem_type arr[], int n) {
  // Algorithm from Numerical Recipes in C (1992)
  // from https://stackoverflow.com/questions/1961173/median-function-in-c-math-library
    int low, high;
    int median;
    int middle, ll, hh;
    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
      if (high <= low) /* One element only */
	return arr[median] ;
      if (high == low + 1) { /* Two elements only */
	if (arr[low] > arr[high])
	  ELEM_SWAP(arr[low], arr[high]) ;
	return arr[median] ;
      }
      /* Find median of low, middle and high items; swap into position low */
      middle = (low + high) / 2;
      if (arr[middle] > arr[high])
	ELEM_SWAP(arr[middle], arr[high]) ;
      if (arr[low] > arr[high])
	ELEM_SWAP(arr[low], arr[high]) ;
      if (arr[middle] > arr[low])
	ELEM_SWAP(arr[middle], arr[low]) ;
      /* Swap low item (now in position middle) into position (low+1) */
      ELEM_SWAP(arr[middle], arr[low+1]) ;
      /* Nibble from each end towards middle, swapping items when stuck */
      ll = low + 1;
      hh = high;
      for (;;) {
	do ll++; while (arr[low] > arr[ll]) ;
	do hh--; while (arr[hh] > arr[low]) ;
	if (hh < ll)
	  break;
	ELEM_SWAP(arr[ll], arr[hh]) ;
      }
      /* Swap middle item (in position low) back into correct position */
      ELEM_SWAP(arr[low], arr[hh]) ;
      /* Re-set active partition */
      if (hh <= median)
	low = ll;
      if (hh >= median)
	high = hh - 1;
    }
    return arr[median] ;
}

static PyObject *fmask_func(PyObject *self, PyObject *args) {
  /* Function called from Python for masking bad pixels, e.g. where the
     pixel intensity is more than n*sigma outside the statistically expected
     range. Comparisons are done with the median value against sigma from 
     the Median absolute deviation, where both quantities are not shifted 
     significantly by a small number of outliers. The value for "n" is 
     provided from parameter esdmul.

     This is an accelerated version of the prevmask_func() below.

  Python parameters: 
    :param float esdmul: multiplier to locate bad pixels, typically 3. 
      A smaller number locates more "bad" pixels.
    :param np.array frame: a 1D data structure of type Bool, where True 
      indicates a pixel that is maxed by the pixel mask and can be ignored.
    :param np.array TA: a 1D data structure of type float with the 2theta 
      value for every pixel.
    :param np.array image: a 1D data structure of type np.int32 with the 
      intensity for each pixel.
    :param np.array TThs: a 1D data structure of type float with the 
      2theta values for each integration step. The pixels within a step
      are expected to have the same intensity.
    :param np.array outMask: a 1D data structure of type bool that is used
      for the output from the function. Pixels to be masked are set to True.
    :param float ttmin: minimum 2theta value for bad pixel scan. Pixels with 
      2theta significantly below this value are not examined. Note that if any 
      pixels in a ring are inside the ttmin to ttmax range (based on TThs), 
      all pixels in that ring are used. 
    :param float ttmax: maximum 2theta value for bad pixel scan. Pixels with 
      2theta significantly above this value are not examined. See ttmin
      for discussion of range expansion. 

    :returns: a count with the number of pixels that are set to True in 
      the current bad pixel scan. 
  */
  double esdmul; // arg 0
  PyArrayObject *frame, *TA, *image, *TThs, *outMask; // args 1-5
  double ttmin, ttmax; // arg 6 & 7
  
  double* TThsP;  // pointer to 2Theta bin array
  char *frameP, *outMaskP; // pointers to args
  double *TAP; // pointer to TA arg
  int *imageP; // pointers to image arg
  int masked = 0; // total masked pixel counter

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "dOOOOOdd", &esdmul,
			&frame, &TA, &image, &TThs, &outMask,
			&ttmin, &ttmax)) {
        return NULL;
  }
 
  // printf("esdmul = %lf 2theta min=%lf, max=%lf\n",esdmul,ttmin,ttmax);

  /* Sanity checks on function arguments */
  if (PyArray_DESCR(frame)->type_num != NPY_BOOL ||
      PyArray_STRIDE(frame, 0) != PyArray_ITEMSIZE(frame) ||
      PyArray_NDIM(frame) != 1)
    return makeTypeException("Arg 1 (frame): expected 1-dimensional array of type bool w/stride 1.");
  if (PyArray_DESCR(TA)->type_num != NPY_DOUBLE ||
      PyArray_STRIDE(TA, 0) != PyArray_ITEMSIZE(TA) ||
      PyArray_NDIM(TA) != 1)
    return makeTypeException("Arg 2 (TA): expected 1-dimensional array of type np.float64 w/stride 1.");
  if (PyArray_DESCR(image)->type_num != NPY_INT32 ||
      PyArray_STRIDE(image, 0) != PyArray_ITEMSIZE(image) ||
      PyArray_NDIM(image) != 1)
    return makeTypeException("Arg 3 (image): expected 1-dimensional array of type np.float64 w/stride 1.");
  if (PyArray_DESCR(TThs)->type_num != NPY_DOUBLE ||
      PyArray_STRIDE(TThs, 0) != PyArray_ITEMSIZE(TThs) ||
      PyArray_NDIM(TThs) != 1)
    return makeTypeException("Arg 4 (TThs): expected 1-dimensional array of type np.float64 w/stride 1.");
  if (PyArray_DESCR(outMask)->type_num != NPY_BOOL ||
      PyArray_STRIDE(outMask, 0) != PyArray_ITEMSIZE(outMask) ||
      PyArray_NDIM(outMask) != 1)
    return makeTypeException("Arg 3 (outMask): expected 1-dimensional array of type bool w/stride 1.");
  // special tests for the output array
  if (!(PyArray_FLAGS(outMask) & NPY_ARRAY_C_CONTIGUOUS) ||
      !(PyArray_FLAGS(outMask) & NPY_ARRAY_ALIGNED) ||
      !(PyArray_FLAGS(outMask) & NPY_ARRAY_WRITEABLE)) {
    return makeTypeException("Arg 3 (outMask): array is not continuous, aligned or writeable");
  }

  /* verify array sizes match */
  int size = PyArray_DIM(image, 0); // size of image, etc. arrays
  if (size != PyArray_DIM(frame, 0))
    return makeTypeException("Arg 1 (frame): array size does not match image.");	
  if (size != PyArray_DIM(TA, 0))
    return makeTypeException("Arg 2 (TA): array size does not match image.");	
  if (size != PyArray_DIM(outMask, 0))
    return makeTypeException("Arg 5 (outMask): array size does not match image.");

  /* prepare to access Python arrays */
  TThsP = (double *) PyArray_DATA(TThs);
  frameP = (char *) PyArray_DATA(frame);
  outMaskP = (char *) PyArray_DATA(outMask);
  TAP = (double *) PyArray_DATA(TA);
  imageP = (int *) PyArray_DATA(image);
  // 2theta ring info 
  int nth = PyArray_DIM(TThs, 0);     // number of 2Theta bin values
  double deltaTT = (TThsP[nth-1] - TThsP[0]) / (nth-1); // 2theta step 
  // printf("TThs [0]=%f, [-1]=%f, delta=%f\n",TThsP[0],TThsP[nth-1],deltaTT);

  // array to count number of pixels in each 2theta ring
  int *size2Th = (int*) malloc(nth * sizeof(int));
  //int size2Th[nth];
  for (int ith=0; ith < nth; ++ith) {
    size2Th[ith] = 0;
  }
  int i2TTh;
  double tth;
  // count the number of pixels in each 2theta ring; no tests to max out speed
  for (int i=0; i < size; ++i) {     // loop over all pixels
    if (TAP[i] < TThsP[0]) continue;
    i2TTh = (int) ((TAP[i] - TThsP[0]) / deltaTT);
    if (i2TTh < 0 || i2TTh >= nth) continue;
    size2Th[i2TTh] += 1;
  }
  // make array of pointers to lists of intensities & positions for each ring
  //int *pixIntP[nth];
  //int *pixLocP[nth];
  int **pixIntP = (int**) malloc(nth * sizeof(int*));
  int **pixLocP = (int**) malloc(nth * sizeof(int*));
  for (int ith=0; ith < nth; ++ith) {
    pixIntP[ith] = (int*) malloc(size2Th[ith] * sizeof(int));
    pixLocP[ith] = (int*) malloc(size2Th[ith] * sizeof(int));
  }

  for (int ith=0; ith < nth; ++ith) {
    size2Th[ith] = 0;
  }
  // loop over pixels and this time put them into the lists; test w/frame mask
  // but include pixels that are in outside integration range in statistics
  for (int i=0; i < size; ++i) {     // loop over all pixels
    if (frameP[i] == 1) continue;    // inside frame mask?
    if (TAP[i] < TThsP[0]) continue;
    i2TTh = (int) ((TAP[i] - TThsP[0]) / deltaTT);
    if (i2TTh < 0 || i2TTh >= nth) continue;
    // put in lists
    pixIntP[i2TTh][size2Th[i2TTh]] = imageP[i];
    pixLocP[i2TTh][size2Th[i2TTh]] = i;
    size2Th[i2TTh] += 1;
  }
  // loop over each ring where 2theta is ~same and where some pixels are in
  // allowed range; compute the median and MAD for the current ring and
  // then remove pixels that are within limits
  for (int ith=0; ith < nth; ++ith) {  // loop over rings 
    if (TThsP[ith]+deltaTT < ttmin || TThsP[ith] > ttmax) continue; // ring not needed
    int count = size2Th[ith];
    if (count < 10) continue;   // too few pixels to do any meaningful exclusion

    int *pixelsInt = pixIntP[ith];
    int *pixelsLoc = pixLocP[ith];

    int median = quick_select_median(pixelsInt, count); // Get median value for pixels in ring

    // Get the Median absolute deviation (MAD, https://en.wikipedia.org/wiki/Median_absolute_deviation)
    int *medianDev = (int*) malloc(count * sizeof(int));
    for (int i=0; i < count; ++i) {
      medianDev[i] = abs(pixelsInt[i] - median);
    }
    double MAD = 1.4826 * quick_select_median(medianDev, count);
    free(medianDev);
    /* 1.4826 scales to normal dist, to match use of median_absolute_deviation
       in pure-Python version of masking, also see 
       https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html
    */
    // printf("%d %f %f %d %f\n",ith,TThsP[ith],TThsP[ith]+deltaTT,median,MAD);
    
    // now look for pixels that are above and below the threshold, but
    // use strict 2theta limits
    double mintt = fmax(ttmin,TThsP[0]);
    double maxtt = fmin(ttmax,TThsP[nth-1]+deltaTT);
    for (int i=0; i < count; ++i) { // loop over pixels in ring
      int j = pixelsLoc[i];
      if (TAP[j] < mintt || TAP[j] > maxtt) continue;
      if ((imageP[j] > median + esdmul*MAD) ||
	  (imageP[j] < median - esdmul*MAD)) {
	  outMaskP[j] = 1;
	  masked += 1;
      }
    } // end loop over pixels in ring
    //printf("%d %f %d %d %f %d\n",ith,TThsP[ith],count,median,MAD,imask);
    free(pixIntP[ith]);
    free(pixLocP[ith]);
  } // end loop over rings
  free(size2Th);
  free(pixIntP);
  free(pixLocP);
  return PyLong_FromLong(masked);
}

static PyObject *prevmask_func(PyObject *self, PyObject *args) {
  /* Function called from Python for masking bad pixels, e.g. where the
     pixel intensity is more than n*sigma outside the statistically expected
     range. Comparisons are done with the median value against sigma from 
     the Median absolute deviation, where both quantities are not shifted 
     significantly by a small number of outliers. The value for "n" is 
     provided from parameter esdmul.

  Python parameters: 
    :param float esdmul: multiplier to locate bad pixels, typically 3. 
      A smaller number locates more "bad" pixels.
    :param np.array frame: a 1D data structure of type Bool, where True 
      indicates a pixel that is maxed by the pixel mask and can be ignored.
    :param np.array TA: a 1D data structure of type float with the 2theta 
      value for every pixel.
    :param np.array image: a 1D data structure of type np.int32 with the 
      intensity for each pixel.
    :param np.array TThs: a 1D data structure of type float with the 
      2theta values for each integration step. The pixels within a step
      are expected to have the same intensity.
    :param np.array outMask: a 1D data structure of type bool that is used
      for the output from the function. Pixels to be masked are set to True.
    :param float ttmin: minimum 2theta value for bad pixel scan. Pixels with 
      2theta significantly below this value are not examined. Note that if any 
      pixels in a ring are inside the ttmin to ttmax range (based on TThs), 
      all pixels in that ring are used. 
    :param float ttmax: maximum 2theta value for bad pixel scan. Pixels with 
      2theta significantly above this value are not examined. See ttmin
      for discussion of range expansion. 

    :returns: a count with the number of pixels that are set to True in 
      the current bad pixel scan. 
  */
  double esdmul; // arg 0
  PyArrayObject *frame, *TA, *image, *TThs, *outMask; // args 1-5
  double ttmin, ttmax; // arg 6 & 7
  
  double* TThsP;  // pointer to 2Theta bin array
  char *frameP, *outMaskP; // pointers to args
  double *TAP; // pointer to TA arg
  int *imageP; // pointers to image arg
  int masked = 0; // total masked pixel counter

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "dOOOOOdd", &esdmul,
			&frame, &TA, &image, &TThs, &outMask,
			&ttmin, &ttmax)) {
        return NULL;
  }
 
  // printf("esdmul = %lf 2theta min=%lf, max=%lf\n",esdmul,ttmin,ttmax);

  /* Sanity checks on function arguments */
  if (PyArray_DESCR(frame)->type_num != NPY_BOOL ||
      PyArray_STRIDE(frame, 0) != PyArray_ITEMSIZE(frame) ||
      PyArray_NDIM(frame) != 1)
    return makeTypeException("Arg 1 (frame): expected 1-dimensional array of type bool w/stride 1.");
  if (PyArray_DESCR(TA)->type_num != NPY_DOUBLE ||
      PyArray_STRIDE(TA, 0) != PyArray_ITEMSIZE(TA) ||
      PyArray_NDIM(TA) != 1)
    return makeTypeException("Arg 2 (TA): expected 1-dimensional array of type np.float64 w/stride 1.");
  if (PyArray_DESCR(image)->type_num != NPY_INT32 ||
      PyArray_STRIDE(image, 0) != PyArray_ITEMSIZE(image) ||
      PyArray_NDIM(image) != 1)
    return makeTypeException("Arg 3 (image): expected 1-dimensional array of type np.float64 w/stride 1.");
  if (PyArray_DESCR(TThs)->type_num != NPY_DOUBLE ||
      PyArray_STRIDE(TThs, 0) != PyArray_ITEMSIZE(TThs) ||
      PyArray_NDIM(TThs) != 1)
    return makeTypeException("Arg 4 (TThs): expected 1-dimensional array of type np.float64 w/stride 1.");
  if (PyArray_DESCR(outMask)->type_num != NPY_BOOL ||
      PyArray_STRIDE(outMask, 0) != PyArray_ITEMSIZE(outMask) ||
      PyArray_NDIM(outMask) != 1)
    return makeTypeException("Arg 3 (outMask): expected 1-dimensional array of type bool w/stride 1.");
  // special tests for the output array
  if (!(PyArray_FLAGS(outMask) & NPY_ARRAY_C_CONTIGUOUS) ||
      !(PyArray_FLAGS(outMask) & NPY_ARRAY_ALIGNED) ||
      !(PyArray_FLAGS(outMask) & NPY_ARRAY_WRITEABLE)) {
    return makeTypeException("Arg 3 (outMask): array is not continuous, aligned or writeable");
  }

  /* verify array sizes match */
  int size = PyArray_DIM(image, 0); // size of image, etc. arrays
  if (size != PyArray_DIM(frame, 0))
    return makeTypeException("Arg 1 (frame): array size does not match image.");	
  if (size != PyArray_DIM(TA, 0))
    return makeTypeException("Arg 2 (TA): array size does not match image.");	
  if (size != PyArray_DIM(outMask, 0))
    return makeTypeException("Arg 5 (outMask): array size does not match image.");

  /* prepare to access Python arrays */
  TThsP = (double *) PyArray_DATA(TThs);
  int nth = PyArray_DIM(TThs, 0);     // number of 2Theta bin values
  double deltaTT = (TThsP[nth-1] - TThsP[0]) / (nth-1); // 2theta step 
  // printf("TThs [0]=%f, [-1]=%f, delta=%f\n",TThsP[0],TThsP[nth-1],deltaTT);
  frameP = (char *) PyArray_DATA(frame);
  outMaskP = (char *) PyArray_DATA(outMask);
  TAP = (double *) PyArray_DATA(TA);
  imageP = (int *) PyArray_DATA(image);

  // create storage for a list of pixel intensities & locations
  // will be set to contents of current 2theta ring
  int *pixelsInt = (int*) malloc(size * sizeof(int));
  int *pixelsLoc = (int*) malloc(size * sizeof(int));
  
  for (int ith=0; ith < nth; ++ith) {  // loop over rings where 2theta is ~same
    if (TThsP[ith]+deltaTT < ttmin || TThsP[ith] > ttmax) continue; // & within 2theta limits
    int count = 0;
    for (int i=0; i < size; ++i) {     // loop over all pixels
      if (frameP[i] == 1) continue;    // inside frame mask?
      if (TAP[i] >= TThsP[ith] && TAP[i] < TThsP[ith]+deltaTT) {
	// Make a list of intensities & their array locations
	// for pixels that are in the current ring
	pixelsInt[count] = imageP[i];
	pixelsLoc[count] = i;
	count += 1;
      }
    }  // end loop over all pixels
    if (count < 10) continue;   // too few pixels to do any meaningful exclusion
    
    // Get median value for selected pixels
    /* int cmpfunc (const void * a, const void * b) { // used in qsort
           return ( *(int*)a - *(int*)b );
	   }
       qsort(pixelsInt, count, sizeof(int), cmpfunc);
       int median = pixelsInt[count/2]; */
    // replacing sort code above with NR call below improves speed by a few %(?)
    int median = quick_select_median(pixelsInt, count);

    // Get the Median absolute deviation (MAD, https://en.wikipedia.org/wiki/Median_absolute_deviation)
    int *medianDev = (int*) malloc(count * sizeof(int));
    for (int i=0; i < count; ++i) {
      medianDev[i] = abs(pixelsInt[i] - median);
    }
    double MAD = 1.4826 * quick_select_median(medianDev, count);
    free(medianDev);
    /* 1.4826 scales to normal dist, to match use of median_absolute_deviation
       in pure-Python version of masking, also see 
       https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.median_abs_deviation.html
    */
    
    // now look for pixels that are above and below the threshold
    for (int i=0; i < count; ++i) {
      int j = pixelsLoc[i];
      if ((imageP[j] > median + esdmul*MAD) ||
	  (imageP[j] < median - esdmul*MAD)) {
	  outMaskP[j] = 1;
	  masked += 1;
      }
    } // end loop over selected pixels
    //printf("%d %f %f %d %d %f %d\n",ith,TThsP[ith],TThsP[ith]+deltaTT,count,median,MAD,imask);
  } // loop over rings
  free(pixelsInt);
  free(pixelsLoc);
  return PyLong_FromLong(masked);
}

// Standard Python package integration code follows ===========================

static PyMethodDef ModuleContents[] = { /* define command(s) supplied here */
    {"prevmask", prevmask_func, METH_VARARGS, "fast masking routine"},
    {"mask", fmask_func, METH_VARARGS, "faster masking routine"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fmaskmodule = { /* define the module */
    PyModuleDef_HEAD_INIT,
    "fmask",
    "Fast Image Masking Module for GSAS-II",
    -1,
    ModuleContents
};

PyMODINIT_FUNC PyInit_fmask(void) { /* initialize the module */
    return PyModule_Create(&fmaskmodule);
    import_array();
}
