## @file
#  @brief The file containing the core definitions for the XRD Fundamental Parameneters Model (FPA) computation in python
#  @mainpage
#  @author Marcus H. Mendenhall (marcus.mendenhall@nist.gov)
#  @date March, 2015, updated August 2018
#  @copyright
#  The "Fundamental Parameters Python Code" ("software") is provided by the National Institute of Standards and Technology (NIST),
#  an agency of the United States Department of Commerce, as a public service.
#  This software is to be used for non-commercial research purposes only and is expressly provided "AS IS."
#  Use of this software is subject to your acceptance of these terms and conditions.
#  \n \n
#  NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION OF LAW,
#  INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  NON-INFRINGEMENT AND DATA ACCURACY.
#  NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED
#  OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED.
#  NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,
#  INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.
#  \n \n
#  NIST SHALL NOT BE LIABLE AND YOU HEREBY RELEASE NIST FROM LIABILITY FOR ANY INDIRECT, CONSEQUENTIAL,
#  SPECIAL, OR INCIDENTAL DAMAGES (INCLUDING DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION,
#  LOSS OF BUSINESS INFORMATION, AND THE LIKE), WHETHER ARISING IN TORT, CONTRACT, OR OTHERWISE,
#  ARISING FROM OR RELATING TO THE SOFTWARE (OR THE USE OF OR INABILITY TO USE THIS SOFTWARE),
#  EVEN IF NIST HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
#  \n \n
#  Software developed by NIST employees is not subject to copyright protection within the United States.
#  By using this software, creating derivative works or by incorporating this software into another product,
#  you agree that you will use the software only for non-commercial research purposes
#  and will indemnify and hold harmless the United States Government
#  for any and all damages or liabilities that arise out of any use by you.
#
#
#  @file
#  This file contains the main computational class FP_profile, which stores cached information to allow it
#  to efficiently recompute profiles when parameters have been modified.  It also includes some
#  helper classes, and a function @ref fourier_line_profile(), which has too many parameters,
#  which computes a profile in a single call.
#
#  @section intro_sec Introduction
#  Compute diffractometer line profile functions by methods
#  from Cheary & Coelho 1998 and Mullen & Cline paper and 'R' package.
#  Accumulate all convolutions in Fourier space, for efficiency,
#  except for axial divergence, which needs to be weighted in real space for I3 integral.
#
# changelog:
#
# 15 August, 2018 -- MHM
# fixed apparent error in flip of eps0 across twotheta=90 around line 549
#
# @file

from __future__ import print_function

import os, sys, math
import numpy
## @brief figure out which FFT package we have, and import it
try:
    from pyfftw.interfaces import numpy_fft, cache
    ## recorded variant of real fft that we will use
    best_rfft=numpy_fft.rfft
    ## recorded variant of inverse real fft that we will use
    best_irfft=numpy_fft.irfft
    ## @brief create a table of nice factorizations for the FFT package
    #
    #  this is built once, and shared by all instances
    #  fftw can handle a variety of transform factorizations
    ft_factors=[
        2*2**i*3**j*5**k for i in xrange(20) for j in range(10) for k in range(8)
        if 2*2**i*3**j*5**k <= 1000000
    ]
    cache.enable()
    cache.set_keepalive_time(1.0)
except ImportError:
    best_rfft=numpy.fft.rfft
    best_irfft=numpy.fft.irfft
    ## @brief create a table of nice factorizations for the FFT package
    #
    # this is built once, and shared by all instances
    # numpy fft is not too efficient for things other than a power of two,
    # although my own measurements says it really does fine.  For now, leave all factors available
    ft_factors=[
        2*2**i*3**j*5**k for i in range(20) for j in range(10) for k in range(8)
        if 2*2**i*3**j*5**k <= 1000000
    ]

ft_factors.sort()
ft_factors=numpy.array(ft_factors, numpy.int)

## @brief used for debugging moments from FP_profile.axial_helper().
moment_list=[]
## @brief if this is \a True, compute and save moment errors
collect_moment_errors=False

## @brief a skeleton class which makes a combined dict and namespace interface for easy pickling and data passing
class profile_data(object):
    ## @brief initialize the class
    #  @param self self
    #  @param kwargs keyword=value list to pre-populate the class
    def __init__(self, **kwargs):
        mydict={}
        mydict.update(kwargs)
        for k,v in mydict.items():
            setattr(self,k,v)
        ## a dictionary which shadows our attributes.
        self.dictionary=mydict

    ## @brief add new symbols to both the attributes and dictionary for the class
    #  @param self self
    #  @param kwargs keyword=value pairs
    def add_symbol(self, **kwargs):
        self.dictionary.update(kwargs)
        for k,v in kwargs.items():
            setattr(self,k,v)

## @brief the main class, which handles a single reflection
#
#  This class is designed to be highly extensible by inheriting new convolvers.  When it is initialized,
#  it scans its namespace for specially formatted names, which can come from mixin classes.
#  If it finds a function name of the form conv_xxx, it will call this funtion to create a convolver.
#  If it finds a name of the form info_xxx it will associate the dictionary with that convolver, which
#  can be used in UI generation, for example.  The class, as it stands, does nothing significant with it.
#  If it finds str_xxx, it will use that function to format a printout of the current state of the
#  convolver conv_xxx, to allow improved report generation for convolvers.
#
#  When it is asked to generate a profile, it calls all known convolvers.  Each convolver returns the
#  Fourier transform of its convolvution.  The transforms are multiplied together, inverse transformed,
#  and after fixing the periodicity issue, subsampled, smoothed and returned.
#
#  If a convolver returns \a None, it is not multipled into the product.
class FP_profile:

    ## @brief the number of histories to cache; can be overridden if memory is an issue.
    max_history_length=5 #max number of histories for each convolver
    
    ## @brief length_scale_m sets scaling for nice printing of parameters.
    #   if the units are in mm everywhere, set it to 0.001, e.g.
    #   convolvers which implement their own str_xxx method may use this to format their results,
    #   especially if 'natural' units are not meters.  Typical is wavelengths and lattices in nm or angstroms,
    #   for example.
    length_scale_m = 1.0
    
    ## @brief initialize the instance
    #
    #  @param anglemode = 'd' if setup will be in terms of a d-spacing,
    #  otherwise 'twotheta' if setup will be at a fixed 2theta value.
    #  @param output_gaussian_smoother_bins_sigma the number of bins for post-smoothing of data. 1.0 is good.
    #  @param oversampling the number of bins internally which will get computed for each bin the the final result.
    def __init__(self, anglemode,
            output_gaussian_smoother_bins_sigma=1.0,
            oversampling=10
        ):
        if anglemode not in ("d","twotheta"):
            raise Exception("invalid angle mode %s, must be 'd' or 'twotheta'"%anglemode)
        ## set to either 'd' for d-spacing based position, or 'twotheta' for angle-based position
        self.anglemode=anglemode
        ## sigma, in units of bins, for the final smoother. \a None means no final smoothing step.
        self.output_gaussian_smoother_bins_sigma=output_gaussian_smoother_bins_sigma
        ## the number of internal bins computed for each bin in the final output. 5-10 is usually plenty.
        self.oversampling=oversampling
        ## List of our convolvers, found by introspection of names beginning with 'conv_'
        self.convolvers=convolvers=[x for x in dir(self) if x.startswith("conv_")]
        ## A dictionary which will store all the parameters local to each convolution
        self.param_dicts=dict([(c,{}) for c in convolvers])
        #add global parameters, associated with no specific convolver
        ## A dictionary of bound functions to call to compute convolutions
        self.convolver_funcs=dict([(x, getattr(self,x)) for x in convolvers])
        ## If \a True, print cache hit information
        self.debug_cache=False

    ## @brief return the name of the function that called this.
    #  @param self self
    #  @return name of calling function
    #
    #  Useful for convolvers to identify themselves
    def get_function_name(self):
        return sys._getframe(1).f_code.co_name

    ## @brief add a numpy array to the list of objects that can be thrown away on pickling.
    #  @param self self
    #  @param b the buffer to add to the list
    #  @return the same buffer, to make nesting easy.
    def add_buffer(self, b):
        self._clean_on_pickle.add(id(b))
        return b
    
    ## @brief move the computation window to a new position, and clear history.
    #  @param self self
    #  @param twotheta_window_center_deg the center position of the middle bin of the window, in degrees
    #  @param twotheta_window_fullwidth_deg the full width of the window, in degrees
    #  @param twotheta_output_points the number of bins in the final output
    #
    #  sets up many arrays
    def set_window(self, twotheta_window_center_deg, twotheta_window_fullwidth_deg,
            twotheta_output_points
        ):
        """move the compute window to a new location and compute grids, 
            without resetting all parameters.  Clears convolution history.
        """
        
        ## the saved width of the window, in degrees
        self.twotheta_window_fullwidth_deg=twotheta_window_fullwidth_deg
        ## the saved center of the window, in degrees
        self.twotheta_window_center_deg=twotheta_window_center_deg
        ## the width of the window, in radians
        self.window_fullwidth=window_fullwidth=twotheta_window_fullwidth_deg*math.pi/180
        ## the center of the window, in radians
        self.twotheta_window_center=twotheta=twotheta_window_center_deg*math.pi/180.0
        ## the number of points to compute in the final results
        self.twotheta_output_points=twotheta_output_points
        ## the number of points in Fourier space to compute
        self.n_omega_points=nn=self.oversampling*twotheta_output_points//2+1
        #build all the arrays
        x=self._clean_on_pickle=set() #keep a record of things we don't keep when pickled
        b=self.add_buffer #shortcut
        
        ## a real-format scratch buffer
        self._rb1=b(numpy.zeros(nn,numpy.float))
        ## a real-format scratch buffer
        self._rb2=b(numpy.zeros(nn,numpy.float))
        ## a real-format scratch buffer
        self._rb3=b(numpy.zeros(nn,numpy.float))
        ## a complex-format scratch buffer
        self._cb1=b(numpy.zeros(nn,numpy.complex))
        ## a scratch buffer used by the axial helper
        self._f0buf=b(numpy.zeros(self.oversampling*twotheta_output_points, numpy.float))
        ## a scratch buffer used for axial divergence
        self._epsb2=b(numpy.zeros(self.oversampling*twotheta_output_points, numpy.float))
        ## the I2+ buffer
        self._I2p=b(numpy.zeros(self.oversampling*twotheta_output_points, numpy.float))
        ## the I2- buffer
        self._I2m=b(numpy.zeros(self.oversampling*twotheta_output_points, numpy.float))
        ## another buffer used for axial divergence
        self._axial=b(numpy.zeros(self.oversampling*twotheta_output_points, numpy.float))
        ## the largest frequency in Fourier space
        omega_max=self.n_omega_points*2*math.pi/window_fullwidth
        #build the x grid and the complex array that is the convolver
        #omega is in inverse radians in twotheta space (!)
        #i.e. if a transform is written
        #I(ds) = integral(A(L) exp(2 pi i L ds) dL where L is a real-space length,
        #and s=2 sin(twotheta/2)/lambda
        #then ds=2*pi*omega*cos(twotheta/2)/lambda (double check this!)
        ## The grid in Fourier space, in inverse radians
        self.omega_vals=b(numpy.linspace(0, omega_max, self.n_omega_points, endpoint=True))
        ## The grid in Fourier space, in inverse degrees
        self.omega_inv_deg=b(self.omega_vals*(math.pi/180))
        ## The grid in real space, in radians, with full oversampling
        self.twothetasamples=b(numpy.linspace(
            twotheta-window_fullwidth/2.0,twotheta+window_fullwidth/2.0,
            self.twotheta_output_points*self.oversampling, endpoint=False))
        ## The grid in real space, in degrees, with full oversampling
        self.twothetasamples_deg=b(numpy.linspace(
            twotheta_window_center_deg-twotheta_window_fullwidth_deg/2.0,
            twotheta_window_center_deg+twotheta_window_fullwidth_deg/2.0,
            self.twotheta_output_points*self.oversampling, endpoint=False))
        ## Offsets around the center of the window, in radians
        self.epsilon=b(self.twothetasamples-twotheta)
        
        ## A dictionary in which we collect recent state for each convolution.
        ## whenever the window gets reset, all of these get cleared
        self.convolution_history=dict([(x, []) for x in self.convolvers])

        ## A dictionary of Lorentz widths, used for de-periodizing the final result.
        self.lor_widths={}

    ## @brief find a bin count close to what we need, which works well for Fourier transforms.
    #  @param self self
    #  @param count a number of bins.
    #  @return a bin count somewhat larger than \a count which is efficient for FFT
    def get_good_bin_count(self, count):
        return ft_factors[ft_factors.searchsorted(count)]

    ## @brief given a window center, and approximate width, and an exact bin spacing,
    #  set up with a really good, slightly wider window for efficient FFT.
    #  @param self self
    #  @param twotheta_window_center_deg exact position of center bin, in degrees
    #  @param twotheta_approx_window_fullwidth_deg approximate desired width
    #  @param twotheta_exact_bin_spacing_deg the exact bin spacing to use
    def set_optimized_window(self, twotheta_window_center_deg,
        twotheta_approx_window_fullwidth_deg,
        twotheta_exact_bin_spacing_deg):
        """pick a bin count which factors cleanly for FFT, and adjust the window width
        to preserve the exact center and bin spacing"""
        bins=self.get_good_bin_count(int(1+
            twotheta_approx_window_fullwidth_deg/twotheta_exact_bin_spacing_deg))
        window_actwidth=twotheta_exact_bin_spacing_deg*bins
        self.set_window(twotheta_window_center_deg=twotheta_window_center_deg,
            twotheta_window_fullwidth_deg=window_actwidth,
            twotheta_output_points=bins)
    
    ## @brief update the dictionary of parameters associated with the given convolver
    #  @param self self
    #  @param convolver the name of the convolver.  name 'global', e.g., attaches to function 'conv_global'
    #  @param kwargs keyword-value pairs to update the convolvers dictionary.
    def set_parameters(self, convolver="global", **kwargs):
        """update the arguments for a specific convolver, by name.  no convolver -> global parameters"""
        self.param_dicts["conv_"+convolver].update(kwargs)
    
    ## @brief get a cached, pre-computed convolver associated with the given parameters,
    #  or a newly zeroed convolver if the cache doesn't contain it. Recycles old cache entries.
    #  @param name the name of the convolver to seek
    #  @param key any hashable object which identifies the parameters for the computation
    #  @param format the type of the array to create, if one is not found.
    #  @return \a flag, which is \a True if valid data were found, or \a False if the returned array is zero,
    #  and \a array, which must be computed by the convolver if \a flag was \a False.
    #
    #  This takes advantage of the mutability of arrays.
    #  When the contents of the array are changed by the convolver, the
    #  cached copy is implicitly updated, so thsat the next time this is called with the same
    #  parameters, it will return the previous array.
    def get_conv(self, name, key, format=numpy.float):
        #either get an old buffer associated with key, with valid convolution kernel,
        #or create a new, empty buffer. If we already have too many buffers, swap one out.
        history=self.convolution_history[name] #previous computed values as a list
        for idx, (k,b) in enumerate(history):
            if k==key:
                history.insert(0,history.pop(idx)) #move to front to mark recently used
                if self.debug_cache: print >> sys.stderr, name, True
                return True, b #True says we got a buffer with valid data
        if len(history)==self.max_history_length:
            buf=history.pop(-1)[1] #re-use oldest buffer
            buf[:]=0
        else:
            buf=numpy.zeros(self.n_omega_points, format)
        history.insert(0,(key,buf))
        if self.debug_cache: print >> sys.stderr, name, False
        return False, buf #False says buffer is empty, need to recompute

    ## @brief scan for functions named conv_xxx, and associated info_xxx entries.
    #  @param self self
    #  @return list of (convolver_xxx, info_xxx) pairs
    def get_convolver_information(self):
        """return a list of convolvers, and what we know about them"""
        info_list=[]
        for k, f in self.convolver_funcs.items():
            info=getattr(self, "info_"+k[5:], {})
            info["docstring"]=f.__doc__
            info_list.append((k,info))

        return info_list

    ## A dictionary of default parameters for the global namespace,
    #  used to seed a GUI which can harvest this for names, descriptions, and initial values
    info_global=dict(
        group_name="Global parameters",
        help="this should be help information",
        param_info=dict(
            twotheta0_deg=("Bragg center of peak (degrees)", 30.0),
            d=("d spacing (m)", 4.00e-10),
            dominant_wavelength=("wavelength of most intense line (m)", 1.5e-10)
        )
    )
    
    ## @brief return a nicely formatted report descibing the current state of this class
    #  @brief self self
    #  @return string of formatted information
    #
    #  this looks for an str_xxx function associated with each conv_xxx name.  If it is found, that
    #  function if called to get the state of conv_xxx.  Otherwise,
    #  this simply formats the dictionary of parameters for the convolver, and uses that.
    def __str__(self):
        #scan for all convolvers, and find if they have str_xxx methods matching
        #if they do, call it for custom output formatting.  Otherwise,
        #just print the dict for each convolver.
        keys=list(self.convolver_funcs.keys())
        keys.sort() #always return info in the same order... maybe later some priority list
        keys.insert(0, keys.pop(keys.index('conv_global'))) #global is always first, anyways!
        strings=["", "***convolver id 0x%08x:"%id(self)]
        for k in keys:
            strfn="str_"+k[5:]
            if hasattr(self, strfn):
                strings.append(getattr(self, strfn)())
            else:
                dd=self.param_dicts["conv_"+k[5:]]
                if dd:
                    strings.append(k[5:]+": "+str(dd))
        return '\n'.join(strings)
        
    ## @brief returns a string representation for the global context.
    #  @param self self
    #  @return report on global parameters.
    def str_global(self):
        self.param_dicts["conv_global"].setdefault("d",0) #in case it's not initialized
        return "global: peak center=%(twotheta0_deg).4f, d=%(d).8g, eq. div=%(equatorial_divergence_deg).3f" % self.param_dicts["conv_global"]
    
    ## @brief the global context isn't really a convolver, returning \a None means ignore result
    #  @param self self
    #  @return \a None, always
    def conv_global(self):
        """a dummy convolver to hold global variables and information"""
        return None
    
    ## @brief the function F0 from the paper.
    #
    #  compute k/sqrt(peakpos-x)+y0 nonzero between outer & inner (inner is closer to peak)
    #  or k/sqrt(x-peakpos)+y0 if reversed (i.e. if outer > peak)
    #  fully evaluated on a specified eps grid, and stuff into destination
    #  @param self self
    #  @param outerbound the edge of the function farthest from the singularity, referenced to epsvals
    #  @param innerbound the edge closest to the singularity, referenced to epsvals
    #  @param epsvals the array of two-theta values or offsets
    #  @param[in,out] destination an array into which final results are summed.  modified in place!
    #  @param peakpos the position of the singularity, referenced to epsvals.
    #  @param y0 the constant offset
    #  @param k the scale factor
    #  @return (\a lower_index, \a upper_index ) python style bounds
    #    for region of \a destination which has been modified.
    def axial_helper(self, outerbound, innerbound, epsvals, destination, peakpos=0, y0=0, k=0):
        if k==0:
            return len(epsvals)//2,len(epsvals)//2+1 #nothing to do, point at the middle
        
        dx=epsvals[1]-epsvals[0] #bin width for normalizer so sum(result*dx)=exact integral
        flip=outerbound > peakpos #flag for whether tail is to the left or right.

        delta1=abs(innerbound-peakpos)
        delta2=abs(outerbound-peakpos)
        #this is the analytic area the function must have, integral(1/sqrt(eps0-eps)) from lower to upper
        exactintegral=2*k*(math.sqrt(delta2)-math.sqrt(delta1))
        exactintegral += y0*(delta2-delta1)
        #exactintegral=max(0,exactintegral) #can never be < 0, beta out of range.
        exactintegral *= 1/dx #normalize so sum is right

        #compute the exact centroid we need for this
        if abs(delta2-delta1) < 1e-12:
            exact_moment1=0
        else:
            exact_moment1=( #simplified from Mathematica FortranForm
                (4*k*(delta2**1.5-delta1**1.5) + 3*y0*(delta2**2-delta1**2))/
                (6.*(2*k*(math.sqrt(delta2)-math.sqrt(delta1)) + y0*(delta2-delta1) ))
            )
            if not flip: exact_moment1=-exact_moment1
        exact_moment1+=peakpos

        # note: because of the way the search is done, this produces a result
        # with a bias of 1/2 channel to the left of where it should be.
        # this is fixed by shifting all the parameters up 1/2 channel
        outerbound+=dx/2; innerbound+=dx/2; peakpos+=dx/2 #fix 1/2 channel average bias from search
        #note: searchsorted(side="left") always returns the position of the bin to the right of the match, or exact bin
        idx0, idx1=epsvals.searchsorted((outerbound, innerbound), side='left')
        
        if abs(outerbound-innerbound)< (2*dx) or abs(idx1-idx0) < 2: #peak has been squeezed out, nothing to do
            #preserve the exact centroid: requires summing into two channels
            #for a peak this narrow, no attempt to preserve the width.
            #note that x1 (1-f1) + (x1+dx) f1 = mu has solution (mu - x1) / dx  = f1
            #thus, we want to sum into a channel that has x1<mu by less than dx, and the one to its right
            idx0 = min(idx0,idx1)-1 #pick left edge and make sure we are past it
            while exact_moment1-epsvals[idx0] > dx:
                #normally only one step max, but make it a loop in case of corner case
                idx0 += 1
            f1=(exact_moment1-epsvals[idx0])/dx
            res=(exactintegral*(1-f1),exactintegral*f1)
            destination[idx0:idx0+2]+=res
            if collect_moment_errors:
                centroid2=(res*epsvals[idx0:idx0+2]).sum()/sum(res)
                moment_list.append((centroid2-exact_moment1)/dx)
            return [idx0, idx0+2] #return collapsed bounds

        if not flip:
            if epsvals[idx0] != outerbound: idx0=max(idx0-1,0)
            idx1=min(idx1+1, len(epsvals))
            sign=1
            deps=self._f0buf[idx0:idx1]
            deps[:]=peakpos
            deps-=epsvals[idx0:idx1]
            deps[-1]=peakpos-min(innerbound, peakpos)
            deps[0]=peakpos-outerbound
        else:
            idx0, idx1 = idx1, idx0
            if epsvals[idx0] != innerbound: idx0=max(idx0-1,0)
            idx1=min(idx1+1, len(epsvals))
            sign=-1
            deps=self._f0buf[idx0:idx1]
            deps[:]=epsvals[idx0:idx1]
            deps-=peakpos
            deps[0]=max(innerbound, peakpos)-peakpos
            deps[-1]=outerbound-peakpos

        dx0=abs(deps[1]-deps[0])
        dx1=abs(deps[-1]-deps[-2])
        
        #make the numerics accurate: compute average on each bin, which is
        #integral of 1/sqrt = 2*sqrt, then difference integral
        intg=numpy.sqrt(deps,deps) #do it in place, return value is actually deps too
        intg *= 2*k*sign

        intg[:-1]-=intg[1:] #do difference in place, running forward to avoid self-trampling
        intg[1:-2] += y0*dx #add constant
        #handle narrowed bins on ends carefully
        intg[0] += y0*dx0
        intg[-2] += y0*dx1
        
        if min(intg[:-1]) < -1e-10*max(intg[:-1]): #intensities are never less than zero!
            print( "bad parameters:", (5*"%10.4e ") %(peakpos, innerbound, outerbound, k, y0))
            print( len(intg), intg[:-1])
            raise ValueError("Bad axial helper parameters")
                
        #now, make sure the underlying area is the exactly correct
        #integral, without bumps due to discretizing the grid.
        intg *= (exactintegral/(intg[:-1].sum()))

        destination[idx0:idx1-1]+=intg[:-1]
        
        ## This is purely for debugging.  If collect_moment_errors is \a True,
        #  compute exact vs. approximate moments.
        if collect_moment_errors:
            centroid2=(intg[:-1]*epsvals[idx0:idx1-1]).sum()/intg[:-1].sum()
            moment_list.append((centroid2-exact_moment1)/dx)

        return [idx0, idx1-1] #useful info for peak position

    ## @brief return the \a I2 function
    #  @param self self
    #  @param Lx length of the xray filament
    #  @param Ls length of the sample
    #  @param Lr length of the receiver slit
    #  @param R diffractometer length, assumed symmetrical
    #  @param twotheta angle, in radians, of the center of the computation
    #  @param beta offset angle
    #  @param epsvals array of offsets from center of computation, in radians
    #  @return ( \a epsvals, \a idxmin, \a idxmax, \a I2p, \a I2m ).
    #   \a idxmin and \a idxmax are the full python-style bounds of the non-zero region of \a I2p and \a I2m.
    #   \a I2p and \a I2m are I2+ and I2- from the paper, the contributions to the intensity.
    def full_axdiv_I2(self, Lx=None, Ls=None, Lr=None, R=None, twotheta=None, beta=None,
            epsvals=None
        ):
        
        from math import sin, cos, tan, pi
        beta1=(Ls-Lx)/(2*R) #Ch&Co after eq. 15abcd
        beta2=(Ls+Lx)/(2*R) #Ch&Co after eq. 15abcd, corrected by KM
        
        eps0=beta*beta*math.tan(twotheta)/2 #after eq. 26 in Ch&Co
        
        if -beta2 <= beta < beta1:
            z0p=Lx/2+beta*R*(1+1/cos(twotheta))
        elif beta1 <= beta <= beta2:
            z0p=Ls/2+beta*R/cos(twotheta)
        
        if -beta2 <= beta <= -beta1:
            z0m=-Ls/2+beta*R/cos(twotheta)
        elif -beta1 < beta <= beta2:
            z0m=-Lx/2+beta*R*(1+1/cos(twotheta))

        epsscale=tan(pi/2-twotheta)/(2*R*R) #=cotan(twotheta)...

        eps1p=(eps0-epsscale*((Lr/2)-z0p)**2) #Ch&Co 18a&18b, KM sign correction
        eps2p=(eps0-epsscale*((Lr/2)-z0m)**2)
        eps2m=(eps0-epsscale*((Lr/2)+z0p)**2) #reversed eps2m and eps1m per KM R
        eps1m=(eps0-epsscale*((Lr/2)+z0m)**2) #flip all epsilons

        if twotheta>pi/2: #this set of inversions from KM 'R' code, simplified here
            eps1p, eps2p, eps1m, eps2m = eps1m, eps2m, eps1p, eps2p

        #identify ranges per Ch&Co 4.2.2 and table 1 and select parameters
        #note table 1 is full of typos, but the minimized
        #tests from 4.2.2 with redundancies removed seem fine.
        if Lr > z0p-z0m:
            if  z0p <= Lr/2 and z0m > -Lr/2: #beam entirely within slit
                rng=1; ea=eps1p; eb=eps2p; ec=eps1m; ed=eps2m
            elif (z0p > Lr/2 and z0m < Lr/2) or (z0m < -Lr/2 and z0p > -Lr/2):
                rng=2; ea=eps2p; eb=eps1p; ec=eps1m; ed=eps2m
            else:
                rng=3; ea=eps2p; eb=eps1p; ec=eps1m; ed=eps2m
        else:
            if z0m < -Lr/2 and z0p > Lr/2: #beam hanging off both ends of slit, peak centered
                rng=1; ea=eps1m; eb=eps2p; ec=eps1p; ed=eps2m
            elif (-Lr/2 < z0m < Lr/2 and z0p > Lr/2) or (-Lr/2 < z0p < Lr/2 and z0m < -Lr/2): #one edge of beam within slit
                rng=2; ea=eps2p; eb=eps1m; ec=eps1p; ed=eps2m
            else:
                rng=3; ea=eps2p; eb=eps1m; ec=eps1p; ed=eps2m

        #now, evaluate function on bounds in table 1 based on ranges
        #note: because of a sign convention in epsilon, the bounds all get switched
        
        def F1(dst, lower, upper, eea, eeb): #define them in our namespace so they inherit ea, eb, ec, ed, etc.
            return self.axial_helper(destination=dst,
                    innerbound=upper, outerbound=lower,
                    epsvals=epsvals, peakpos=eps0,
                    k=math.sqrt(abs(eps0-eeb))-math.sqrt(abs(eps0-eea)), y0=0
                )
        def F2(dst, lower, upper, eea):
            return self.axial_helper(destination=dst,
                    innerbound=upper, outerbound=lower,
                    epsvals=epsvals, peakpos=eps0,
                    k=math.sqrt(abs(eps0-eea)), y0=-1
                )
        def F3(dst, lower, upper, eea):
            return self.axial_helper(destination=dst,
                    innerbound=upper, outerbound=lower,
                    epsvals=epsvals, peakpos=eps0,
                    k=math.sqrt(abs(eps0-eea)), y0=+1
                )
        def F4(dst, lower, upper, eea):
            #just like F2 but k and y0 negated
            return self.axial_helper(destination=dst,
                    innerbound=upper, outerbound=lower,
                    epsvals=epsvals, peakpos=eps0,
                    k=-math.sqrt(abs(eps0-eea)), y0=+1
                )

        I2p=self._I2p
        I2p[:]=0
        I2m=self._I2m
        I2m[:]=0
        
        indices=[]
        if rng==1:
            indices+=F1(dst=I2p, lower=ea, upper=eps0, eea=ea, eeb=eb)
            indices+=F2(dst=I2p, lower=eb, upper=ea,   eea=eb)
            indices+=F1(dst=I2m, lower=ec, upper=eps0, eea=ec, eeb=ed)
            indices+=F2(dst=I2m, lower=ed, upper=ec,   eea=ed)
        elif rng==2:
            indices+=F2(dst=I2p, lower=ea, upper=eps0, eea=ea)
            indices+=F3(dst=I2m, lower=eb, upper=eps0, eea=ea)
            indices+=F1(dst=I2m, lower=ec, upper=eb  , eea=ec, eeb=ed)
            indices+=F2(dst=I2m, lower=ed, upper=ec  , eea=ed)
        elif rng==3:
            indices+=F4(dst=I2m, lower=eb, upper=ea  , eea=ea)
            indices+=F1(dst=I2m, lower=ec, upper=eb  , eea=ec, eeb=ed)
            indices+=F2(dst=I2m, lower=ed, upper=ec  , eea=ed)

        idxmin=min(indices)
        idxmax=max(indices)

        return epsvals, idxmin, idxmax, I2p,I2m

    ## @brief carry out the integral of \a I2 over \a beta and the Soller slits.
    #  @param self self
    #  @param Lx length of the xray filament
    #  @param Ls length of the sample
    #  @param Lr length of the receiver slit
    #  @param R the (assumed symmetrical) diffractometer radius
    #  @param twotheta angle, in radians, of the center of the computation
    #  @param epsvals array of offsets from center of computation, in radians
    #  @param sollerIdeg the full-width (both sides) cutoff angle of the incident Soller slit
    #  @param sollerDdeg the full-width (both sides) cutoff angle of the detector Soller slit
    #  @param nsteps the number of subdivisions for the integral
    #  @param axDiv not used
    #  @return the accumulated integral, a copy of a persistent buffer \a _axial
    def full_axdiv_I3(self, Lx=None, Ls=None, Lr=None, R=None,
        twotheta=None,
        epsvals=None, sollerIdeg=None, sollerDdeg=None, nsteps=10, axDiv=""):
        from math import sin, cos, tan, pi
        
        beta2=(Ls+Lx)/(2*R) #Ch&Co after eq. 15abcd, corrected by KM

        if sollerIdeg is not None:
            solIrad=sollerIdeg*pi/180/2
            #solIfunc=numpy.frompyfunc(lambda x: max(0, 1.0-abs(x/solIrad)),1,1)
            #note: solIfunc is only called with a scalar, no need for numpy-ized version, really
            def solIfunc(x):
                return numpy.clip(1.0-abs(x/solIrad),0,1)
            beta2=min(beta2, solIrad) #no point going beyond Soller
        else:
            def solIfunc(x):
                return numpy.ones_like(x)
        if sollerDdeg is not None:
            solDrad=sollerDdeg*pi/180/2
            #solDfunc=numpy.frompyfunc(lambda x: max(0, 1.0-abs(x/solDrad)),1,1)
            def solDfunc(x):
                return numpy.clip(1.0-abs(x/solDrad),0,1)
        else:
            def solDfunc(x):
                return numpy.ones_like(x)

        accum=self._axial
        accum[:]=0
        
        if twotheta > pi/2:
            tth1=pi-twotheta
        else:
            tth1=twotheta

        for iidx in range(nsteps):
            beta=beta2*iidx/float(nsteps)
            
            eps, idxmin, idxmax, I2p, I2m=self.full_axdiv_I2(
                Lx=Lx,
                Lr=Lr,
                Ls=Ls,
                beta=beta, R=R,
                twotheta=twotheta, epsvals=epsvals
            )
            
            eps0=beta*beta*math.tan(twotheta)/2 #after eq. 26 in Ch&Co
            
            gamma0=beta/cos(tth1)
            deps=self._f0buf[idxmin:idxmax]
            deps[:]=eps0
            deps-=epsvals[idxmin:idxmax]
            deps *= 2*math.tan(twotheta)
            #check two channels on each end for negative argument.
            deps[-1]=max(deps[-1],0)
            deps[0]=max(deps[0],0)
            if len(deps) >= 2:
                deps[-2]=max(deps[-2],0)
                deps[1]=max(deps[1],0)
            
            gamarg=numpy.sqrt(deps, deps) #do sqrt in place for speed
            #still need to convert these to in-place
            gamp=gamma0+gamarg
            gamm=gamma0-gamarg
            
            if iidx==0 or iidx==nsteps-1: weight=1.0 #trapezoidal rule weighting
            else: weight=2.0

            #sum into the accumulator only channels which can be non-zero
            #do scaling in-place to save a  lot of slow array copying
            I2p[idxmin:idxmax] *= solDfunc(gamp)
            I2p[idxmin:idxmax] *=(weight*solIfunc(beta))
            accum[idxmin:idxmax]+=I2p[idxmin:idxmax]
            I2m[idxmin:idxmax] *= solDfunc(gamm)
            I2m[idxmin:idxmax] *=(weight*solIfunc(beta))
            accum[idxmin:idxmax]+=I2m[idxmin:idxmax]

        #keep this normalized
        K = 2 * R*R * abs(tan(twotheta))
        accum *= K
        
        return accum
    
    ## @brief compute the Fourier transform of the axial divergence component
    #  @param self self
    #  @return the transform buffer, or \a None if this is being ignored
    def conv_axial(self):
        me=self.get_function_name() #the name of this convolver,as a string
        if self.param_dicts[me].get("axDiv",None) is None:
            return None
        kwargs={}
        kwargs.update(self.param_dicts[me]) #get all of our parameters
        kwargs.update(self.param_dicts["conv_global"])
        if "equatorial_divergence_deg" in kwargs:
            del kwargs["equatorial_divergence_deg"] #not used
        
        flag, axfn = self.get_conv(me, kwargs, numpy.complex)
        if flag: return axfn #already up to date if first return is True
        
        xx=type("data",(), kwargs)
        if xx.axDiv!="full" or xx.twotheta0_deg==90.0: #no axial divergence, transform of delta fn
            axfn[:]=1
            return axfn
        else:
            axbuf=self.full_axdiv_I3(
                nsteps=xx.n_integral_points,
                epsvals=self.epsilon,
                Lx=xx.slit_length_source,
                Lr=xx.slit_length_target,
                Ls=xx.length_sample,
                sollerIdeg=xx.angI_deg,
                sollerDdeg=xx.angD_deg,
                R=xx.diffractometer_radius,
                twotheta=xx.twotheta0
            )
        axfn[:]=best_rfft(axbuf)
        
        return axfn

    ## @brief compute the Fourier transform of the rectangular tube tails function
    #  @param self self
    #  @return the transform buffer, or \a None if this is being ignored
    def conv_tube_tails(self):
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs={}
        kwargs.update(self.param_dicts[me]) #get all of our parameters
        if not kwargs:
            return None #no convolver
        #we also need the diffractometer radius from the global space
        kwargs["diffractometer_radius"]=self.param_dicts["conv_global"]["diffractometer_radius"]
        flag, tailfn = self.get_conv(me, kwargs, numpy.complex)
        if flag: return tailfn #already up to date
        
        #tube_tails is (main width, left width, right width, intensity),
        #so peak is raw peak + tophat centered at (left width+ right width)
        #with area intensity*(right width-left width)/main_width
        #check this normalization!
        #note: this widths are as defined by Topas... really I think it should be
        #x/(2*diffractometer_radius) since the detector is 2R from the source,
        #but since this is just a fit parameter, we'll defin it as does Topas
        xx=type("data",(), kwargs) #allow dotted notation
        
        tail_eps=(xx.tail_right-xx.tail_left)/xx.diffractometer_radius
        main_eps=xx.main_width/xx.diffractometer_radius
        tail_center=(xx.tail_right+xx.tail_left)/xx.diffractometer_radius/2.0
        tail_area=xx.tail_intens*(xx.tail_right-xx.tail_left)/xx.main_width

        cb1=self._cb1
        rb1=self._rb1
        
        cb1.real=0
        cb1.imag = self.omega_vals
        cb1.imag *= tail_center #sign is consistent with Topas definition
        numpy.exp(cb1, tailfn) #shifted center, computed into tailfn
        
        rb1[:]=self.omega_vals
        rb1 *= (tail_eps/2/math.pi)
        rb1=numpy.sinc(rb1)
        tailfn *= rb1
        tailfn *= tail_area #normalize
        
        rb1[:]=self.omega_vals
        rb1 *= (main_eps/2/math.pi)
        rb1=numpy.sinc(rb1)
        tailfn += rb1 #add central peak
        return tailfn
    
    ## @brief a utility to compute a transformed tophat function and save it in a convolver buffer
    #  @param self self
    #  @param name the name of the convolver cache buffer to update
    #  @param width the width in 2-theta space of the tophat
    #  @return the updated convolver buffer, or \a None if the width was \a None
    def general_tophat(self, name="", width=None):
        """handle all centered top-hats"""
        if width is None:
            return #no convolver
        flag, conv = self.get_conv(name, width, numpy.float)
        if flag: return conv #already up to date
        rb1=self._rb1
        rb1[:]=self.omega_vals
        rb1 *= (width/2/math.pi)
        conv[:]=numpy.sinc(rb1)
        return conv

    ## A dictionary of default parameters for conv_emissions,
    #  used to seed a GUI which can harvest this for names, descriptions, and initial values
    info_emission=dict(
        group_name="Incident beam and crystal size",
        help="this should be help information",
        param_info=dict(
            emiss_wavelengths=("wavelengths (m)", (1.58e-10,)),
            emiss_intensities=("relative intensities", (1.00,)),
            emiss_lor_widths=("Lorenztian emission fwhm (m)",(1e-13,)),
            emiss_gauss_widths=("Gaussian emissions fwhm (m)", (1e-13,)),
            crystallite_size_gauss=("Gaussian crystallite size fwhm (m)", 1e-6),
            crystallite_size_lor=("Lorentzian crystallite size fwhm (m)", 1e-6),
        )
    )
    
    ## @brief format the emission spectrum and crystal size information
    #  @param self self
    #  @return the formatted information
    def str_emission(self):
        dd=self.param_dicts["conv_emission"]
        if not dd: return "No emission spectrum"
        dd.setdefault("crystallite_size_lor",1e10)
        dd.setdefault("crystallite_size_gauss",1e10)
        dd.setdefault("strain_lor",0)
        dd.setdefault("strain_gauss",0)
        xx=type("data",(), dd)
        spect=numpy.array((
            xx.emiss_wavelengths, xx.emiss_intensities,
            xx.emiss_lor_widths, xx.emiss_gauss_widths))
        spect[0]*=1e10*self.length_scale_m #convert to Angstroms, like Topas
        spect[2]*=1e13*self.length_scale_m #milli-Angstroms
        spect[3]*=1e13*self.length_scale_m #milli-Angstroms
        nm=1e9*self.length_scale_m
        items=["emission and broadening:"]
        items.append("spectrum=\n"+str(spect.transpose()))
        items.append("crystallite_size_lor (nm): %.5g" % (xx.crystallite_size_lor*nm))
        items.append("crystallite_size_gauss (nm): %.5g" % (xx.crystallite_size_gauss*nm))
        items.append("strain_lor: %.5g" % xx.strain_lor)
        items.append("strain_gauss: %.5g" % xx.strain_gauss)
        return '\n'.join(items)
    
    ## @brief compute the emission spectrum and (for convenience) the particle size widths
    #  @param self self
    #  @return the convolver for the emission and particle sizes
    #  @note the particle size and strain stuff here is just to be consistent with \a Topas
    #  and to be vaguely efficient about the computation, since all of these
    #  have the same general shape.
    def conv_emission(self):
        """handle emission spectrum and crystal size together, since it makes sense"""
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs={}
        kwargs.update(self.param_dicts[me]) #get all of our parameters
        kwargs.update(self.param_dicts["conv_global"])
        #if the crystallite size and strain parameters are not set, set them to
        #values that make their corrections disappear
        kwargs.setdefault("crystallite_size_lor",1e10)
        kwargs.setdefault("crystallite_size_gauss",1e10)
        kwargs.setdefault("strain_lor",0)
        kwargs.setdefault("strain_gauss",0)

        #convert arrays to lists for key checking
        key={}
        key.update(kwargs)
        for k,v in key.items():
            if hasattr(v,'tolist'):
                key[k]=v.tolist()
    
        flag, emiss = self.get_conv(me, key, numpy.complex)
        if flag: return emiss #already up to date
        
        xx=type("data", (), kwargs) #make it dot-notation accessible
        
        epsilon0s = 2 * numpy.arcsin(xx.emiss_wavelengths / (2.0 * xx.d)) - xx.twotheta0
        theta=xx.twotheta0/2
        ## Emission profile FWHM + crystallite broadening (scale factors are Topas choice!) (Lorentzian)
        #note:  the strain broadenings in Topas are expressed in degrees 2theta, must convert to radians(theta) with pi/360
        widths = (
            (xx.emiss_lor_widths/xx.emiss_wavelengths) * math.tan(theta) +
            xx.strain_lor*(math.pi/360) * math.tan(theta) +
            (xx.emiss_wavelengths / (2*xx.crystallite_size_lor * math.cos(theta)))
        )
        #save weighted average width for future reference in periodicity fixer
        self.lor_widths[me]=sum(widths*xx.emiss_intensities)/sum(xx.emiss_intensities)
        #gaussian bits add in quadrature
        gfwhm2s=(
            ((2*xx.emiss_gauss_widths/xx.emiss_wavelengths) * math.tan(theta))**2 +
            (xx.strain_gauss*(math.pi/360) * math.tan(theta))**2 +
            (xx.emiss_wavelengths / (xx.crystallite_size_gauss * math.cos(theta)))**2
        )
        
        # note that the Fourier transform of a lorentzian with FWHM 2a
        # is exp(-abs(a omega))
        #now, the line profiles in Fourier space have to have phases
        #carefully handled to put the lines in the right places.
        #note that the transform of f(x+dx)=exp(i omega dx) f~(x)
        omega_vals=self.omega_vals
        for wid, gfwhm2, eps, intens in zip(widths, gfwhm2s, epsilon0s, xx.emiss_intensities):
            xvals=numpy.clip(omega_vals*(-wid),-100,0)
            sig2=gfwhm2/(8*math.log(2.0)) #convert fwhm**2 to sigma**2
            gxv=numpy.clip((sig2/-2.0)*omega_vals*omega_vals,-100,0)
            emiss += numpy.exp(xvals+gxv+complex(0,-eps)*omega_vals)*intens
        return emiss
    
    ## @brief compute the convolver for the flat-specimen correction
    #  @param self self
    #  @return the convolver
    def conv_flat_specimen(self):
        """handle flat specimen"""
        me=self.get_function_name() #the name of this convolver,as a string
        equatorial_divergence_deg=self.param_dicts["conv_global"].get("equatorial_divergence_deg",None)
        if not equatorial_divergence_deg: return None
        twotheta0=self.param_dicts["conv_global"]["twotheta0"]
        key=(twotheta0, equatorial_divergence_deg)
        flag, conv = self.get_conv(me, key, numpy.complex)
        if flag: return conv #already up to date
        
        #Flat-specimen error, from Cheary, Coelho & Cline 2004 NIST eq. 9 & 10
        #compute epsm in radians from eq. divergence in degrees
        #to make it easy to use the axial_helper to compute the function
        epsm = ((equatorial_divergence_deg*math.pi/180)**2)/math.tan(twotheta0/ 2.0)/2.0
        eqdiv=self._epsb2
        eqdiv[:]=0
        dtwoth=(self.twothetasamples[1]-self.twothetasamples[0])
        idx0, idx1=self.axial_helper(destination=eqdiv,
            outerbound=-epsm,
            innerbound=0,
            epsvals=self.epsilon,
            peakpos=0, k=dtwoth/(2.0*math.sqrt(epsm)))

        conv[:]=best_rfft(eqdiv)
        conv[1::2] *= -1 #flip center
        return conv
    
    ## @brief compute the sample transparency correction, including the finite-thickness version
    #  @param self self
    #  @return the convolver
    def conv_absorption(self):
        """handle transparency"""
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs={}
        kwargs.update(self.param_dicts[me]) #get all of our parameters
        if not kwargs: return None
        kwargs["twotheta0"]=self.param_dicts["conv_global"]["twotheta0"]
        kwargs["diffractometer_radius"]=self.param_dicts["conv_global"]["diffractometer_radius"]

        flag, conv = self.get_conv(me, kwargs, numpy.complex)
        if flag: return conv #already up to date
        xx=type("data", (), kwargs) #make it dot-notation accessible
        
        #absorption, from Cheary, Coelho & Cline 2004 NIST eq. 12,
        #EXCEPT delta = 1/(2*mu*R) instead of 2/(mu*R)
        #from Mathematica, unnormalized transform is
        #(1-exp(epsmin*(i w + 1/delta)))/(i w + 1/delta)
        delta = math.sin(xx.twotheta0)/(2*xx.absorption_coefficient*xx.diffractometer_radius)
        #arg=(1/delta)+complex(0,-1)*omega_vals
        cb=self._cb1
        cb.imag=self.omega_vals
        cb.imag*=-1
        cb.real=1/delta
        numpy.reciprocal(cb,conv) #limit for thick samples=1/(delta*arg)
        conv *= 1.0/delta #normalize
        if kwargs.get("sample_thickness", None) is not None: #rest of transform of function with cutoff
            epsmin = -2.0*xx.sample_thickness*math.cos(xx.twotheta0/2.0)/xx.diffractometer_radius
            cb*=epsmin
            numpy.expm1(cb,cb)
            cb*=-1
            conv *= cb
        return conv
    
    ## @brief compute the peak shift due to sample displacement and the \a 2theta zero
    #  @param self self
    #  @return the convolver
    def conv_displacement(self):
        """handle displacements and peak offset from center"""
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs=self.param_dicts[me]
        twotheta0=self.param_dicts["conv_global"]["twotheta0"]
        diffractometer_radius=self.param_dicts["conv_global"]["diffractometer_radius"]
        specimen_displacement=kwargs.get("specimen_displacement", 0.0)
        if specimen_displacement is None: specimen_displacement=0.0
        zero_error_deg=kwargs.get("zero_error_deg", 0.0)
        if zero_error_deg is None: zero_error_deg=0.0
                
        flag, conv = self.get_conv(me,
            (twotheta0, diffractometer_radius,
            specimen_displacement, zero_error_deg),
            numpy.complex)
        if flag: return conv#already up to date

        delta=-2*math.cos(twotheta0/2.0)*specimen_displacement/diffractometer_radius
        conv.real=0
        conv.imag=self.omega_vals
        #convolver *= numpy.exp(complex(0,-delta-zero_error_deg*pi/180.0)*omega_vals)
        conv.imag *= (-delta-zero_error_deg*math.pi/180.0-(twotheta0-self.twotheta_window_center))
        numpy.exp(conv,conv)
        return conv

    ## @brief compute the rectangular convolution for the receiver slit or SiPSD pixel size
    #  @param self self
    #  @return the convolver
    def conv_receiver_slit(self):
        """ receiver slit width or si psd pixel width"""
        me=self.get_function_name() #the name of this convolver,as a string
        #The receiver slit convolution is a top-hat of angular half-width
        #a=(slit_width/2)/diffractometer_radius
        #which has Fourier transform of sin(a omega)/(a omega)
        #NOTE! numpy's sinc(x) is sin(pi x)/(pi x), not sin(x)/x
        if self.param_dicts[me].get("slit_width",None) is None:
            return None
        
        epsr=(self.param_dicts["conv_receiver_slit"]["slit_width"]/
            self.param_dicts["conv_global"]["diffractometer_radius"])
        return self.general_tophat(me, epsr)

    ## @brief compute the convolver for the integral of defocusing of the face of an Si PSD
    #  @param self self
    #  @return the convolver
    def conv_si_psd(self):
        #omega offset defocussing from Cheary, Coelho & Cline 2004 eq. 15
        #expressed in terms of a Si PSD looking at channels with vertical offset
        #from the center between psd_window_lower_offset and psd_window_upper_offset
        #do this last, because we may ultimately take a list of bounds,
        #and return a list of convolutions, for efficiency
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs={}
        kwargs.update(self.param_dicts[me]) #get all of our parameters
        if not kwargs: return None
        kwargs.update(self.param_dicts["conv_global"])

        flag, conv = self.get_conv(me, kwargs, numpy.float)
        if flag: return conv #already up to date
        
        xx=type("data",(), kwargs)
        
        if not xx.equatorial_divergence_deg or not xx.si_psd_window_bounds:
            #if either of these is zero or None, convolution is trivial
            conv[:]=1
            return conv
        
        psd_lower_window_pos, psd_upper_window_pos = xx.si_psd_window_bounds
        dthl=psd_lower_window_pos/xx.diffractometer_radius
        dthu=psd_upper_window_pos/xx.diffractometer_radius
        alpha=xx.equatorial_divergence_deg*math.pi/180
        argscale=alpha/(2.0*math.tan(xx.twotheta0/2))
        from scipy.special import sici #get the sine and cosine integral
        #WARNING si(x)=integral(sin(x)/x), not integral(sin(pi x)/(pi x))
        #i.e. they sinc function is not consistent with the si function
        #whence the missing pi in the denominator of argscale
        rb1=self._rb1
        rb2=self._rb2
        rb3=self._rb3
        rb1[:]=self.omega_vals
        rb1 *= argscale*dthu
        sici(rb1,conv,rb3) #gets both sine and cosine integral, si in conv
        if dthl:  #no need to do this if the lower bound is 0
            rb1[:]=self.omega_vals
            rb1 *= argscale*dthl
            sici(rb1,rb2,rb3) #gets both sine and cosine integral, si in rb2
            conv-=rb2
        conv[1:] /= self.omega_vals[1:]
        conv *= 1/argscale
        conv[0] = dthu-dthl #fix 0/0 with proper area
        return conv
    
    ## @brief compute the convolver to smooth the final result with a Gaussian before downsampling.
    #  @param self self
    #  @return the convolver
    def conv_smoother(self):
        #create a smoother for output result, independent of real physics, if wanted
        me=self.get_function_name() #the name of this convolver,as a string
        if not self.output_gaussian_smoother_bins_sigma: return # no smoothing
        flag, buf=self.get_conv(me, self.output_gaussian_smoother_bins_sigma,
            format=numpy.float)
        if flag: return buf #already computed
        buf[:]=self.omega_vals
        buf*=(self.output_gaussian_smoother_bins_sigma*(
            self.twothetasamples[1]-self.twothetasamples[0]))
        buf *= buf
        buf *=-0.5
        numpy.exp(buf, buf)
        return buf
    
    ## @brief execute all the convolutions; if convolver_names is None, use everything
    #    we have, otherwise, use named convolutions.
    #  @param self self
    #  @param convolver_names a list of convolvers to select. If \a None, use all found convolvers.
    #  @param compute_derivative if \a True, also return d/dx(function) for peak position fitting
    #  @return a profile_data object with much information about the peak
    def compute_line_profile(self, convolver_names=None, compute_derivative=False, return_convolver=False):
        
        #create a function which is the Fourier transform of the
        #combined convolutions of all the factors
        from numpy import sin, cos, tan, arcsin as asin, arccos as acos
        from math import pi

        anglemode=self.anglemode #="d" if we are using 'd' space, "twotheta" if using twotheta
        
        ## the rough center of the spectrum, used for things which need it. Copied from global convolver.
        self.dominant_wavelength=dominant_wavelength=self.param_dicts["conv_global"].get("dominant_wavelength",None)
        
        if anglemode=="twotheta":
            twotheta0_deg=self.param_dicts["conv_global"]["twotheta0_deg"]
            twotheta0=math.radians(twotheta0_deg)
            d=dominant_wavelength/(2*sin(twotheta0/2.0))
        else:
            d=self.param_dicts["conv_global"]["d"]
            twotheta0 = 2 * asin(dominant_wavelength/ (2.0 * d))
            twotheta0_deg=math.degrees(twotheta0)
        
        self.set_parameters(d=d,twotheta0=twotheta0,twotheta0_deg=twotheta0_deg) #set these in global namespace

        if convolver_names is None:
            convolver_names=self.convolver_funcs.keys() #get all names
        
        #run through the name list, and call the convolver to harvest its result
        conv_list=[self.convolver_funcs[x]() for x in convolver_names]

        #now, multiply everything together
        convolver=self._cb1 #get a complex scratch buffer
        convolver[:]=1 #initialize
        for c in conv_list: #accumulate product
            if c is not None: convolver *= c

        if convolver[1].real > 0: #recenter peak!
            convolver[1::2]*=-1
        
        peak=best_irfft(convolver)

        #now, use the trick from Mendenhall JQSRT Voigt paper to remove periodic function correction
        #JQSRT 105 number 3 July 2007 p. 519 eq. 7
        emiss_intensities=self.param_dicts["conv_emission"]["emiss_intensities"]
        correction_width=2*sum(self.lor_widths.values()) #total lor widths, created by the various colvolvers
        
        d2p=2.0*pi/self.window_fullwidth
        alpha=correction_width/2.0 #be consistent with convolver
        mu=(peak*self.twothetasamples).sum()/peak.sum() #centroid
        dx=self.twothetasamples-mu
        eps_corr1=(math.sinh(d2p*alpha)/self.window_fullwidth)/(math.cosh(d2p*alpha)-numpy.cos(d2p*dx))
        eps_corr2=(alpha/pi)/(dx*dx+alpha*alpha)
        corr=(convolver[0].real/numpy.sum(eps_corr2))*(eps_corr1-eps_corr2)
        peak -= corr
        
        peak*=self.window_fullwidth/(self.twotheta_output_points/self.oversampling) #scale to area
        
        if compute_derivative:
            #this is useful 
            convolver *= self.omega_vals
            convolver *= complex(0,1)
            deriv = best_irfft(convolver)
            deriv*=self.window_fullwidth/(self.twotheta_output_points/self.oversampling)
            deriv=deriv[::self.oversampling]
        else:
            deriv=None
        
        result=profile_data(twotheta0_deg=twotheta0*180/math.pi,
            twotheta=self.twothetasamples[::self.oversampling],
            omega_inv_deg=self.omega_inv_deg[:self.twotheta_output_points//2+1],
            twotheta_deg=self.twothetasamples_deg[::self.oversampling],
            peak=peak[::self.oversampling],
            derivative=deriv
        )

        if return_convolver:
            result.add_symbol(convolver=convolver[:self.twotheta_output_points//2+1])
            
        return result
    
    ## @brief do some cleanup to make us more compact;
    #  Instance can no longer be used after doing this, but can be pickled.
    #  @param self self
    def self_clean(self):
        clean=self._clean_on_pickle
        pd=dict()
        pd.update(self.__dict__)        
        for thing in pd.keys():
            x=getattr(self, thing)
            if id(x) in clean:
                delattr(self,thing)
        #delete extra attributes cautiously, in case we have already been cleaned
        for k in ('convolver_funcs','convolvers','factors','convolution_history'):
            if pd.pop(k,None) is not None:
                delattr(self,k)

    ## @brief return information for pickling
    #  @param self self
    #  @return dictionary of sufficient information to reanimate instance.
    #
    #  Removes transient data from cache of shadow copy so resulting object is fairly compact.
    #  This does not affect the state of the actual instance.
    def __getstate__(self):
        #do some cleanup on state before we get pickled
        #(even if main class not cleaned)
        clean=self._clean_on_pickle
        pd=dict()
        pd.update(self.__dict__)
        for thing in pd.keys():
            x=getattr(self, thing)
            if id(x) in clean:
                del pd[thing]
        #delete extra attributes cautiously, in case we have alread been cleaned
        for k in ('convolver_funcs','convolvers','factors','convolution_history'):
            pd.pop(k,None)
        return pd

    ## @brief reconstruct class from pickled information
    #  @param self an empty class instance
    #  @param setdict dictionary from FP_profile.__getstate__()
    #
    #  This rebuilds the class instance so it is ready to use on unpickling.
    def __setstate__(self, setdict):
        self.__init__(anglemode=setdict["anglemode"],
            output_gaussian_smoother_bins_sigma=setdict["output_gaussian_smoother_bins_sigma"],
            oversampling=setdict["oversampling"]
        )
        for k,v in setdict.items():
            setattr(self,k,v)
        self.set_window(
            twotheta_window_center_deg=self.twotheta_window_center_deg,
            twotheta_window_fullwidth_deg=self.twotheta_window_fullwidth_deg,
            twotheta_output_points=self.twotheta_output_points
        )
        self.lor_widths=setdict["lor_widths"] #override clearing of this by set_window

##
#  @brief A single function call interface to an instance of the class.
#  This is only intended for debugging and simple calculation.  Since
#  it re-creates the class each time, there is no cacheing of results.  It is
#  fairly inefficient to do stuff this way.
#  @param d d-spacing, if present.  Either a d-spacing and wavelength or an angle must be provided
#  @param twotheta_deg the peak center, if a d-spacing is not given.
#  @return an object with much information about a line profile.
def fourier_line_profile(
        d=None,  twotheta_deg=None,
        crystallite_size_lor=None,
        crystallite_size_gauss=None,
        strain_lor=0.0, strain_gauss=0.0,
        wavelength=None,
        slit_width=None,
        slit_length_target=None,
        slit_length_source=None,
        length_sample=None,
        angI_deg = None,
        angD_deg = None,
        axDiv = "simple",
        diffractometer_radius=None,
        si_psd_window_bounds=None,
        equatorial_divergence_deg=None, absorption_coefficient=None,
        sample_thickness=None,
        specimen_displacement=None, zero_error_deg=None, 
        target_width=None, specimen_tilt=None, defocus_delta_omega_deg=None,
        mat_wavelengths=None, mat_lor_widths=None, mat_gauss_widths=None, mat_intensities=None,
        tube_tails=None,
        window_fullwidth_deg=1.,
        twotheta_output_points=1000,
        n_integral_points=10,
        output_gaussian_smoother_bins_sigma=1.0,
        oversampling=10
    ):
    if twotheta_deg is not None: #give priority to angle
        anglemode="twotheta"
        d=0.0
    else:
        anglemode="d"
        twotheta_deg = 2 * math.degrees(math.asin(wavelength/ (2.0 * d)))
    
    p=FP_profile(anglemode=anglemode,
        output_gaussian_smoother_bins_sigma=output_gaussian_smoother_bins_sigma,
        oversampling=oversampling
    )
    #put the compute window in the right place, using old convention
    #centering window on line
    p.set_window(
        twotheta_output_points=twotheta_output_points,
        twotheta_window_center_deg=twotheta_deg,
        twotheta_window_fullwidth_deg=window_fullwidth_deg,
    )

    #set parameters which are shared by many things
    p.set_parameters(d=d,twotheta0_deg=twotheta_deg, dominant_wavelength=wavelength,
        equatorial_divergence_deg=equatorial_divergence_deg,
        diffractometer_radius=diffractometer_radius)
    #set parameters for each convolver
    p.set_parameters(convolver="emission",
        emiss_wavelengths=mat_wavelengths,
        emiss_intensities=mat_intensities,
        emiss_gauss_widths=mat_gauss_widths,
        emiss_lor_widths=mat_lor_widths,
        crystallite_size_gauss=crystallite_size_gauss,
        crystallite_size_lor=crystallite_size_lor,
        strain_lor=strain_lor, strain_gauss=strain_gauss
    )
    p.set_parameters(convolver="displacement",
        zero_error_deg=zero_error_deg, specimen_displacement=specimen_displacement
    )
    if axDiv is not None:
        p.set_parameters(convolver="axial",
            angI_deg=angI_deg, angD_deg=angD_deg,
            axDiv=axDiv, slit_length_source=slit_length_source,
            slit_length_target=slit_length_target,
            length_sample=length_sample,
            n_integral_points=n_integral_points,
        )
    if absorption_coefficient is not None:
        p.set_parameters(convolver="absorption",
            absorption_coefficient=absorption_coefficient,
            sample_thickness=sample_thickness,
        )
    if si_psd_window_bounds is not None:
        p.set_parameters(convolver="si_psd",
         si_psd_window_bounds=si_psd_window_bounds
        )
    p.set_parameters(convolver="receiver_slit",
     slit_width=slit_width
    )
    if tube_tails is not None:
        main_width, tail_left, tail_right, tail_intens=tube_tails
        p.set_parameters(convolver="tube_tails",
            main_width=main_width, tail_left=tail_left,
            tail_right=tail_right, tail_intens=tail_intens
        )

    res=p.compute_line_profile()
    res.add_symbol(profile_class=p)

    return res

if __name__=='__main__' and sys.argv.pop(-1)=='plot':
    
    ## fixed parameters
    mat_wavelengths = numpy.array((1.540596e-10, 1.540596e-10*1.001))
    mat_lor_widths = numpy.array((1.5e-13/5, 1.5e-13/5))
    mat_gauss_widths = numpy.array((1.5e-13/5, 1.5e-13/5))
    mat_intensities = numpy.array((1.0,0.0))
    collect_moment_errors=True
    
    try:
        from scardi_leoni_polydisperse import scardi_leoni_lognormal_polydisperse
        ## @brief subclass with scardi & leoni mixin for testing here;
        ## this shows how to trivially extend the base class with new modules
        class FP_with_poly(FP_profile, scardi_leoni_lognormal_polydisperse):
            pass
        useclass=FP_with_poly
    except:
        useclass=FP_profile

    p=useclass(anglemode="twotheta",
            output_gaussian_smoother_bins_sigma=1.0,
            oversampling=2
        )
    p.debug_cache=False
    #set parameters for each convolver
    p.set_parameters(convolver="emission",
        emiss_wavelengths=mat_wavelengths,
        emiss_intensities=mat_intensities,
        emiss_gauss_widths=mat_gauss_widths,
        emiss_lor_widths=mat_lor_widths,
        crystallite_size_gauss=6000e-9,
        crystallite_size_lor=6000e-9
    )
    p.set_parameters(convolver="axial",
        axDiv="full", slit_length_source=15e-3,
        slit_length_target=10e-3,
        length_sample=12e-3,
        n_integral_points=10
    )
    p.set_parameters(convolver="absorption",
        absorption_coefficient=800.0*100, #like LaB6, in m^(-1)
        sample_thickness=1e-3,
    )
    p.set_parameters(convolver="si_psd",
     si_psd_window_bounds=(0.,2.5e-3)
    )
    p.set_parameters(convolver="receiver_slit",
     slit_width=75e-6,
    )

    p.set_parameters(convolver="tube_tails",
        main_width=200e-6, tail_left=-1e-3,tail_right=1e-3, tail_intens=0.001
    )

    if useclass != FP_profile: #we have scardi & leoni
        from scardi_leoni_polydisperse import H_params_hex
        p.set_parameters(convolver="scardi_leoni_lognormal_polydisperse",
            hkl=(1,0,0), shape_function=H_params_hex,
            crystallite_size=200.e-9, sigma=0.5, flat_is_100=True, 
            crystallite_height_over_edge=2., lattice_c_over_a=1.,
        )

    import time
    t_start=time.time()
    
    profiles=[]
    for ww in (20,):
        for twotheta_x in range(20,161,20):
            #put the compute window in the right place and clear all histories
            p.set_window(twotheta_output_points=5000,
                twotheta_window_center_deg=twotheta_x,
                twotheta_window_fullwidth_deg=ww,
            )
            #set parameters which are shared by many things
            p.set_parameters(twotheta0_deg=twotheta_x,
                dominant_wavelength=mat_wavelengths[0],
                diffractometer_radius=217.5e-3)
        
            for anglescale in (1.,2.,3.):
                p.set_parameters(equatorial_divergence_deg=0.5*anglescale)
                p.set_parameters(convolver="axial",
                    angI_deg=2.5*anglescale, angD_deg=2.5*anglescale,
                )
                profiles.append((twotheta_x, 0.5*anglescale, 2.5*anglescale, p.compute_line_profile()))

            print(p,file=sys.stderr)
    print("execution time=", time.time()-t_start, file=sys.stderr)

    from matplotlib import pyplot as plt
    
    datasets=[]
    for idx, (twotheta, eqdiv, soller, result) in enumerate(profiles):
        #normalize both to integral=sum/dx
        dx1=(result.twotheta_deg[1]-result.twotheta_deg[0])
        peaknorm=40*result.peak/dx1
        plt.plot(result.twotheta_deg, peaknorm, label="%.1f : %.1f" % (eqdiv, soller) )
    
    plt.legend()
    plt.show()
    
    hh=numpy.histogram(moment_list, bins=20)
    print("error histogram=", hh)
    xl=(hh[1][1:]+hh[1][:-1])*0.5 #bin centers
    print("error mean, rms=", (xl*hh[0]).sum()/hh[0].sum(), math.sqrt((xl**2*hh[0]).sum()/hh[0].sum()))
