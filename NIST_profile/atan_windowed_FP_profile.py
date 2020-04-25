# -*- coding: utf-8 -*-
"""
Code to convert FP_Profile to use a windowed spectrum,
and to generate real-space and Fourier IPF
MHM October 2019
"""

from __future__ import print_function

import numpy
import math
import sys

from . import profile_functions_class

class FP_atan_windowed_convolver:

    ## @brief compute the emission spectrum with a window function
    #  @param self self
    #  @return the convolver for the emission
    # requires an array of emiss_wavelengths, emiss_lor_widths, emiss_intensities,
    # and a function emiss_passband_func which takes an argument in wavelength which 
    # will be offset by emiss_window_offset before being applied.
    # @note this tramples the regular conv_emission, so this class must be to the left
    #    of the main class in the final class creation e.g.:
    #    class my_fp(FP_atan_windowed_convolver, profile_functions_class.FP_profile):
    #        pass

    def conv_emission(self):
        """handle Lorentzian-sum emission spectrum with window function"""
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs={}
        kwargs.update(self.param_dicts[me]) #get all of our parameters
        kwargs.update(self.param_dicts["conv_global"])
 
        #convert arrays to lists for key checking
        key={}
        key.update(kwargs)
        for k,v in key.items():
            if hasattr(v,'tolist') and hasattr(v,'__len__') and len(v) > 1:
                key[k]=tuple(v.tolist())
    
        flag, emiss = self.get_conv(me, key, numpy.complex)
        if flag: return emiss #already up to date
        
        xx=type("data", (), kwargs) #make it dot-notation accessible
    
        theta=xx.twotheta0/2
        mono_width = (xx.b_mono/xx.diffractometer_radius) * math.tan(xx.two_theta_mono*math.pi/360)
        dispersion=(2*math.tan(theta)+mono_width)
        relative_passwidth = xx.ibm_source_width/(xx.a_mono*math.tan(xx.two_theta_mono*math.pi/360)) #FW delta-lambda/lambda eq. 6 in paper
        passwidth=dispersion*relative_passwidth/2 #to match Topas macro notation
        #and compute wavelength-scaled passband parameters from this
        center_wavelength=xx.dominant_wavelength #fiducial wavelength from instrument global for mistune, etc.
        center_angle=math.asin(center_wavelength/(2.0 * xx.d)) #nominal bragg angle

        self.abs_window_pm=relative_passwidth*center_wavelength*1e12

        X=self.twothetasamples-2*center_angle #to match notation in Topas macro
        real_emiss=numpy.zeros_like(X) #emission in real space

        for wid, lam0, intens in zip(xx.emiss_lor_widths, xx.emiss_wavelengths, xx.emiss_intensities):
            #wid comes in as an FWHM, convert to HWHM
            xvals=X-dispersion*(lam0-center_wavelength)/center_wavelength
            hw=dispersion*(wid/center_wavelength)/2
            real_emiss += (hw*intens/math.pi)/(xvals*xvals+hw*hw)

        left_passband  = numpy.arctan((X+passwidth-xx.passband_mistune*passwidth)/(xx.passband_shoulder*passwidth))
        right_passband = numpy.arctan((X-passwidth-xx.passband_mistune*passwidth)/(xx.passband_shoulder*passwidth))

        self.raw_emission_shape=numpy.array(real_emiss)
        
        self.passband_window=(left_passband-right_passband)
        real_emiss *= (left_passband-right_passband)
        real_emiss /= real_emiss.sum() #normalize result
        
        self.full_emission_shape=real_emiss
        
        emiss[:]=profile_functions_class.best_rfft(real_emiss)
        
        return emiss

    ## @brief compute the emission spectrum and (for convenience) the particle size widths
    #  @param self self
    #  @return the convolver for the emission and particle sizes
    #  @note the particle size and strain stuff here is just to be consistent with \a Topas
    #  and to be vaguely efficient about the computation, since all of these
    #  have the same general shape.
    def conv_crystallite_size(self):
        """handle crystallite size broadening separately from windowed emission"""
        me=self.get_function_name() #the name of this convolver,as a string
        kwargs={}
        kwargs.update(self.param_dicts["conv_emission"]) #steal our parameters from conv_emission
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
            if hasattr(v,'tolist') and hasattr(v,'__len__') and len(v) > 1:
                key[k]=tuple(v.tolist()) #tuples are better keys
    
        flag, emiss = self.get_conv(me, key, numpy.complex)
        if flag: return emiss #already up to date
        
        xx=type("data", (), kwargs) #make it dot-notation accessible
        
        theta=xx.twotheta0/2
        ## Emission profile FWHM + crystallite broadening (scale factors are Topas choice!) (Lorentzian)
        #note:  the strain broadenings in Topas are expressed in degrees 2theta, must convert to radians(theta) with pi/360
        widths = (
            xx.strain_lor*(math.pi/360) * math.tan(theta) +
            (xx.emiss_wavelengths / (2*xx.crystallite_size_lor * math.cos(theta)))
        )
        #save weighted average width for future reference in periodicity fixer
        self.lor_widths[me]=sum(widths*xx.emiss_intensities)/sum(xx.emiss_intensities)
        #gaussian bits add in quadrature
        gfwhm2s=(
            (xx.strain_gauss*(math.pi/360) * math.tan(theta))**2 +
            (xx.emiss_wavelengths / (xx.crystallite_size_gauss * math.cos(theta)))**2
        )
        
        # note that the Fourier transform of a lorentzian with FWHM 2a
        # is exp(-abs(a omega))
        omega_vals=self.omega_vals
        for wid, gfwhm2, intens in zip(widths, gfwhm2s, xx.emiss_intensities):
            xvals=numpy.clip(omega_vals*(-wid),-100,0)
            sig2=gfwhm2/(8*math.log(2.0)) #convert fwhm**2 to sigma**2
            gxv=numpy.clip((sig2/-2.0)*omega_vals*omega_vals,-100,0)
            emiss += numpy.exp(xvals+gxv)*intens
        return emiss


#cu_ka_spectdata=numpy.array(( #each cluster is wavelength/m, intensity, fwhm/m, from Cu kalpha paper
#        (0.15405925, 3.91, 0.0436e-3), #ka11
#        (0.15410769, 0.474, 0.0558e-3), #ka12
#        (0.15443873, 1.53, 0.0487e-3), #ka21
#        (0.15446782, 0.754, 0.0630e-3), #ka22
#    ))*(1e-9,1,1e-9)
# replaced due to Sphinx problem with scaled values:
cu_ka_spectdata=numpy.array(( #each cluster is wavelength/m, intensity, fwhm/m, from Cu kalpha paper
        (0.15405925e-9, 3.91, 0.0436e-12), #ka11
        (0.15410769e-9, 0.474, 0.0558e-12), #ka12
        (0.15443873e-9, 1.53, 0.0487e-12), #ka21
        (0.15446782e-9, 0.754, 0.0630e-12), #ka22
    ))

if __name__ == "__main__":
    ## fixed parameters

    class FP_windowed(FP_atan_windowed_convolver, profile_functions_class.FP_profile):
        pass

    def fixphase(xfrm, second_try_bin=None):
        """adjust phases to a transform is close to centered. Warning: done in place!"""
        phase=xfrm[1].conjugate()/abs(xfrm[1])
        phase0=phase
        for i in range(1,len(xfrm)):
            xfrm[i]*=phase
            phase *= phase0
        if second_try_bin:
            zz=xfrm[second_try_bin]
            theta=-math.atan2(zz.imag, zz.real)/second_try_bin #better average phase shift
            phase=complex(math.cos(theta), math.sin(theta))
            phase0=phase
            for i in range(1,len(xfrm)):
                xfrm[i]*=phase
                phase *= phase0
    
        return xfrm
    
    from diffraction_constants import omega_length_scale_deg
    
    from matplotlib import pyplot as plt
    plt.style.use('posters_and_pubs')
    
    emiss_wavelengths, emiss_intensities, emiss_lor_widths = cu_ka_spectdata.transpose()
    
    p=FP_windowed(anglemode="twotheta",
            output_gaussian_smoother_bins_sigma=1.0,
            oversampling=2
        )

    #parameters from dbd_ibm_setup_atan_20181116.inc
    #    prm !pixedge  0.08562_0.00399 'in mm
    #    hat  = Rad pixedge/Rs; :  0.0168`_0.00016 num_hats 1 'rectangular pixels
    #
    #    #ifndef REFINE_SOURCE_WIDTH
    #    prm !ibm_source_width  0.04828_0.00032 min 0.0002 max 1
    #    prm !passband_shoulder  0.25772_0.00982 min 0.01 max 1.5  'width ratio of Lorentzian component to flattop
    #    #endif
    #
    #    prm !two_theta_mono 27.27 'degrees twotheta, of the Ge crystal for a peak at CuKa1
    #    prm !b_mono 217.5 'happens to be the same as the DBD radius
    #
    
    #set parameters for each convolver
    p.set_parameters(convolver="emission",
        emiss_wavelengths=emiss_wavelengths,
        emiss_intensities=emiss_intensities,
        emiss_lor_widths=emiss_lor_widths,
        emiss_gauss_widths=[0,0,0,0],
        ibm_source_width=0.05102e-3, #51 um
        passband_shoulder=0.08711,
        two_theta_mono=27.27,
        a_mono=119.e-3,
        b_mono=217.5e-3,
        passband_mistune=-0.14496,
    )
    
    p.set_parameters(convolver="axial",
        axDiv="full", slit_length_source=8e-3,
        slit_length_target=12e-3,
        length_sample=15e-3,
        n_integral_points=20,
        angI_deg=2.73, angD_deg=2.73,
    )
    
    p.set_parameters(convolver="absorption",
        absorption_coefficient=50.0*100*10000, #like SRM 1979, in m^(-1)
    )

    p.set_parameters(convolver="receiver_slit",
     slit_width=0.11562e-3, #convert mm to m
    )

    import time
    t_start=time.time()
    
    profiles=[]
    for twotheta_x in (31.769, 66.37, 135,):
#    for dstar in (4., 7., 11.):
#        twotheta_x=360/math.pi*math.asin(0.15409*dstar/2)
        #put the compute window in the right place and clear all histories
        if twotheta_x > 90: window=2.
        else: window=1.   
        p.set_window(twotheta_output_points=2000,
            twotheta_window_center_deg=twotheta_x,
            twotheta_window_fullwidth_deg=window,
        )
        #set parameters which are shared by many things
        p.set_parameters(twotheta0_deg=twotheta_x,
            dominant_wavelength=emiss_wavelengths[0],
            diffractometer_radius=217.5e-3)
    
        p.set_parameters(equatorial_divergence_deg=1.05)
        
        result=p.compute_line_profile(return_convolver=True)
        print(p,file=sys.stderr)
        abs_window_pm=p.abs_window_pm
        
        plt.figure("real {0:.0f}".format(twotheta_x))
        pk=result.peak/result.peak.max()
        emiss=p.full_emission_shape[::2]/p.full_emission_shape.max()
        raw_emiss=p.raw_emission_shape[::2]/p.raw_emission_shape.max()
        passband=p.passband_window[::2]/p.passband_window.max()
        plt.plot(result.twotheta_deg, pk)
        plt.plot(result.twotheta_deg, emiss)
        plt.plot(result.twotheta_deg, raw_emiss)
        plt.plot(result.twotheta_deg, passband)
        plt.xlabel(r"$2\theta$ / degree")
        plt.savefig("realspace_window_{abs_window_pm:.2f}_tth_{twotheta_x:.0f}.pdf".format(**locals()),
                    bbox_inches='tight', transparent=True, frameon=False)
        #plt.legend()
        #plt.show()
        dstarplot=numpy.array((result.twotheta_deg, 2*numpy.sin(result.twotheta/2)/0.15409, pk)).transpose()
        numpy.savetxt("dstar_{abs_window_pm:.2f}_tth_{twotheta_x:.0f}.dat".format(**locals()),
            dstarplot)
        
        omega=result.omega_inv_deg*omega_length_scale_deg(0.15409, twotheta_x)
        omega_max=numpy.searchsorted(omega, 300)
        fft=numpy.array(result.convolver)
        fft /= fft[0]
        fft[1::2]*=-1 #center it
        fixphase(fft, second_try_bin=5) #adjust to best centering
        numpy.savetxt("fourier_{abs_window_pm:.2f}_tth_{twotheta_x:.0f}.dat".format(**locals()),
            numpy.transpose((omega[:omega_max], numpy.abs(fft[:omega_max])))
        )
        print(fft[:10], file=sys.stderr)
        plt.figure("fourier log {0:.0f}".format(twotheta_x))
        plt.semilogy(omega[:omega_max], numpy.abs(fft[:omega_max]))      
        plt.xlabel(r"$\omega$ / nm")
        plt.savefig("ft_log_window_{abs_window_pm:.2f}_tth_{twotheta_x:.0f}.pdf".format(**locals()),
                    bbox_inches='tight', transparent=True, frameon=False)
        plt.figure("fourier lin {0:.0f}".format(twotheta_x))
        plt.plot(omega[:omega_max], numpy.abs(fft[:omega_max]))       
        plt.xlabel(r"$\omega$ / nm")
        plt.savefig("ft_lin_window_{abs_window_pm:.2f}_tth_{twotheta_x:.0f}.pdf".format(**locals()),
                    bbox_inches='tight', transparent=True, frameon=False)
        #plt.legend()
        #plt.show()
