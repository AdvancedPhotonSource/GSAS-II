# 
# Copyright v1.0, 1.2, 1.3: 2019, John Badger.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

# Version 1.2 is intended to be runnable under Python2 and Python3
# Version 1.3 includes an optional 'glue' term for extended structures
#
# this version modified by R B. Von Dreele for inclusion in GSAS-II

from __future__ import division, print_function
import math
import sys
import os
import copy
import random
import time
import cProfile,pstats
import io as StringIO
import numpy as np
nxs = np.newaxis

def G2shapes(Profile,ProfDict,Limits,data,dlg=None):    

########## FUNCTIONS ########################

    # NEW 1.1 Calculate intensity from P(r) dropping all scaling factors
    
    def ft_to_intensity(aList_q,aList_i_calc,aList_r,aList_pr_model,nbeads):
    
        num_q = len(aList_q)
        num_r = len(aList_r)
    
        count_q = 0
        while count_q < num_q:
    
            q = float(aList_q[count_q])
    
            # Sets P(r=0) =0.0. Later scaling includes a baseline term. 
            integral = 0.0
    
            count_r = 1
            while count_r < num_r:
    
                r = float(aList_r[count_r])
                pr = float(aList_pr_model[count_r])
                qr = q*r
                integral = integral + pr*math.sin(qr)/qr
    
                count_r = count_r + 1
    
            aList_i_calc.append(integral)
    
            count_q = count_q + 1
    
        return;
    
    # NEW 1.1 Scale and Compare I and Ic. Includes a baseline correction 
    
    def score_Ic(aList_q,aList_i,aList_i_sd,aList_i_calc):
    
        num_q = len(aList_q)
        
        idif = 0.0
        isum = 0.0
        sd_sq = 0.0
        chi_sq = 0.0
    
        # Least squares scale for calculated I onto observed I
    
        S = 0.0
        Sx = 0.0
        Sy = 0.0
        Sxx = 0.0
        Sxy = 0.0
    
        count = 0
        while count < num_q:
            x = float(aList_i_calc[count])
            y = float(aList_i[count])
    
            S = S + 1.0
            Sx = Sx + x
            Sy = Sy + y
            Sxx = Sxx + x*x
            Sxy = Sxy + x*y
    
            count = count + 1
    
        delta = S*Sxx - Sx*Sx
        a = (Sxx*Sy - Sx*Sxy)/delta
        b = (S*Sxy - Sx*Sy)/delta
    
        # Apply scale and find statistics
    
        i = 0
        while i < num_q:
            iobs = float(aList_i[i])
            sd = float(aList_i_sd[i])
    
            aList_i_calc[i] = b*float(aList_i_calc[i]) + a        
            icalc = aList_i_calc[i]
            
            idif = idif + abs(iobs - icalc)
            isum = isum + iobs  + icalc
    
            dif = iobs - icalc
            dif_sq = dif*dif
            sd_sq = sd*sd
    
            chi_sq = chi_sq + dif_sq/sd_sq
                
            i = i  + 1
    
        rvalue = 2.0*idif/isum
        chi_sq = chi_sq/num_q
    
        return (chi_sq,rvalue);
    
    # NEW 1.1 Write original and calculated data.
    
    def write_all_data(file_intensity,aList_q,aList_i,aList_i_calc,aString):
    
        num_q = len(aList_q)
    
        file = open(file_intensity,'w',)
        aString = '# ' + aString + '\n'
        file.write(aString)
    
        i = 0
        while i < num_q:
    
            q = aList_q[i]
            intensity = aList_i[i]
            intensity_calc = aList_i_calc[i]
    
            aString = str(q) + ' ' + str(intensity) + ' ' + str(intensity_calc) + '\n'
    
            file.write(aString)
    
            i = i + 1
    
        file.close()
    
        return;    
    
    # NEW 1.1 Read intensity data from GNOM output file
    
    def read_i(aList_q,aList_i,aList_i_sd,inFile,angstrom_scale):
    
        scale_units = 1.0/angstrom_scale
        
        Q,Io,wt,Ic,Ib,Ifb = Profile[:6]
        Qmin = Limits[1][0]
        Qmax = Limits[1][1]
        wtFactor = ProfDict['wtFactor']
        Back,ifBack = data['Back']
        Ibeg = np.searchsorted(Q,Qmin)
        Ifin = np.searchsorted(Q,Qmax)+1    #include last point
        aList_q += list(Q[Ibeg:Ifin]*scale_units)
        aList_i += list(Io[Ibeg:Ifin])
        aList_i_sd += list(1./np.sqrt(wtFactor*wt[Ibeg:Ifin]))
    
#        file = open(inFile,'r')
#        allLines = file.readlines()
#        file.close()
#    
#        parse_data = 'no'
#        for eachLine in allLines:
#    
#            if parse_data == 'yes':
#            
#                aList = eachLine.split()
#                num_params = len(aList)
#    
#                if num_params == 5:
#    
#                    q = float(aList[0]) * scale_units
#                    if q > 0.0:
#                        i = float(aList[1])
#                        i_sd = float(aList[2])
#                        aList_q.append(q)
#                        aList_i.append(i)
#                        aList_i_sd.append(i_sd)
#    
#                else:
#                    if num_params == 0 and len(aList_q) > 0:
#                        parse_data = 'no'
#    
#            if eachLine.find('S          J EXP       ERROR       J REG       I REG') > -1:
#                parse_data = 'yes'
#                    
        return;
    
    # Read PDB for starting point
    
    def read_pdb(aList_beads_x,aList_beads_y,aList_beads_z,pdbfile_in):
    
        xmean = 0.0
        ymean = 0.0
        zmean = 0.0
        nbeads = 0
    
        file = open(pdbfile_in,'r')
        allLines = file.readlines()
        file.close()
    
        for eachLine in allLines:
    
            tag = eachLine[0:6]
            tag = tag.strip()
    
            if tag == 'ATOM':
    
                atom_name = eachLine[13:16]
                atom_name = atom_name.strip()
         
                if atom_name == 'CA':
    
                    x = float(eachLine[30:38])
                    y = float(eachLine[38:46])
                    z = float(eachLine[46:54])
    
                    xmean = xmean + x
                    ymean = ymean + y
                    zmean = zmean + z
    
                    nbeads = nbeads + 1
    
        xmean = xmean/float(nbeads)
        ymean = ymean/float(nbeads)
        zmean = zmean/float(nbeads)
    
        for eachLine in allLines:
    
            tag = eachLine[0:6]
            tag = tag.strip()
    
            if tag == 'ATOM':
    
                atom_name = eachLine[13:16]
                atom_name = atom_name.strip()
         
                if atom_name == 'CA':
    
                    x = float(eachLine[30:38]) - xmean
                    y = float(eachLine[38:46]) - ymean
                    z = float(eachLine[46:54]) - zmean
                    aList_beads_x.append(x)
                    aList_beads_y.append(y)
                    aList_beads_z.append(z)           
    
        return;
    
#    # Write P(r) with SD and calculated value.
#    
#    def pr_writer(aList_pr,aList_r,aList_pr_model,outfile_pr):
#    
#        num_pr = len(aList_pr)
#    
#        file = open(outfile_pr,'w')
#        file.write('#\n')
#    
#        i = 0
#        while i < num_pr:
#        
#            r = float(aList_r[i])
#            pr = float(aList_pr[i])
#            pr_calc = float(aList_pr_model[i])
#            aString = str(r) + ' ' + str(pr) + ' ' + str(pr_calc) + '\n'
#            file.write(aString)
#    
#            i = i + 1
#    
#        file.close()
#    
#        return;
    
    # Write a set of points as a pseudo-PDB file
    
    def pdb_writer(aList_x_write,aList_y_write,aList_z_write,out_file,aString):
        
        atom_number = 0
        res_number = 0
        dummy_b = 0
        num_atoms = len(aList_x_write)
    
        file = open(out_file,'w')
        file.write(aString)
        file.write('\n')
    
        i = 0
        while i < num_atoms:
    
            x = float(aList_x_write[i])
            y = float(aList_y_write[i])
            z = float(aList_z_write[i]) 
            x = '%.3f'%(x)
            y = '%.3f'%(y)
            z = '%.3f'%(z)
            x = x.rjust(8)
            y = y.rjust(8)
            z = z.rjust(8)
            atom_number = atom_number + 1
            res_number = res_number + 1
            atom_number_str = str(atom_number)
            atom_number_str = atom_number_str.rjust(5)
            res_number_str = str(res_number)
            res_number_str = res_number_str.rjust(4)
            bfactor = str(dummy_b) + '.00'
            bfactor = bfactor.rjust(6)
    
            if res_number < 10000:
                aLine_out = 'ATOM  ' + atom_number_str + '  CA  ALA A' + res_number_str + '    ' + x + y + z + '  1.00' + bfactor + '\n'
            else:
                res_number_str = str(res_number - 9999)
                aLine_out = 'ATOM  ' + atom_number_str + '  CA  ALA B' + res_number_str + '    ' + x + y + z + '  1.00' + bfactor + '\n'            
            file.write(aLine_out)
    
            i = i + 1
    
        file.close()
    
        return;
    
    # Evaluate local bead densities and point density on a notional grid
    
    def set_box(aList_beads_x,aList_beads_y,aList_beads_z,\
                aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                aList_box_score,box_step,dmax,rsearch):
    
        dmax_over2 = dmax/2.0
        search_sq = rsearch**2
        num_beads = len(aList_beads_x)
    
        count_x = -dmax_over2
        while count_x < dmax_over2:
            count_y = -dmax_over2
            while count_y < dmax_over2:
                count_z = -dmax_over2
                while count_z < dmax_over2:
    
                    count_local = 0
                    i = 0
                    while i < num_beads:
    
                        dx = float(aList_beads_x[i]) - count_x
                        dy = float(aList_beads_y[i]) - count_y
                        dz = float(aList_beads_z[i]) - count_z
    
                        dsq = dx*dx + dy*dy + dz*dz
    
                        if dsq < search_sq:
                            count_local = count_local + 1
    
                        i = i + 1
    
                    if count_local > 1:
                        aList_box_x_all.append(count_x)
                        aList_box_y_all.append(count_y)
                        aList_box_z_all.append(count_z)
                        aList_box_score.append(count_local)
    
                    count_z = count_z + box_step
                count_y = count_y + box_step
            count_x = count_x + box_step
    
        return;
    
     # Evaluate local bead densities and point density on a notional grid - fast version
    
    def set_box_fast(aList_beads_x,aList_beads_y,aList_beads_z,\
                aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                aList_box_score,box_step,dmax,rsearch):
    
        dmax_over2 = dmax/2.0
        num_box = int(dmax/box_step)
        search_sq = rsearch**2
        
        XYZ = np.meshgrid(np.linspace(-dmax_over2,dmax_over2,num_box),
            np.linspace(-dmax_over2,dmax_over2,num_box),
            np.linspace(-dmax_over2,dmax_over2,num_box))
        XYZ = np.array([XYZ[0].flatten(),XYZ[1].flatten(),XYZ[2].flatten()]).T
        xyz = np.array((aList_beads_y,aList_beads_x,aList_beads_z)).T
        for XYZi in XYZ:
            dsq = np.sum((xyz-XYZi)**2,axis=1)
            count = int(np.sum(np.where(dsq<search_sq,1,0)))
            if count>1:
                aList_box_x_all.append(XYZi[0])
                aList_box_y_all.append(XYZi[1])
                aList_box_z_all.append(XYZi[2])
                aList_box_score.append(count)
        return;
    
   # Establish a volume
    
    def set_vol(aList_box_x_all,aList_box_y_all,aList_box_z_all,aList_box_score,\
                aList_box_x,aList_box_y,aList_box_z,vol_target,box_pt_vol):
    
        num_pts = len(aList_box_score)
        num_tries = int(max(aList_box_score))
        density_thresh = max(aList_box_score) - 1.0
        vol = vol_target + 1.0
        
        i = 0
        while i < num_tries:
    
            density_thresh = density_thresh - 1.0
            num_box_pts = 0.0
    
            j = 0
            while j < num_pts:
                density = float(aList_box_score[j])
                if density >= density_thresh:
                    num_box_pts = num_box_pts + 1.0
                j = j + 1
    
            vol = num_box_pts*box_pt_vol
            if vol < vol_target:
                density_use = density_thresh
     
            i = i + 1
    
        #
    
        num_box_pts1 = 0.0
        num_box_pts2 = 0.0
        density_thresh1 = density_use
        density_thresh2 = density_use - 1.0
        i = 0
        while i < num_pts:
            density_value = float(aList_box_score[i])
            if density_value >= density_thresh1:
                num_box_pts1 = num_box_pts1 + 1.0
            if density_value >= density_thresh2:
                num_box_pts2 = num_box_pts2 + 1.0            
            i = i + 1
            
        vol1 = num_box_pts1*box_pt_vol    
        vol2 = num_box_pts2*box_pt_vol
        delta1 = abs(vol1-vol_target)
        delta2 = abs(vol2-vol_target)
    
        if delta1 < delta2:
            density_thresh = density_thresh1
        else:
            density_thresh = density_thresh2
        
        i = 0
        while i < num_pts:
            
            density_value = float(aList_box_score[i])
            if density_value >= density_thresh:
                aList_box_x.append(aList_box_x_all[i])
                aList_box_y.append(aList_box_y_all[i])
                aList_box_z.append(aList_box_z_all[i])
                    
            i = i + 1
    
        return;
    
    # Find beads that are not in allowed volume
    
    def disallowed_beads(aList_beads_x,aList_beads_y,aList_beads_z,aList_contacts,\
                        aList_box_x,aList_box_y,aList_box_z,rsearch):
    
        num_beads = len(aList_beads_x)
        num_boxes = len(aList_box_x)
        contact_limit_sq = rsearch**2
    
        count = 0
        while count < num_beads:
    
            x_bead = float(aList_beads_x[count])
            y_bead = float(aList_beads_y[count])
            z_bead = float(aList_beads_z[count])
    
            inbox = 'no'
            i = 0
            while i < num_boxes:
    
                x_box = float(aList_box_x[i])
                y_box = float(aList_box_y[i])            
                z_box = float(aList_box_z[i])
                dsq = (x_bead - x_box)**2 + (y_bead - y_box)**2 + (z_bead - z_box)**2
    
                if dsq < contact_limit_sq:
                    inbox = 'yes'
                    i = num_boxes
    
                i = i + 1
    
            if inbox == 'no':
                aList_contacts.append(count)
    
            count = count + 1
    
        return;
    
    # Compute a P(r)
    
    def calc_pr(aList_beads_x,aList_beads_y,aList_beads_z,aList_pr_model,hist_grid):
    
        num_hist = len(aList_pr_model)
        count = 0
        while count < num_hist:
            aList_pr_model[count] = 0.0
            count = count + 1
    
        nbeads = len(aList_beads_x)
        max_dr = (float(num_hist)-1.0)*hist_grid
    
        i = 0
        while i < nbeads:
    
            j = 0
            while j < i:
    
                dr = get_dr(aList_beads_x[i],aList_beads_y[i],aList_beads_z[i],\
                            aList_beads_x[j],aList_beads_y[j],aList_beads_z[j])
                dr = min(dr,max_dr)
    
                # Find pointers and do un-interpolation
                
                dr_grid = dr/hist_grid
                int_dr_grid = int(dr_grid)
    
                int_dr_grid = min(int_dr_grid,num_hist-2)
                ip_low = int_dr_grid
                ip_high = ip_low + 1
                
                ip_high_frac = dr_grid - float(int_dr_grid)
                ip_low_frac = 1.0 - ip_high_frac
    
                aList_pr_model[ip_low] = float(aList_pr_model[ip_low]) + ip_low_frac
                aList_pr_model[ip_high] = float(aList_pr_model[ip_high]) + ip_high_frac
    
                j = j + 1
            i = i + 1    
    
        return;
    
    # Score for rms difference between observed and model histograms
    
    def pr_dif(aList_pr,aList_pr_model,skip):
    
        num_hist = len(aList_pr)
        delta_hist_sum = 0.0
        delta_hist_sum_sq = 0.0
        hist_sum = 0.0
    
        i = skip
        while i < num_hist:
            
            model = float(aList_pr_model[i])
            data = float(aList_pr[i])
            delta_hist = abs(data - model)
            delta_hist_sum = delta_hist_sum + delta_hist
            hist_sum = hist_sum + data
    
            delta_hist_sum_sq = delta_hist_sum_sq + delta_hist*delta_hist
            
            i = i + 1
    
        mean_hist_sum = hist_sum/(num_hist - skip)
        delta_hist_sum_sq = delta_hist_sum_sq/(num_hist - skip)
        delta_hist_sum_sq = math.sqrt(delta_hist_sum_sq)/mean_hist_sum    
    
        return delta_hist_sum_sq;
    
    # Statistics for fractional difference between observed and model histograms
    
    def pr_rfactor(aList_pr,aList_pr_sd,aList_pr_model,skip):
    
        num_hist = len(aList_pr)
        delta_hist_sum = 0.0
        hist_sum = 0.0
    
        i = skip
        while i < num_hist:
            
            model = float(aList_pr_model[i])
            exp = float(aList_pr[i])
            delta_hist = exp - model
            delta_hist_sum = delta_hist_sum + abs(delta_hist)
            hist_sum = hist_sum + exp
            
            i = i + 1
    
        delta_hist_sum = delta_hist_sum/hist_sum
        
        return delta_hist_sum;
    
    # Compute the VDW energy for a interaction
    
    def vdw_energy(econ12,econ6,e_width,dr,bead_sep3):
    
        if dr > bead_sep3:
            vdw = 0.0
        else:
            dr_e6 = dr**6
            dr_e12 = dr_e6**2
            vdw = econ12/dr_e12 - 2.0*econ6/dr_e6
            vdw = max(vdw,e_width)  
    
        return vdw;
    
    def vdw_energies(econ12,econ6,e_width,drs,bead_sep3):

        drs_e6 = drs**6
        drs_e12 = drs_e6**2
        vdws = econ12/drs_e12 - 2.0*econ6/drs_e6
        vdws = np.where(drs>bead_sep3,0.,vdws)
        vdws = np.where(vdws>e_width,vdws,e_width)
        return vdws
        
    # Set a random distribution of beads in a box with maximum extent dmax
    
    def random_beads(aList_beads_x,aList_beads_y,aList_beads_z,\
                     nbeads,dmax,aList_symm,bias_z):
    
        half_side = 0.5*dmax
    
        scale_xy = 1.0 - bias_z
        scale_z = 1.0 + bias_z
        x_range = scale_xy * half_side
        y_range = scale_xy * half_side
        z_range = scale_z * half_side
    
        num_ops = len(aList_symm)
    
        i = 0
        while i < nbeads:
    
            triangle = random.triangular(-0.7,0.7,0.0)
            x = triangle*x_range
            triangle = random.triangular(-0.7,0.7,0.0)
            y = triangle*y_range
            triangle = random.triangular(-0.7,0.7,0.0)
            z = triangle*z_range
                
            aList_beads_x.append(x)
            aList_beads_y.append(y)
            aList_beads_z.append(z)
    
            j = 0
            while j < num_ops:
                aList_s = aList_symm[j]
                m11 = float(aList_s[0])
                m12 = float(aList_s[1])
                m21 = float(aList_s[2])
                m22 = float(aList_s[3])
                
                xs = m11*x + m12*y
                ys = m21*x + m22*y
                zs = z
                aList_beads_x.append(xs)
                aList_beads_y.append(ys)
                aList_beads_z.append(zs)
    
                j = j + 1
    
            i = i + num_symm
    
        return;
    
    # Read experimentalal P(r) from GNOM output file
    
    def read_pr(aList_r,aList_pr,aList_pr_sd,aList_pr_model,\
                aList_pr_model_test,aList_pr_model_test2,inFile):
    
        angstrom_scale = 1.0
        Bins,Dbins,BinMag = data['Pair']['Distribution']
        
        aList_r += list(Bins)
        aList_pr += list(BinMag)
        aList_pr_sd += list(np.ones_like(Bins)/100.)
        aList_pr_model += list(np.zeros_like(Bins))
        aList_pr_model_test += list(np.zeros_like(Bins))
        aList_pr_model_test2 += list(np.zeros_like(Bins))
            
#        file = open(inFile,'r')
#        allLines = file.readlines()
#        file.close()
#    
#        parse_data = 'no'
#        for eachLine in allLines:
#    
#            if parse_data == 'yes':
#            
#                aList = eachLine.split()
#                num_params = len(aList)
#    
#                if num_params == 3:
#                    r = float(aList[0])
#                    pr = float(aList[1])
#                    pr_sd = float(aList[2])
#    
#                    aList_pr.append(pr)
#                    aList_pr_sd.append(pr_sd)
#                    aList_r.append(r)
#                    aList_pr_model.append(0.0)
#                    aList_pr_model_test.append(0.0)
#                    aList_pr_model_test2.append(0.0)
#    
#            if eachLine.find('R          P(R)      ERROR') > -1:
#                parse_data = 'yes'
#                    
        num_hist = len(aList_r)
        hist_grid = float(aList_r[1]) - float(aList_r[0])
    
    
        # Heuristic for checking units
#        test_r = max(aList_r)
#        if test_r < 30.0:
#    
#            aString = 'P(r)appears to be in nm. Converting to Angstrom units'
#            print (aString)
#            angstrom_scale = 10.0
#    
#            i = 0
#            while i < num_hist:
#                aList_r[i] = angstrom_scale * aList_r[i]
#                i = i + 1
#    
#            hist_grid = angstrom_scale * hist_grid
#    
#        i = 0
#        while i < num_hist:
#            r = float(aList_r[i])
#            r_calc = float(i)*hist_grid
#    
#            if abs(r - r_calc) > 0.15:
#                aString = 'Input P(r) grid is irregular! Exiting'
#                print (aString)
#                time.sleep(4)
#                sys.exit(1)
#    
#            i = i + 1
#    
        dmax = aList_r[num_hist-1]
    
        # Pad histogram by 5 Angstrom
    
        ipad = int(5.0/hist_grid)
        
        i = 0
        while i < ipad:
            r = dmax + float(i)*hist_grid
            aList_pr.append(0.0)
            aList_pr_sd.append(0.0)
            aList_r.append(r)
            aList_pr_model.append(0.0)
            aList_pr_model_test.append(0.0)
            aList_pr_model_test2.append(0.0)
            i = i + 1
    
        return (dmax,hist_grid,num_hist,angstrom_scale);
    
    # Scale P(r) onto model P(r) assuming same grid
    
    def scale_pr(aList_pr,aList_pr_sd,aList_pr_model):
    
        num_hist = len(aList_pr)
        total_dr = 0.0
        total_pr = 0.0
    
        i = 0
        while i < num_hist:
            total_dr = total_dr + float(aList_pr_model[i])
            total_pr = total_pr + float(aList_pr[i])
            i = i + 1
    
        scale = total_dr/total_pr
        
        i = 0
        while i < num_hist:
            aList_pr[i] = scale*float(aList_pr[i]) 
            aList_pr_sd[i] = scale * float(aList_pr_sd[i])
            i = i + 1
    
        return;
    
    # Return a non-zero distance between two coordinates
    
    def get_dr(x1,y1,z1,x2,y2,z2):
    
        x1 = float(x1)
        y1 = float(y1)
        z1 = float(z1)
        x2 = float(x2)
        y2 = float(y2)
        z2 = float(z2)    
        dr = (x1 - x2)**2 + (y1-y2)**2 + (z1-z2)**2
        dr = max(dr,0.1)
        dr = math.sqrt(dr)
    
        return dr;
    
    # Return non-zero distances one coordinate & the rest
    
    def get_drs(xyz,XYZ):
        
        return  np.sqrt(np.sum((XYZ-xyz)**2,axis=1))
    
    # Find center of beads within a radii
    
    def center_beads(x,y,z,aList_beads_x,aList_beads_y,aList_beads_z,radii_1,radii_2):
    
        num_beads = len(aList_beads_x)
    
#        xsum = 0.0
#        ysum = 0.0
#        zsum = 0.0
#        count_beads = 0.0
#    
#        i = 0
#        while i < num_beads:
#            
#            dr = get_dr(x,y,z,aList_beads_x[i],aList_beads_y[i],aList_beads_z[i])
#    
#            if dr > radii_1 and dr < radii_2:
#                count_beads = count_beads + 1.0
#                xsum = xsum + float(aList_beads_x[i])
#                ysum = ysum + float(aList_beads_y[i])
#                zsum = zsum + float(aList_beads_z[i])           
#    
#            i = i + 1
#    
#        if count_beads > 0.0:
#            xsum = xsum/count_beads
#            ysum = ysum/count_beads
#            zsum = zsum/count_beads
#            delta = (xsum - x)**2 + (ysum - y)**2 + (zsum - z)**2
#            delta = math.sqrt(delta)
#        else:
#            delta = 0.0
            
        XYZ = np.array([aList_beads_x,aList_beads_y,aList_beads_z]).T
        xyz = np.array([x,y,z])
        drs = get_drs(xyz,XYZ)
        sumXYZ = np.array([XYZ[i] for i in range(num_beads) if radii_1<drs[i]<radii_2])
        count_beads = sumXYZ.shape[0]
        
        delta = 0.0
        if count_beads:
            delta = np.sqrt(np.sum(((np.sum(sumXYZ,axis=0)/count_beads)-xyz)**2))
            
        return delta;
    
    # Obtain mean of total VDW energy 
    
    def get_total_energy(aList_beads_x,aList_beads_y,aList_beads_z,econ12,econ6,bead_sep3):
    
        nbeads = len(aList_beads_x)
        vdw_all = 0.0
        
        
        i = 0
        while i < nbeads:
            xyz = np.array([aList_beads_x[i],aList_beads_y[i],aList_beads_z[i]])
            XYZ = np.array([aList_beads_x[:i],aList_beads_y[:i],aList_beads_z[:i]]).T
            drs = get_drs(xyz,XYZ)
            vdws = vdw_energies(econ12,econ6,e_width,drs,bead_sep3)
            vdw_all += np.sum(vdws)
            i += 1
        
#        i = 0
#        while i < nbeads:
#            j = 0
#            while j < i:
#                dr = get_dr(aList_beads_x[i],aList_beads_y[i],aList_beads_z[i],\
#                            aList_beads_x[j],aList_beads_y[j],aList_beads_z[j])
#                vdw = vdw_energy(econ12,econ6,e_width,dr,bead_sep3)            
#                vdw_all = vdw_all + vdw
#                j = j + 1
#            i = i + 1
            
        vdw_all = vdw_all/float(nbeads)
    
        return vdw_all;
    
    # Energy minimize
    
    def e_min(aList_beads_x,aList_beads_y,aList_beads_z,bead_sep,bead_sep3,aList_symm):
    
        eps = bead_sep/(2.0**(1.0/6.0))
        eps12 = eps**12
        eps6 = eps**6
        step_max = bead_sep
        scale = 0.0
        icount = -1
    
        nbeads = len(aList_beads_x)
        num_ops = len(aList_symm)
        num_cycles = 51
    
        i = 0
        while i < num_cycles:
    
            icount = icount + 1
    
            aList_beads_x_new = []
            aList_beads_y_new = []
            aList_beads_z_new = []
    
            sum_forces_scale = 0.0
    
            k = 0
            while k < nbeads - num_ops:
     
                xold = float(aList_beads_x[k])
                yold = float(aList_beads_y[k])
                zold = float(aList_beads_z[k])
    
                
#                fxyz = np.zeros(3)
#                XYZ = np.array([aList_beads_x,aList_beads_y,aList_beads_z]).T
#                xyz = np.array([xold,yold,zold])
#                drs = get_drs(xyz,XYZ)
#                drs = np.where(drs>eps,drs,eps)
#                drs_sq = drs*drs
#                drs12 = drs_sq**6
#                drs6 = drs_sq**3
#                dxs = xyz-XYZ
#                forces = np.where(drs<bead_sep3,(1.0/drs_sq)*(eps12/drs12 - 0.5*eps6/drs6),0.0)
#                fxyz = np.sum(forces[:,nxs]*dxs,axis=0)
#                sum_forces_scale = np.sum(np.abs(fxyz))
#                
#                xstep = scale*fxyz[0]
#                ystep = scale*fxyz[1]
#                zstep = scale*fxyz[2]
                
                
                
                fx = 0.0
                fy = 0.0
                fz = 0.0
                
                j = 0
                while j < nbeads:
    
                    xj = aList_beads_x[j]
                    yj = aList_beads_y[j]
                    zj = aList_beads_z[j]
                    
                    dr = get_dr(xold,yold,zold,xj,yj,zj)
    
                    # Truncate very steep
                    dr = min(eps,dr)
    
                    if dr < bead_sep3:
                        dr_sq = dr*dr
                        dr12 = dr_sq**6
                        dr6 = dr_sq**3
    
                        dx = xold - xj
                        dy = yold - yj
                        dz = zold - zj
    
                        force = (1.0/dr_sq)*(eps12/dr12 - 0.5*eps6/dr6)
                        fx = fx + force*dx
                        fy = fy + force*dy
                        fz = fz + force*dz
    
                        sum_forces_scale = sum_forces_scale + abs(fx) + abs(fy) + abs(fz)
    
                    j = j + 1
    
                #
                xstep = scale*fx
                ystep = scale*fy
                zstep = scale*fz
                
                if xstep > 0.0:
                    xstep = min(xstep,step_max)
                else:
                    xstep = max(xstep,-step_max)
                    
                if ystep > 0.0:
                    ystep = min(ystep,step_max)
                else:
                    ystep = max(ystep,-step_max)
                    
                if zstep > 0.0:
                    zstep = min(zstep,step_max)
                else:
                    zstep = max(zstep,-step_max)
                
                xtest = xold + xstep
                ytest = yold + ystep
                ztest = zold + zstep
                aList_beads_x_new.append(xtest)
                aList_beads_y_new.append(ytest)
                aList_beads_z_new.append(ztest)
    
                # Apply shifs to symm positions
                l = 0
                while l < num_ops:
                    aList_s = aList_symm[l]
                    m11 = float(aList_s[0])
                    m12 = float(aList_s[1])
                    m21 = float(aList_s[2])
                    m22 = float(aList_s[3])
                
                    xs = m11*xtest + m12*ytest
                    ys = m21*xtest + m22*ytest
                    zs = ztest
    
                    aList_beads_x_new.append(xs)
                    aList_beads_y_new.append(ys)
                    aList_beads_z_new.append(zs)
    
                    l = l + 1
                
                #
    
                k = k + num_ops + 1
    
            # Apply shifted positions after first cycle
            if i > 0:
                
                m = 0
                while m < nbeads:
                    aList_beads_x[m] =  aList_beads_x_new[m] 
                    aList_beads_y[m] =  aList_beads_y_new[m]
                    aList_beads_z[m] =  aList_beads_z_new[m]
                    m = m + 1
    
            #
    
            mean_force = (num_ops+1)*sum_forces_scale/(nbeads*3.0)
            scale = bead_sep/mean_force
    
            vdw_all = get_total_energy(aList_beads_x,aList_beads_y,aList_beads_z,econ12,econ6,bead_sep3)
    
            if icount == 0:
                aString = 'Emin cycle: ' + str(i) + ' Energy: ' + str('%4.2f'%(vdw_all))
                print (aString)
                icount = -10
    
            if vdw_all < 0.0:
                i = num_cycles
    
            i = i + 1
    
        return;
    
    # Set up symmetry operators for rotational symmetry
    
    def make_symm(aList_symm,num_symm):
    
        angle_step = 360.0/float(num_symm)
    
        i = 1
        while i < num_symm:
            theta = float(i) * angle_step
            theta = math.radians(theta)
            cos_theta = math.cos(theta)
            sin_theta = math.sin(theta)
            aList_s = [cos_theta,sin_theta,-sin_theta,cos_theta]
            aList_symm.append(aList_s)
            i = i + 1
    
        return aList_symm;
    
    # Set up a shift vector in P(r) for a change in bead position
    
    def pr_shift_atom(aList_pr_model_test2,x1,y1,z1,\
                      aList_beads_x,aList_beads_y,aList_beads_z,hist_grid,ii):
    
        num_hist = len(aList_r)
        max_dr = (float(num_hist)-1.0)*hist_grid
#        num_beads = len(aList_beads_x)
        
#        aList_pr_model_test2 = num_hist*[0.0,]
#        
#        i = 0
#        while i < num_hist:
#            aList_pr_model_test2[i] = 0.0
#            i = i + 1
            
        XYZ = np.array([aList_beads_x,aList_beads_y,aList_beads_z]).T
        xyz = np.array([x1,y1,z1])
        drs = get_drs(xyz,XYZ)
        drs_grid = np.where(drs>max_dr,max_dr,drs)/hist_grid
        int_drs_grid = np.array(drs_grid,dtype=np.int)
        int_drs_grid = np.where(int_drs_grid>num_hist-2,num_hist-2,int_drs_grid)                
        ip_lows = int_drs_grid
        ip_highs = ip_lows + 1            
        ip_high_fracs = drs_grid - int_drs_grid
        ip_low_fracs = 1.0 - ip_high_fracs
        for ip_low in ip_lows:
            aList_pr_model_test2[ip_low] += ip_low_fracs[ip_low]
        for ip_high in ip_highs:
            aList_pr_model_test2[ip_high] += ip_high_fracs[ip_high]
    
    
#        i = 0
#        while i < num_beads:
#    
#            if i != ii:
#                x2 = float(aList_beads_x[i])
#                y2 = float(aList_beads_y[i])
#                z2 = float(aList_beads_z[i])
#                dr = get_dr(x1,y1,z1,x2,y2,z2)
#                dr = min(dr,max_dr)
#                dr_grid = dr/hist_grid
#                int_dr_grid = int(dr_grid)
#                int_dr_grid = max(int_dr_grid,0)
#                int_dr_grid = min(int_dr_grid,num_hist-2)                
#                ip_low = int_dr_grid
#                ip_high = ip_low + 1            
#                ip_high_frac = dr_grid - float(int_dr_grid)
#                ip_low_frac = 1.0 - ip_high_frac
#                aList_pr_model_test2[ip_low] = float(aList_pr_model_test2[ip_low]) + ip_low_frac
#                aList_pr_model_test2[ip_high] = float(aList_pr_model_test2[ip_high]) + ip_high_frac
#    
#            i = i + 1
#    
        return;
    
    # Recenter set of beads to origin
    
    def recenter_pdb(aList_beads_x,aList_beads_y,aList_beads_z):
    
        nbeads = len(aList_beads_x)
        xsum = 0.0
        ysum = 0.0
        zsum = 0.0
    
        i = 0
        while i < nbeads:
            xsum = xsum + float(aList_beads_x[i])
            ysum = ysum + float(aList_beads_y[i])
            zsum = zsum + float(aList_beads_z[i])
            i = i + 1
    
        xmean = xsum/float(nbeads)
        ymean = ysum/float(nbeads)
        zmean = zsum/float(nbeads)
    
        i = 0
        while i < nbeads:
            aList_beads_x[i] = float(aList_beads_x[i]) - xmean
            aList_beads_y[i] = float(aList_beads_y[i]) - ymean
            aList_beads_z[i] = float(aList_beads_z[i]) - zmean
            i = i + 1
    
        return;
    
    #############
    # EXECUTION #
    #############
    
    #profiling start
    pr = cProfile.Profile()
    pr.enable()
    time0 = time.time()
    
    version_aString = 'Program: SHAPES version 1.3'
    
    print (version_aString)
    aString = 'Author: John Badger'
    print (aString)
    aString = 'Copyright: 2019, John Badger'
    print (aString)
    aString = 'License: GNU GPLv3'
    print (aString)
    
    localtime = time.asctime( time.localtime(time.time()) )
    aString = 'Starting time: ' + str(localtime) + '\n'
    print (aString)
    
#    aList_summary = []
#    aList_summary.append(version_aString)
#    aList_summary.append(str(localtime))
    
    ######################
    # Start up parmeters #
    ######################
#        data['Shapes'] = {'outName':'','NumAA':100,'Niter':1,'AAscale':1.0,'Symm':1,'bias-z':0.0,
#            'inflateV':1.0,'AAglue':0.0}
        
    nbeads = 0
    num_sols = 1
    num_aa = 1.0
    num_symm = 1
    bias_z = 0.0
    inflate = 1.0
    prefix = ''
    surface_scale = 0.0
    starting_pdb = 'no'
    inFile = 'none'
    pdbfile_in  = 'none'
    shapeDict = data['Shapes']
    prefix = shapeDict['outName']
    nbeads = shapeDict['NumAA']
    num_sols = shapeDict['Niter']
    num_aa = shapeDict['AAscale']
    num_symm = shapeDict['Symm']
    bias_z = shapeDict['bias-z']
    inflate = shapeDict['inflateV']
    surface_scale = shapeDict['AAglue']
    pdbOut = shapeDict['pdbOut']
    box_step = shapeDict.get('boxStep',4.0)
    Phases = []
    Patterns = []
    PRcalc = []
    
#    # Parse
#    
#    if os.path.exists('shapes_ip.txt'):
#        file = open('shapes_ip.txt','r')
#        allLines = file.readlines()
#        file.close()
#    else:
#        aString = 'The local parameter file shapes_ip.txt was not found ! Exiting'
#        print (aString)
#        time.sleep(4)
#        sys.exit(1)
#    
#    for eachLine in allLines:
#    
#        aList = eachLine.split()
#        
#        num_params = len(aList)
#        if num_params > 0:
#    
#            if aList[0] == 'num_amino_acids':
#                nbeads = int(aList[1])
##            elif aList[0] == 'input_pr':
##                inFile = str(aList[1])
#            elif aList[0] == 'num_solns':
#                num_sols = int(aList[1])
#            elif aList[0] == 'num_aa_scale':
#                num_aa = float(aList[1])
#            elif aList[0] == 'symm':
#                num_symm = int(aList[1])
#            elif aList[0] == 'bias_z':
#                bias_z = float(aList[1])
#            elif aList[0] == 'inflate_vol':
#                inflate = float(aList[1])
#            elif aList[0] == 'pdb_start':
#                pdbfile_in = str(aList[1])
#            elif aList[0] == 'id':
#                prefix = str(aList[1]) + '_'
#            elif aList[0] == 'glue':
#                surface_scale = float(aList[1])
            
    
    # Check inputs
    
    if num_sols > 0:
        aString = 'Number of runs: ' + str(num_sols)
        print (aString)
    else:
        aString = 'Zero reconstruction runs specified! Exiting'
        print (aString)
        time.sleep(4)
        sys.exit(1)
    
    #
    if nbeads == 0:
        if os.path.exists(pdbfile_in):
            aString = 'Will use CA atoms from ' + str(pdbfile_in) + ' as the initial bead distribution.'
            print (aString)
            starting_pdb = 'yes'
        else:
            aString = 'Zero amino acid count specified and no starting file found. Exiting'
            print (aString)
            time.sleep(4)
            sys.exit(1)
    else:
        aString = 'Number of amino acids: ' + str(nbeads)
        print (aString)
    
    #
#    if os.path.exists(inFile):
#        aString = 'Input P(r) file name: ' + str(inFile)
#        print (aString)
#    else:
#        aString = 'P(r) input file not found. Exiting'
#        print (aString)
#        time.sleep(4)
#        sys.exit(1)
    
    #
    if num_aa == 0.0:
        aString = 'Scale for amino acid count to particle number cannot be zero! Exiting'
        print (aString)
        time.sleep(4)
        sys.exit(1)
    else:
        aString = 'Scale aa to bead count: ' + str(num_aa)
        print (aString)
        
    #
    if num_symm == 0:
        aString = 'Rotational symmetry cannot be zero! Set to 1 for no symmetry. Exiting'
        print (aString)
        time.sleep(4)
        sys.exit(1)
    else:
        aString = 'Point symmetry: ' + str(num_symm)
        print (aString)
    
    #
    if bias_z > 0.2:
        aString = 'Max bias on Z axis for initial particle distribution is 0.2 (rods). Reset to 0.2.'
        print (aString)
        bias_z = 0.2
    elif bias_z < -0.2:
        aString = 'Min bias on Z axis for initial particle distribution is -0.2 (disks). Reset to -0.2.'
        print (aString)
        bias_z = -0.2
    else:
        aString = 'Z-axis bias: ' + str(bias_z)
        print (aString)
    
    #
    if inflate < 0.0:
        aString = 'Inflation of PSV cannot be less than zero! Exiting'
        print (aString)
        time.sleep(4)
        sys.exit(1)
    elif inflate > 2.0:
        aString = 'Inflation of PSV cannt be greater than 2.0! Exiting'
        print (aString)
        time.sleep(4)
        sys.exit(1)    
    else:
        aString = 'PSV inflation factor: ' + str(inflate)
        print (aString)
    
    #
    if surface_scale > 0.0:
        aString = 'Cavity weight: ' + str(surface_scale)
        print (aString)
    
    ########## UNIVERSAL CONSTANTS ######################
    
    # No of macrocycles (gives extra cycles at constant volume after num_contract)
    niter = 160
    
    # No of contraction cycles
    num_contract = 140
    
    # Number of cycles at each fixed volume 
    num_micro_cyc = 10
    
    # Final quench
    num_sa_max = niter - num_micro_cyc
    
    # Initial scale for P(r) shifts versus E shifts
    hscale = 3000.0
    
    # Standard deviation for annealing acceptance (cf well-depth of -1 unit for two beads)
    sd_mc = float(num_symm) * 2.0
    
    # Fiddle factor for keeping the accessible, molecular volume larger than PSV 
    scale_vol = 1.15
    
    # Standard amino acid volume MW = 110.0 x 1.21 i.e. mean mw x mw-to-vol-scale
    vol_bead = 133.1 
    
    # Bead separation for best packing ~5.6 (I think)
    #- 75% better than rectangular grid 5.1 for this amino acid vol
    bead_sep = 5.6 
    
    # Usually num_aa is unity. Adjust parameters otherwise
    if num_aa != 1 and nbeads != 0:
        nbeads = int(nbeads*num_aa)
        vol_bead = vol_bead / num_aa
        bead_sep = (vol_bead * 4/3)**(1.0/3.0)
    
    # Increase bead separation for inflated volumes
    bead_sep = bead_sep * inflate**(1.0/3.0)
    
    # Partial specific volumes at start and end
    
    if starting_pdb == 'yes':
        nmols_vol_start = 1.1 * inflate
    else:
        nmols_vol_start = 2.0 * inflate
    
    nmols_vol_end = 1.0 * inflate
    nmols_vol_subtract = nmols_vol_start - nmols_vol_end
    
    # Box parametere
#    box_step = 4.0      #5.0
    box_pt_vol = box_step*box_step*box_step
    
    # Energy parameters - flat bottomed VDW (2.0A for a 5.6A bead separation)
    
    well_width = 0.36*bead_sep
    econ12 = bead_sep**12
    econ6 = bead_sep**6
    r_width = bead_sep + well_width
    r_width6 = r_width**6
    r_width12 = r_width6**2
    e_width = econ12/r_width12 - 2.0*econ6/r_width6
    bead_sep3 = 3.0*bead_sep
    abs_e_width = abs(e_width)
    
    # Range for box identification (might need to increase for poor data)
    rsearch = (bead_sep + 0.5*well_width)*1.5
    
    # Search range for optional cavity inhibition energy term
    radii_1 = 1.5*bead_sep 
    radii_2 = 4.0*bead_sep
    
    # Setup symmetry operators
    
    aList_symm = []
    aList_symm = make_symm(aList_symm,num_symm)
    num_ops = len(aList_symm)
    
    # Read experimental histogram
    
    aList_r = []
    aList_pr = []
    aList_pr_sd = []
    aList_pr_model = []
    aList_pr_model_test = []
    aList_pr_model_test2 = []
    
    (dmax,hist_grid,num_hist_in,angstrom_scale) = read_pr(aList_r,aList_pr,aList_pr_sd,\
                                            aList_pr_model,aList_pr_model_test,\
                                            aList_pr_model_test2,inFile)
    
#    dmax_over2 = dmax/2.0       
    num_hist = len(aList_r)
    
    aString = 'Number of points read from P(r): ' + str(num_hist_in)
    print (aString)
    aString = 'Grid sampling: ' + str(hist_grid) + ' Dmax: ' + str(dmax)
    print (aString)
    
    # Skip over initial points in scoring
    
    skip = r_width/float(num_hist)
    skip = int(skip) + 2
    
    # Read intensity data that was used for P(r)
    
    aList_q = []
    aList_i = []
    aList_i_sd = []
    
    read_i(aList_q,aList_i,aList_i_sd,inFile,angstrom_scale)
    
    num_intensities = len(aList_q)
    aString = 'Number of intensity data points read: ' + str(num_intensities)
    print (aString)
    
    #########################
    # CYCLE OVER SOLUTIONS  #
    #########################
    
    i_soln = 0
    while i_soln < num_sols:
    
        file_no = str(i_soln + 1)
    
        aString = '\nReconstruction trial: ' + str(file_no)
        print (aString)
                         
        aString  = 'Trial:' + file_no
#        aList_summary.append(aString)
    
        file_beads = prefix + 'beads_' + file_no + '.pdb'
#        file_pr = prefix + 'pr_calc_' + file_no + '.dat'
        file_psv = prefix + 'psv_shape_' + file_no + '.pdb'
#        file_intensity = prefix + 'intensity_' + file_no + '.dat'
    
        # Setup initial bead distribution
    
        aList_beads_x = []
        aList_beads_y = []
        aList_beads_z = []
    
        # Re-initialize standard deviation for annealing acceptance
        sd_mc = float(num_symm) * 2.0
    
        # Set random bead distribution
    
        if starting_pdb == 'yes':
            read_pdb(aList_beads_x,aList_beads_y,aList_beads_z,pdbfile_in)
            nbeads = len(aList_beads_x)
            num_symm = 1
            aString = 'Number of CA sites read: ' + str(nbeads)
            print (aString)
            aString = 'Symmetry set to 1 (required)'
            print (aString)
            aString = 'Input center was shifted to the origin'
            print (aString)
        else:
            random_beads(aList_beads_x,aList_beads_y,aList_beads_z,nbeads,dmax,aList_symm,bias_z)
            nbeads = len(aList_beads_x)
            aString = 'Number of beads randomly placed: ' + str(nbeads)
            print (aString)
    
        # Protein partial specific volume
        psv_vol = float(nbeads)*vol_bead
    
        # Histogram of inter-bead distance 
    
        calc_pr(aList_beads_x,aList_beads_y,aList_beads_z,aList_pr_model,hist_grid)
    
        # Scale experimental P(r) and model histogram  
    
        scale_pr(aList_pr,aList_pr_sd,aList_pr_model)
    
        # Minimize initial energy using expanded VDW
    
        if starting_pdb != 'yes':
            aString = 'Minimize energy of initial positions'
            print (aString)
            bead_sep_e = 1.35*bead_sep
            bead_sep3_e = 3.0*bead_sep_e
            e_min(aList_beads_x,aList_beads_y,aList_beads_z,bead_sep_e,bead_sep3_e,aList_symm)
        else:
            aString = 'Skipping energy minimization of initial positions'
            print (aString)
    
        # Get the initial score between observed and calculated P(r)
    
        hist_score_best = pr_dif(aList_pr,aList_pr_model,skip)
        aString = 'Initial rms P(r): ' + str('%4.3f'%(hist_score_best))
        print (aString)
        
        aList_i_calc = []
        ft_to_intensity(aList_q,aList_i_calc,aList_r,aList_pr_model,nbeads)
        (chi_sq,rvalue) = score_Ic(aList_q,aList_i,aList_i_sd,aList_i_calc)
        aString = 'Initial Rvalue: ' + str('%4.3f'%(rvalue)) + ' CHI-squared: ' + str('%4.3f'%(chi_sq))
        print (aString)
    
        ###########################
        # Iterate                 #
        ###########################
    
        num_boxes = 0
        count_boxing = 0
        fraction_psv = 0
        success_rate_all = 0.0
        box_iter = num_micro_cyc - 1
    
        sum_delta_pack = 0.0
        
        count_it = 0
        while count_it < niter:
    
            success = 0
            count_hist_yes = 0
            sum_e = 0.0
            sum_h = 0.0
    
            # Find populated volumes and fix solution every 10 macrocycles
    
            box_iter = box_iter + 1
        
            if box_iter == num_micro_cyc:
    
                box_iter = 0
                count_boxing = count_boxing + 1
    
                if count_it < num_contract - 1:
                    scale = float(count_it)/float(num_contract)            
                else:
                    scale = 1.0
    
                # Establish confinement volume using a local average
    
                aList_box_x_all = []
                aList_box_y_all = []
                aList_box_z_all = []
                aList_box_score = []
    
                recenter_pdb(aList_beads_x,aList_beads_y,aList_beads_z)
    
                # Adaptive masking was not helpful
                # rsearch_use = (2.0 - scale)*rsearch
    
                set_box_fast(aList_beads_x,aList_beads_y,aList_beads_z,\
                        aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                        aList_box_score,box_step,dmax,rsearch)
    
                aList_box_x = []
                aList_box_y = []
                aList_box_z = []
    
                psv_ratio = nmols_vol_start - scale*nmols_vol_subtract
                vol_target = scale_vol*(psv_ratio*psv_vol)
    
                set_vol(aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                        aList_box_score,aList_box_x,aList_box_y,aList_box_z,\
                        vol_target,box_pt_vol)
    
                num_boxes = len(aList_box_x)
                fraction_psv = float(num_boxes)*box_pt_vol/psv_vol
    
                # Find beads that are ouside the allowed volume
    
                aList_contacts = []
                disallowed_beads(aList_beads_x,aList_beads_y,aList_beads_z,aList_contacts,\
                        aList_box_x,aList_box_y,aList_box_z,rsearch)
                num_outof_box = len(aList_contacts)
    
                aString = 'Target volume: ' + str('%4.2f'%(scale_vol*psv_ratio)) + ' Actual volume: ' + \
                          str('%4.2f'%(fraction_psv)) + ' Beads outside volume: ' + str(num_outof_box)
                print (aString)
    
                # Recalculate P(r) and rescore for reliability
    
                calc_pr(aList_beads_x,aList_beads_y,aList_beads_z,aList_pr_model,hist_grid)
                hist_score_best = pr_dif(aList_pr,aList_pr_model,skip)
                
#                aList_i_calc = []
#                ft_to_intensity(aList_q,aList_i_calc,aList_r,aList_pr_model,nbeads)
#                (chi_sq,rvalue) = score_Ic(aList_q,aList_i,aList_i_sd,aList_i_calc)
    
                # Reset SA deviation if mean success rate over last trials is under 0.1
    
                mean_success_rate = float(success_rate_all)/float(num_micro_cyc)
    
                if count_it < num_contract and count_boxing != 1:
    
                    if mean_success_rate < 0.1:            
                        sd_mc = 1.3*sd_mc
                        aString = 'Raising allowed energy deviation to ' + str('%4.2f'%(sd_mc))
                        print (aString)
    
                    if mean_success_rate > 0.2:
                        sd_mc = 0.7*sd_mc
                        aString = 'Reducing allowed energy deviation to ' + str('%4.2f'%(sd_mc))
                        print (aString)                
    
                success_rate_all = 0.0
                
            # Make one macrocycle that is a run over the nbeads 
    
            ii = 0
            while ii < nbeads:
    
                # Initialize
    
                energy_old = 0.0
                energy_new = 0.0
    
                i = 0
                while i < num_hist:
                    aList_pr_model_test[i]  = 0.0
                    i = i + 1
    
                # Select a target bead and make trial shift
    
                xold = float(aList_beads_x[ii])
                yold = float(aList_beads_y[ii])
                zold = float(aList_beads_z[ii])
    
                ibox = random.randint(0,num_boxes-1)
                xtest = float(aList_box_x[ibox]) + random.uniform(-rsearch,rsearch)
                ytest = float(aList_box_y[ibox]) + random.uniform(-rsearch,rsearch)      
                ztest = float(aList_box_z[ibox]) + random.uniform(-rsearch,rsearch)
    
                # Calculate and capture symmetry mates
    
                aList_temp_save_x = []
                aList_temp_save_y = []
                aList_temp_save_z = []
                aList_symm_x = []
                aList_symm_y = []
                aList_symm_z = []
            
                l = 0
                while l < num_ops:
                    aList_s = aList_symm[l]
                    m11 = float(aList_s[0])
                    m12 = float(aList_s[1])
                    m21 = float(aList_s[2])
                    m22 = float(aList_s[3])
                
                    xs = m11*xtest + m12*ytest
                    ys = m21*xtest + m22*ytest
                    zs = ztest
    
                    aList_symm_x.append(xs)
                    aList_symm_y.append(ys)
                    aList_symm_z.append(zs)
    
                    ipt = ii + l + 1
                    aList_temp_save_x.append(aList_beads_x[ipt])
                    aList_temp_save_y.append(aList_beads_y[ipt])
                    aList_temp_save_z.append(aList_beads_z[ipt])
    
                    l = l + 1
    
                # Get initial VDW energy for interactions of this bead with all others
    
                i = 0
                while i < nbeads:
                    if i != ii:
                        x = float(aList_beads_x[i])
                        y = float(aList_beads_y[i])
                        z = float(aList_beads_z[i])
                        dr_old = get_dr(xold,yold,zold,x,y,z)
                        vdw_old = vdw_energy(econ12,econ6,e_width,dr_old,bead_sep3)
                        energy_old = energy_old + num_symm*vdw_old
                    i = i + 1
    
                # Get initial contributions to P(r)  
    
                aList_pr_model_test2 = num_hist*[0.0,]
                pr_shift_atom(aList_pr_model_test2,xold,yold,zold,aList_beads_x,\
                              aList_beads_y,aList_beads_z,hist_grid,ii)
    
                i = 0
                while i < num_hist:
                    aList_pr_model_test[i] = aList_pr_model_test2[i]
                    i = i + 1
    
                # Get VDW energy for interactions of the shifted bead with all others
    
                l = 0
                while l < num_ops:       
                    ipt = ii + l + 1
                    aList_beads_x[ipt] = aList_symm_x[l]
                    aList_beads_y[ipt] = aList_symm_y[l]
                    aList_beads_z[ipt] = aList_symm_z[l]
                    l = l + 1
    
                i = 0
                while i < nbeads:
                    if i != ii:
                        x = float(aList_beads_x[i])
                        y = float(aList_beads_y[i])
                        z = float(aList_beads_z[i])
                        dr_new = get_dr(xtest,ytest,ztest,x,y,z)
                        vdw_new = vdw_energy(econ12,econ6,e_width,dr_new,bead_sep3)
                        energy_new = energy_new + num_symm*vdw_new
                    i = i + 1
    
                # Get cavity energy difference
    
                delta_old = center_beads(xold,yold,zold,aList_beads_x,\
                                         aList_beads_y,aList_beads_z,radii_1,radii_2)
                delta_new = center_beads(xtest,ytest,ztest,aList_beads_x,\
                                         aList_beads_y,aList_beads_z,radii_1,radii_2)
    
                delta_pack = num_symm*surface_scale*(delta_new - delta_old)/(radii_1 + radii_2)
                sum_delta_pack = sum_delta_pack + abs(delta_pack)
    
                # Get shifted contributions to P(r)  
    
                aList_pr_model_test2 = num_hist*[0.0,]
                pr_shift_atom(aList_pr_model_test2,xtest,ytest,ztest,aList_beads_x,\
                              aList_beads_y,aList_beads_z,hist_grid,ii)
    
                # Get net shift in contribution to P(r)
            
                i = 0
                while i < num_hist:
                    aList_pr_model_test[i] = aList_pr_model_test2[i] - aList_pr_model_test[i]
                    i = i + 1
    
                # Get statistic for agreement with P(r) after accumulating shift vectors
    
                i = 0
                while i < num_hist:
                    aList_pr_model_test[i] = float(aList_pr_model[i]) + num_symm*float(aList_pr_model_test[i])
                    i = i + 1
    
                hist_score = pr_dif(aList_pr,aList_pr_model_test,skip)
                
#                aList_i_calc = []
#                ft_to_intensity(aList_q,aList_i_calc,aList_r,aList_pr_model,nbeads)
#                (chi_sq,rvalue) = score_Ic(aList_q,aList_i,aList_i_sd,aList_i_calc)
    
                # Scoring shifts
    
                delta_h = hist_score - hist_score_best
                delta_e = energy_new - energy_old + delta_pack
    
                # Recalibrate scale so impact of energy and P(r) is equal on plausible shifts
            
                if delta_e < abs_e_width:
                    sum_e = sum_e + abs(delta_e)
                    sum_h = sum_h + abs(delta_h)
    
                # Monitor potential moves based of P(r)
    
                if delta_h < 0.0:
                    count_hist_yes = count_hist_yes + 1.0            
    
                # Acceptance and update
    
                score = delta_e + delta_h*hscale
    
                if count_it < num_sa_max:
                    barrier = abs(random.gauss(0.0,sd_mc))
                else:
                    barrier = 0.0
     
                if score < barrier:
    
                    success = success + 1.0
                    hist_score_best = hist_score
    
                    # Update model but symmetry positions that were already put in
                
                    aList_beads_x[ii] = xtest
                    aList_beads_y[ii] = ytest
                    aList_beads_z[ii] = ztest
    
                    # Update P(r)
    
                    i = 0
                    while i < num_hist:
                        aList_pr_model[i] = aList_pr_model_test[i]
                        i = i + 1
    
                else:
    
                    # Revert to original model
                
                    aList_beads_x[ii] = xold
                    aList_beads_y[ii] = yold
                    aList_beads_z[ii] = zold
    
                    l = 0
                    while l < num_ops:
                        ipt = ii + l + 1
                        aList_beads_x[ipt] = aList_temp_save_x[l]
                        aList_beads_y[ipt] = aList_temp_save_y[l]
                        aList_beads_z[ipt] = aList_temp_save_z[l]              
                        l = l + 1     
                #
    
                ii = ii + num_symm
    
            # Get energy statistics at end of macrocycle
    
            vdw_all = get_total_energy(aList_beads_x,aList_beads_y,aList_beads_z,econ12,econ6,bead_sep3)
    
            # Rescale and convergence statistics
    
            if sum_h > 0.0:
                hscale = sum_e/sum_h
            
            count_hist_yes = count_hist_yes*float(num_symm)/float(nbeads)
            success_rate = success*float(num_symm)/float(nbeads)
            success_rate_all = success_rate_all + success_rate

            if not (count_it+1)%10:   
                aString = 'Cycle ' + str(count_it+1) + ' Moves ' + str('%.2f'%(success_rate)) + \
                          ' Possibles ' + str('%.2f'%(count_hist_yes)) + ' rms P(r) '+ str('%4.3f'%(hist_score)) + \
                          ' Energy ' + str('%4.2f'%(vdw_all))
                print (aString)
#                print('Rvalue: %4.3f CHI-squared: %4.3f'%(rvalue,chi_sq))

            if dlg:
                dlg.Update(count_it+1,newmsg='Cycle no.: '+str(count_it)+' of 160')
    
            # Debug statitics. Weight of 10 gives about 1.0
            #sum_delta_pack = sum_delta_pack/float(nbeads)
            #print (sum_delta_pack)
            #
    
            count_it = count_it + 1
    
        #####################
        # ANALYZE AND WRITE #
        #####################
    
        aString = '\nFinal model statistics'
        print (aString)
    
        calc_pr(aList_beads_x,aList_beads_y,aList_beads_z,aList_pr_model,hist_grid)
        hist_score_best = pr_dif(aList_pr,aList_pr_model,skip)
    
        # P(r) fitting statistics
        delta_hist_sum = pr_rfactor(aList_pr,aList_pr_sd,aList_pr_model_test,skip)
    
        aString = 'Delta P(r): ' + str('%4.3f'%(delta_hist_sum))
        print (aString)
    
        # Get final energy 
        vdw_all = get_total_energy(aList_beads_x,aList_beads_y,aList_beads_z,econ12,econ6,bead_sep3)
    
        aString = 'VDW energy: ' + str('%4.2f'%(vdw_all))
        print (aString)
    
        Phases.append([file_beads.split('.')[0],aList_beads_x,aList_beads_y,aList_beads_z])
        # Write out beads as pseudo a PDB file
        if pdbOut:
            pdb_writer(aList_beads_x,aList_beads_y,aList_beads_z,file_beads,aString)
    
        # Calculate and write final PSV shape
    
        aList_box_x_all = []
        aList_box_y_all = []
        aList_box_z_all = []
        aList_box_score = []
    
        set_box_fast(aList_beads_x,aList_beads_y,aList_beads_z,\
                aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                aList_box_score,box_step,dmax,rsearch)
    
        aList_box_x = []
        aList_box_y = []
        aList_box_z = []
        psv_vol_use = psv_vol*inflate
     
        set_vol(aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                aList_box_score,aList_box_x,aList_box_y,aList_box_z,\
                psv_vol_use,box_pt_vol)
    
        num_boxes = len(aList_box_x)
        fraction_psv = num_boxes*box_pt_vol/psv_vol
    
        # Correct final volume if initial estimate is too small
    
        fraction_psv_use = num_boxes*box_pt_vol/psv_vol_use
        
        if fraction_psv_use < 1.0:
            aList_box_x = []
            aList_box_y = []
            aList_box_z = []
            vol_use = 1.05*psv_vol_use
            set_vol(aList_box_x_all,aList_box_y_all,aList_box_z_all,\
                    aList_box_score,aList_box_x,aList_box_y,aList_box_z,\
                    vol_use,box_pt_vol)
    
            num_boxes = len(aList_box_x)
            fraction_psv = float(num_boxes)*box_pt_vol/psv_vol    
    
        aString = 'Final PSV of protein envelope: ' + str('%4.2f'%(fraction_psv))
        print (aString)
    
        # Write input and model P(r)
#        pr_writer(aList_pr,aList_r,aList_pr_model,file_pr)
        PRcalc.append([aList_r,aList_pr,copy.copy(aList_pr_model),delta_hist_sum])
    
        # Calculate comparison versus intensities
    
        if num_intensities > 0:
    
            # calculate intensity
            aList_i_calc = []
#            num_beads = len(aList_box_x)
            ft_to_intensity(aList_q,aList_i_calc,aList_r,aList_pr_model,nbeads)
    
            # Scale and obtain statistics
            (chi_sq,rvalue) = score_Ic(aList_q,aList_i,aList_i_sd,aList_i_calc)
            
            aString = 'Rvalue: ' + str('%4.3f'%(rvalue)) + ' CHI-squared: ' + str('%4.3f'%(chi_sq))
            print (aString)
    
            # Write output intensity file
            Patterns.append([aList_q,aList_i,aList_i_calc,rvalue])
#            write_all_data(file_intensity,aList_q,aList_i,aList_i_calc,aString)
    
#            aString = 'Output intensity file: ' + str(file_intensity)
#            print (aString)
            
        else:
    
            chi_sq = 'NA'
    
        # Write final volume
    
        delta_hist_sum = '%4.3f'%(delta_hist_sum)
        vdw_all = '%4.2f'%(vdw_all)
        fraction_psv = '%4.2f'%(fraction_psv)
        chi_sq = '%4.3f'%(chi_sq)
    
        aString = 'REMARK     P(r) dif:'+ str(delta_hist_sum) + ' E:'\
                  + str(vdw_all) + ' CHIsq:' + str(chi_sq) + \
                  ' PSV:' + str(fraction_psv)
    
        Phases.append([file_psv.split('.')[0],aList_box_x,aList_box_y,aList_box_z])
        if pdbOut:
            pdb_writer(aList_box_x,aList_box_y,aList_box_z,file_psv,aString)
    
        # Write Summary
    
        aString = 'P(r) dif:' + str(delta_hist_sum) + ' E:' \
                         + str(vdw_all) + ' CHISQ:' + str(chi_sq) + ' PSV:' + str(fraction_psv)
#        aList_summary.append(aString)               
    
        i_soln = i_soln + 1
    
    #########################################
    #### END OF LOOP OVER MULTI-SOLUTIONS ###
    #########################################
    
    #end profiling            
    pr.disable()
    s = StringIO.StringIO()
    sortby = 'tottime'
    ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
    print('Profiler of function calculation; top 50% of routines:')
    ps.print_stats(.5)
    print(s.getvalue())
    print('%s%.3f'%('Run time = ',time.time()-time0))
            
    localtime = time.asctime( time.localtime(time.time()) )
    
    aString = '\nCompletion time: ' + str(localtime)
    print (aString)
                         
#    aList_summary.append(str(localtime))
#    
#    # Create summary
#    
#    aFile_log = prefix + 'shapes_summary.log'
#    num_lines = len(aList_summary)
#    
#    file = open(aFile_log,'w')
#    i = 0
#    while i < num_lines:
#        aString = str(aList_summary[i])
#        file.write(aString)
#        file.write('\n')
#        i = i + 1
#    file.close()

    return Phases,Patterns,PRcalc
    
    
        
