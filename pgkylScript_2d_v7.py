#.Make plots from Gkyl data.
#.Manaure Francisquez (base) and Lucio Milanese (updates and extensions).
#.Spring 2019.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
import postgkyl as pg
import numpy as np
import adios as ad
import sys
sys.path.insert(0, '/Users/rlw/Dropbox/ForBB/LucioFiles/pgkyl_script/')
from shutil import copyfile
import pgkylUtil as pgu
import os
import subprocess
import shlex
import re
directories = ['.']

#rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})
#params = {'backend': 'pdf',
#          'text.usetex': True,
#          'axes.unicode_minus': True}
#matplotlib.rcParams.update(params)

fileNames   = ['Biskamp_newNoise_FixedJ_Tratio100']    #.Root name of files to process.

outDir     = './gkyl_plots_v7/'               #.Output directory where figures are saved.
fourier_transform = True
auto_loading      = False
#creating the directory for plots if it does not exist yet
if not os.path.exists('gkyl_plots_v7'):
    os.makedirs('gkyl_plots_v7')


def main_plotting(fileroot, max_output):

    fileName = fileroot
    iFrame     = 0     #.Initial frame (0 is the t=0 frame).
    fFrame     = max_output    #.Final frame.

    polyOrder  = 2
    basisType  = 'ms'
    m_ion = 25
    vTe0 = 0.02
    alpha = 0.00

    #.Component of the quantity we wish to extract from data file.
    #.For field files this specifies the field component (e.g. Ex,
    #.Ey, Ez, Bx, By, or Bz) while for Mi1 it specifies the vector
    #.component of the momentum density.
    compZero = 0

    #..................... NO MORE USER INPUTS BELOW (maybe) ....................#

    nFrames = fFrame-iFrame+1    #.Number of frames plotted.

    #.Extract grid details from one of the data files for each species
    fName_elc = fileName+'_elc_'+str(iFrame)+'.bp'
    fName_ion = fileName+'_ion_'+str(iFrame)+'.bp'

    # getGrid data
    x_elc, _, nx, lx, _ = pgu.getGrid(fName_elc,polyOrder,basisType,location='center')
    x_ion, _, _, _, _ = pgu.getGrid(fName_ion,polyOrder,basisType,location='center')
            #x_e, gridDim, nx, lx, dx = pgu.getGrid(fName,polyOrder,basisType,location='center')
        	#    x_e = np.array(x_e)  #RLW: Not needed

    #Store needed data from getGrid
    nz = nx[0]
    ny = nx[1]
    lz = lx[0]  #get physical z location of center of box
    ly = lx[1]
	#information about the X grid (same for both species)
    points_z = np.linspace(-lz/2, lz/2, nz)
    points_y = np.linspace(-ly/2, ly/2, ny)
    #information about the V grid for the electrons
    vz_elc_min = x_elc[2][0]
    vz_elc_max = x_elc[2][-1]
    vy_elc_min = x_elc[3][0]
    vy_elc_max = x_elc[3][-1]
    #information about the V grid for the ions
    vz_ion_min = x_ion[2][0]
    vz_ion_max = x_ion[2][-1]
    vy_ion_min = x_ion[3][0]
    vy_ion_max = x_ion[3][-1]

    # getInterpData data
    elcd = np.squeeze(pgu.getInterpData(fName_elc,polyOrder,basisType,comp=compZero))
    iond = np.squeeze(pgu.getInterpData(fName_ion,polyOrder,basisType,comp=compZero))

    #obtaining information about the grid in velocity space for electrons
    z0_elc = int(elcd.shape[0]/2)  #rlw: get z coordinate of center of box (it is 48)
    y0_elc = int(elcd.shape[1]/2)
    vzsize_elc = int(elcd.shape[2])
    vysize_elc = int(elcd.shape[3])
    velocitiesz_elc = np.linspace(vz_elc_min, vz_elc_max, vzsize_elc)
    velocitiesy_elc = np.linspace(vy_elc_min, vy_elc_max, vysize_elc)
    #obtaining information about the grid in velocity space for ions
    z0_ion = int(iond.shape[0]/2)
    y0_ion = int(iond.shape[1]/2)
    vzsize_ion = int(iond.shape[2])
    vysize_ion = int(iond.shape[3])
    velocitiesz_ion = np.linspace(vz_ion_min, vz_ion_max, vzsize_ion)
    velocitiesy_ion = np.linspace(vy_ion_min, vy_ion_max, vysize_ion)



    #Setting up the grids in V space for plotting
    Vz_elc, Vy_elc = np.meshgrid(velocitiesz_elc,velocitiesy_elc,indexing='ij')
    Vz_ion, Vy_ion = np.meshgrid(velocitiesz_ion,velocitiesy_ion,indexing='ij')
    #Setting up the grids in X space for plotting
    grid_z, grid_y = np.meshgrid(points_z, points_y,indexing='ij')


	#initialize all arrays containing time-dependent quantities...that are going to be filled in during loop over frames
    times = np.zeros(nFrames)
    eField_boxavg_z = np.zeros(nFrames)
    Rei = np.zeros(nFrames)
    J_boxavg_z = np.zeros(nFrames)
    eField_fluct_squared   = np.zeros(nFrames)
    eField_fluct_squared_byT = np.zeros(nFrames)
    E_over_J_rolling = np.zeros(nFrames)
    R_over_J_rolling = np.zeros(nFrames)
    energy_e_tot = np.zeros(nFrames)
    energy_b_tot = np.zeros(nFrames)
    E_boxavg_over_J_boxavg = np.zeros(nFrames)
    Tion_boxavg = np.zeros(nFrames)
    Telc_boxavg = np.zeros(nFrames)
    nu_Sagdeev = np.zeros(nFrames)

#    Tion_boxavg_over_elc = np.zeros(nFrames)

#    arrays the used to be initialized that did not need to be (though it may be slower)
#    currentDen = np.zeros((nz,ny))
#    Den        = np.zeros((nz,ny))
#    ek_t = np.zeros((nFrames,int(nz/2+1),int(ny/2+1)))

    #this is the electron distribution function at a specific location in space as a function of time
#    elcd_x0t = np.zeros((nFrames,vzsize_e,vysize_e))  #Never used!  Probably meant to be elcd_cut
#    iond_x0t = np.zeros((nFrames,vzsize_i,vysize_i))
    for nFr in range(iFrame,fFrame+1):
        fName_elc = fileName+'_elc_'+str(nFr)+'.bp'
        fNameM0_elc = fileName+'_elc_M0_'+str(nFr)+'.bp'
        fNameM1_elc = fileName+'_elc_M1i_'+str(nFr)+'.bp'        #.Complete file name.
        fName_vT_elc = fileName+'_elc_vthSq_'+str(nFr)+'.bp'
        fName_u_elc = fileName+'_elc_u_'+str(nFr)+'.bp'

#        fNameM2_elc = fileName+'_elc_intM2Thermal_'+str(nFr)+'.bp'
#        for JJ: what is in the vthSqCross file?
#        for JJ: what is in the vthSq file?_uCross, u_, _intM2Flow, _intL2
#       intM2Flow = \int dx n*(u^2)
#       intL2 = \int dx \int dv f^2
#       vthSq = v_t^2
#       u = u, the mean flow
#       cross ones are special for cross-collisions. you probably don’t need them but their formulas are in https://gkeyll.readthedocs.io/en/latest/dev/collisionmodels.html (2nd and 3rd equations below Dougherty collisions if you are using LBO collisions)

        fName_ion = fileName+'_ion_'+str(nFr)+'.bp'
        fNameM0_ion = fileName+'_ion_M0_'+str(nFr)+'.bp'
        fNameM1_ion = fileName+'_ion_M1i_'+str(nFr)+'.bp'        #.Complete file name.
        fName_vT_ion = fileName+'_ion_vthSq_'+str(nFr)+'.bp'
        fName_u_ion = fileName+'_ion_u_'+str(nFr)+'.bp'

        fName_field = fileName+'_field_'+str(nFr)+'.bp'

        elcd = np.squeeze(pgu.getInterpData(fName_elc,polyOrder,basisType))
        elcM0 = np.squeeze(pgu.getInterpData(fNameM0_elc,polyOrder,basisType,comp=0))	# JJ: are we sure that we have to specify the keyword comp?
        elcM1_z = np.squeeze(pgu.getInterpData(fNameM1_elc,polyOrder,basisType,comp=0))
        elcM1_y = np.squeeze(pgu.getInterpData(fNameM1_elc,polyOrder,basisType,comp=1))
        elc_u_z = np.squeeze(pgu.getInterpData(fName_u_elc,polyOrder,basisType,comp=0))
        elc_u_y = np.squeeze(pgu.getInterpData(fName_u_elc,polyOrder,basisType,comp=1))
        elcM1_raw = np.squeeze(pgu.getRawData(fNameM1_elc))
        elc_vTsq = np.squeeze(pgu.getInterpData(fName_vT_elc, polyOrder, basisType))
        elc_vT = np.sqrt(elc_vTsq)
        elcT = elc_vTsq/((vTe0)**2)  #This is an array of temperature values: one for each spatial location
        elcT_boxavg = np.average(elcT)
#        elcM2 = np.squeeze(pgu.getRawData(fNameM2_elc))


        iond = np.squeeze(pgu.getInterpData(fName_ion,polyOrder,basisType))
        ionM0 = np.squeeze(pgu.getInterpData(fNameM0_ion,polyOrder,basisType,comp=0))	# JJ: are we sure that we have to specify the keyword comp
        ionM1_z = np.squeeze(pgu.getInterpData(fNameM1_ion,polyOrder,basisType,comp=0))
        ionM1_y = np.squeeze(pgu.getInterpData(fNameM1_ion,polyOrder,basisType,comp=1))
        ion_u_z = np.squeeze(pgu.getInterpData(fName_u_ion,polyOrder,basisType,comp=0))
        ion_u_y = np.squeeze(pgu.getInterpData(fName_u_ion,polyOrder,basisType,comp=1))
        ionM1_raw = np.squeeze(pgu.getRawData(fNameM1_ion))
        ion_vTsq = np.squeeze(pgu.getInterpData(fName_vT_ion, polyOrder, basisType))
        ionT = m_ion * ion_vTsq/((vTe0)**2)  #This is an array of temperature values: one for each spatial location
        ionT_boxavg = np.average(ionT)



        eField_z = np.squeeze(pgu.getInterpData(fName_field,polyOrder,basisType,comp=0)) # JJ: can this be turned into 1 call without comp specified?
        eField_y = np.squeeze(pgu.getInterpData(fName_field,polyOrder,basisType,comp=1))
#        e_raw = np.squeeze(pgu.getRawData(fName_field))
        #fName     = fileName + '_ion_M1i_'+str(nFr)+'.bp'    #.Complete file name.
        #fName_den = fileName + '_ion_M0_'+str(nFr)+'.bp'    #.Complete file name.

        #ionM1 = np.squeeze(pgu.getInterpData(fName,polyOrder,basisType,comp=compZero))
       # ionM0 = np.squeeze(pgu.getInterpData(fName_den,polyOrder,basisType,comp=compZero))

        #ionM1_raw = np.squeeze(pgu.getRawData(fName))

        #fName = fileName+'_elc_'+str(nFr)+'.bp'


       # fName = fileName+'_ion_'+str(nFr)+'.bp'
       # iond = np.squeeze(pgu.getInterpData(fName,polyOrder,basisType))

        #compute box-averaged distribution function
        elcd_box_avg = np.zeros((vzsize_elc,vysize_elc)) #needs to be initialized b/c it is referenced a few lines down
        iond_box_avg = np.zeros((vzsize_ion,vysize_ion))

        #for k in range(nz):
        #    for j in range(ny):
        #        elcd_box_avg[:,:] += elcd[k,j,:,:]
        #        iond_box_avg[:,:] += iond[k,j,:,:]
        #elcd_box_avg = elcd_box_avg/(ny*nz)
        #iond_box_avg = iond_box_avg/(ny*nz)
        elcd_box_avg = np.average(elcd, axis = (0,1)) #average over both spatial dimensions
        iond_box_avg = np.average(iond, axis = (0,1))

        elcd_cut = elcd[z0_elc,y0_elc,:,:]
        iond_cut = iond[z0_ion,y0_ion,:,:]

       # eName = fileName+'_field_'+str(nFr)+'.bp'

        #temperature plotting
        #Ti_fname = fileName + '_ion_' + 'intM2Thermal_' + str(nFr) + '.bp'
        #ionM2 = np.squeeze(pgu.getRawData(Ti_fname))
        #        ionM2 = np.squeeze(pgu.getRawData(Ti_fname,polyOrder,basisType,comp=compZero))
        #Jincreasing_simplest_randomkicksfunJ_Tratio100_standard2D_RLW_LMM_edit3_diag_test_ion_intM2Thermal

        #loading electric field

        #split into spatial average and fluctuating parts
        boxavg_eField_z = np.average(eField_z)
        boxavg_eField_y = np.average(eField_y)
        eField_fluct_z = eField_z - boxavg_eField_z
        eField_fluct_y = eField_y - boxavg_eField_y

        boxavg_uElc_z = np.average(elc_u_z)
        boxavg_uElc_y = np.average(elc_u_y)

        #e_raw will contain Ex, Ey, Ez

#        e_z_raw_cell_average = eField_fluct_z#0.5*e_raw[:,:,0]  #rlw: why multiply by 0.5?
#        e_y_raw_cell_average = eField_fluct_y#0.5*e_raw[:,:,8]
#        nz_avg          = e_z_raw_cell_average.shape[0]
#        ny_avg          = e_z_raw_cell_average.shape[1]

        #compute energy in electric and magnetic fields field
        Tion_boxavg[nFr-iFrame] = ionT_boxavg
        Telc_boxavg[nFr-iFrame] = elcT_boxavg
        eField_fluct_squared[nFr-iFrame] = (np.sum(eField_fluct_z**2 + eField_fluct_y**2) /(nz*ny) ) / vTe0**2
        eField_fluct_squared_byT[nFr-iFrame] = (np.sum(eField_fluct_z**2 + eField_fluct_y**2) /(nz*ny) ) / (vTe0*elcT_boxavg)**2
        #this would be the cell-averaged values (?)
#        energy_e_tot[nFr-iFrame] = 0.5*(np.sum(e_raw[:,:,0]**2)+ np.sum(e_raw[:,:,8]**2) + np.sum(e_raw[:,:,16]**2))
#        energy_b_tot[nFr-iFrame] = 0.5*(np.sum(e_raw[:,:,24]**2)+ np.sum(e_raw[:,:,32]**2) + np.sum(e_raw[:,:,40]**2))

        #.Compute the kz=0 component of the electric field
        eField_boxavg_z[nFr-iFrame] = boxavg_eField_z
        #LMM (should this maybe not be normalized?)

        #.Compute current density.
        Jz     = ionM1_z     - elcM1_z
        Jy     = ionM1_y     - elcM1_y
        Den            = ionM0      #Mass density - ought to subtract the electrons (?) LMM
        J_fluct_z = Jz - np.sum(Jz)/(nz*ny)
        J_fluct_y = Jy - np.sum(Jy)/(nz*ny)
#        currentDen_raw_z = 0.5*(ionM1_raw[:,:,0] - elcM1_raw[:,:,0])
#        currentDen_raw_y = 0.5*(ionM1_raw[:,:,8] - elcM1_raw[:,:,8])  #Is 8 the right component??

        #updating J_boxavg_z
        J_boxavg_z[nFr-iFrame] = np.sum(Jz)/(nz*ny)

                #.Extract the time from file.
        hF         = ad.file(fName_elc)
        times[nFr-iFrame] = hF['time'].read()
        hF.close()
        time = float('%.3g' % times[nFr-iFrame])  #fix me please

        fignum = str(nFr).zfill(4)
        

        #Distribution function plots [Dplots]
        fig, axs = plt.subplots(2,3,figsize=(40, 20), facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace = .5, wspace =.1)
        axs = axs.ravel()

        pos0 = axs[0].pcolormesh(Vz_elc, Vy_elc, elcd_cut)
        axs[0].set_xlabel(r'$v_z/c$', fontsize=30)
        axs[0].set_ylabel(r'$v_y/c$', fontsize=30, labelpad=-1)
        axs[0].set_title(r'$F_e(z_0,y_0,v_z,v_y),$' + rf'$J_z$ = {np.average(Jz)}/vTe0'+ r' [$c_{s0}$]', fontsize=26)
        axs[0].tick_params(labelsize = 26)
        cbar = fig.colorbar(pos0, ax=axs[0])
        cbar.ax.tick_params(labelsize=22)

        pos1 = axs[1].pcolormesh(Vz_ion, Vy_ion, iond_cut)
        axs[1].set_xlabel(r'$v_z/c$', fontsize=30)
        axs[1].set_ylabel(r'$v_y/c$', fontsize=30, labelpad=-1)
        axs[1].set_title(r'$F_i(z_0,y_0,v_z,v_y),$' + rf't = {time}'+ r' [$\omega_{pe}^{-1}$]', fontsize=26)
        axs[1].tick_params(labelsize = 26)
        cbar = fig.colorbar(pos1, ax=axs[1])
        cbar.ax.tick_params(labelsize=22)

        pos2 = axs[2].pcolormesh(grid_y, grid_z, Jz)
        axs[2].set_xlabel(r'$y \ [d_e]$', fontsize=30)
        axs[2].set_ylabel(r'$z \ [d_e]$', fontsize=30, labelpad=-1)
        axs[2].set_title(r'$J(z,y)$', fontsize=30)
        axs[2].tick_params(labelsize = 26)
        cbar = fig.colorbar(pos2, ax=axs[2])
        cbar.ax.tick_params(labelsize=22)

        pos3 = axs[3].pcolormesh(Vz_elc, Vy_elc, elcd_box_avg)
        axs[3].scatter(boxavg_uElc_z, boxavg_uElc_y, s = 60)
        axs[3].scatter(np.squeeze(Vz_elc[np.where(elcd_box_avg==np.max(elcd_box_avg))]),np.squeeze(Vy_elc[np.where(elcd_box_avg==np.max(elcd_box_avg))]),s = 40, marker = 'x', alpha = 1)
        axs[3].set_xlabel(r'$v_z/c$', fontsize=30)
        axs[3].set_ylabel(r'$v_y/c$', fontsize=30, labelpad=-1)
        axs[3].set_title(r'$<F_e(v_z,v_y)>_{z,y},$' + rf't = {time}'+ r' [$\omega_{pe}^{-1}$]', fontsize=26)
        axs[3].tick_params(labelsize = 26)
        cbar = fig.colorbar(pos3, ax=axs[3])
        cbar.ax.tick_params(labelsize=22)

        pos4 = axs[4].pcolormesh(Vz_ion, Vy_ion, iond_box_avg)
        axs[4].set_xlabel(r'$v_z/c$', fontsize=30)
        axs[4].set_ylabel(r'$v_y/c$', fontsize=30, labelpad=-1)
        axs[4].set_title(r'$<F_i(v_z,v_y)>_{z,y},$' + rf't = {time}'+ r' [$\omega_{pe}^{-1}$]', fontsize=26)
        axs[4].tick_params(labelsize = 26)
        cbar = fig.colorbar(pos4, ax=axs[4])
        cbar.ax.tick_params(labelsize=22)

        pos5 = axs[5].pcolormesh(grid_y, grid_z, Den)
        axs[5].set_xlabel(r'$y \ [d_e]$', fontsize=30)
        axs[5].set_ylabel(r'$z \ [d_e]$', fontsize=30, labelpad=-1)
        axs[5].set_title('$\\rho (z,y)$', fontsize=30)
        axs[5].tick_params(labelsize = 26)
        cbar = fig.colorbar(pos5, ax=axs[5])
        cbar.ax.tick_params(labelsize=22)

        fig.tight_layout()
        plt.savefig(outDir+fileName+rf'_diagnostics_{fignum}.png', bbox_inches='tight')
        plt.close()

        #computing the fourier transform of the electric field
        if(fourier_transform):
#            e_z_raw_cell_average_k = np.zeros((nz_avg, ny_avg),dtype=complex)
#            J_z_k_raw              = np.zeros((nz_avg, ny_avg),dtype=complex)
            eField_fluct_z_k = np.fft.fftn(eField_fluct_z/vTe0)
            eField_fluct_y_k = np.fft.fftn(eField_fluct_y/vTe0)
            #print('fft_freq', np.fft.fftfreq(128) )
            J_fluct_z_k              = np.fft.fftn(J_fluct_z)
            J_fluct_y_k              = np.fft.fftn(J_fluct_y)

            #e_z_raw_cell_average_k_y_int and J_z_k_raw_y_int are integrated over y
#            e_z_raw_cell_average_k_y_int = np.zeros(nz_avg)
#            J_z_k_raw_y_int                 = np.zeros(nz_avg)

#            e_z_raw_cell_average_k[0,0] = 0 #test: try to remove box-averaged component
#            e_y_raw_cell_average_k[0,0] = 0
#            J_z_k_raw[0,0] = 0
#            J_y_k_raw[0,0] = 0
            JdotE_k = np.abs(np.transpose(np.fft.fftshift(J_fluct_z_k*eField_fluct_z_k + J_fluct_y_k*eField_fluct_y_k) ) )
            eField_fluct_square_K = np.abs(np.transpose(np.fft.fftshift(eField_fluct_z_k**2 + eField_fluct_y_k**2) ) )
            #integrating over y direction (RLW: can we do this in a vectorized way?)
            #for j in range(ny_avg):
            #    for i in range(nz_avg):
#            for j in range(ny):
#                for i in range(nz):
#                    e_z_raw_cell_average_k_y_int[i] += np.abs(e_z_raw_cell_average_k[i,j])
#                    J_z_k_raw_y_int[i]              += np.abs(J_z_k_raw[i,j])


            #ek_t[nFr,:] = ek
            fignum = str(nFr).zfill(4)

            #z_plot_2d_sp   = np.linspace(-int(nz_avg/2), int(nz_avg/2-1), nz_avg)/lz
            #y_plot_2d_sp   = np.linspace(-int(ny_avg/2), int(ny_avg/2-1), ny_avg)/ly

            kz_plot_2d_sp   = 2.0*3.14159*vTe0*np.linspace(-int(nz/2), int(nz/2-1), nz)/lz
            ky_plot_2d_sp   = 2.0*3.14159*vTe0*np.linspace(-int(ny/2), int(ny/2-1), ny)/ly
            K_z, K_y = np.meshgrid(kz_plot_2d_sp, ky_plot_2d_sp, indexing = 'xy')


            #Plotting FFT of electric field and J

            fig, axs = plt.subplots(2,2,figsize=(20, 20), facecolor='w', edgecolor='k')
            fig.subplots_adjust(hspace = .5, wspace =.1)
            axs = axs.ravel()
            pos0 = axs[0].pcolormesh(Vz_elc, Vy_elc, elcd_box_avg)
            axs[0].scatter(boxavg_uElc_z, boxavg_uElc_y, s = 60)
            axs[0].scatter(np.squeeze(Vz_elc[np.where(elcd_box_avg==np.max(elcd_box_avg))]),np.squeeze(Vy_elc[np.where(elcd_box_avg==np.max(elcd_box_avg))]),s = 40, marker = 'x', alpha = 1)
            axs[0].set_xlabel(r'$v_z/c$', fontsize=30)
            axs[0].set_ylabel(r'$v_y/c$', fontsize=30, labelpad=-1)
            axs[0].set_title(r'$<F_e(v_z,v_y)>_{z,y},$' + rf't = {time}'+ r' [$\omega_{pe}^{-1}$]', fontsize=26)
            axs[0].tick_params(labelsize = 26)
            cbar = fig.colorbar(pos0, ax=axs[0])
            cbar.ax.tick_params(labelsize=22)

            pos1 = axs[1].pcolormesh(Vz_ion, Vy_ion, iond_box_avg)
            axs[1].set_xlabel(r'$v_z/c$', fontsize=30)
            axs[1].set_ylabel(r'$v_y/c$', fontsize=30, labelpad=-1)
            axs[1].set_title(r'$<F_e(v_z,v_y)>_{z,y},$' + rf't = {time}'+ r' [$\omega_{pe}^{-1}$]', fontsize=26)
            axs[1].tick_params(labelsize = 26)
            cbar = fig.colorbar(pos1, ax=axs[1])
            cbar.ax.tick_params(labelsize=22)
#            pos0 = axs[0].plot(e_z_raw_cell_average_k_y_int)
#            axs[0].set_xlabel(r'$k_z \lambda_{De0}$', fontsize=30)
#            axs[0].set_ylabel(r'$E(k_z)$', fontsize=30, labelpad=-1)
#            axs[0].set_title(rf't = {time}'+ r'[$\omega_{pe}^{-1}$]', fontsize=30)
#            axs[0].tick_params(labelsize = 26)
#            axs[0].set_yscale('symlog') #plot on log scale and take absolute value of number

#            pos1 = axs[1].plot(J_z_k_raw_y_int)
#            axs[1].set_xlabel(r'$k_z d_e$', fontsize=30)
#            axs[1].set_ylabel(r'$J(k_z)$', fontsize=30, labelpad=-1)
#            axs[1].set_title(rf't = {time}'+ r'[$\omega_{pe}^{-1}$]', fontsize=30)
#            axs[1].tick_params(labelsize = 26)
#            axs[1].set_yscale('symlog')

            pos2 = axs[2].contourf(K_z, K_y, eField_fluct_square_K)
            axs[2].set_xlabel(r'$k_z \lambda_{De0}$', fontsize=30)
            axs[2].set_ylabel(r'$k_y \lambda_{De0}$', fontsize=30, labelpad=-1)
            axs[2].set_title(rf'$|\delta E^2|_k$, t = {time}'+ r'[$\omega_{pe}^{-1}$]', fontsize=30)
            axs[2].tick_params(labelsize = 26)
            cbar = fig.colorbar(pos2, ax=axs[2])
            cbar.ax.tick_params(labelsize=22)

            pos2 = axs[3].contourf(K_z, K_y, JdotE_k)
            axs[3].set_xlabel(r'$k_z \lambda_{De0}$', fontsize=30)
            axs[3].set_ylabel(r'$k_y \lambda_{De0}$', fontsize=30, labelpad=-1)
            axs[3].set_title(rf'$(\delta J \cdot \delta E)_k$, t = {time}'+ r'[$\omega_{pe}^{-1}$]', fontsize=30)
            axs[3].tick_params(labelsize = 26)
            cbar = fig.colorbar(pos3, ax=axs[3])
            cbar.ax.tick_params(labelsize=22)
#            pos3 = axs[3].contourf(K_z, K_y, JdotE_k )
#            axs[3].set_xlabel(r'$k_z d_e$', fontsize=30)
#            axs[3].set_ylabel(r'$k_y d_e$', fontsize=30, labelpad=-1)
#            axs[3].set_title(rf'$(\delta J \cdot \delta E)_k$, t = {time}'+ r'[$\omega_{pe}^{-1}$]', fontsize=30)
#            axs[3].tick_params(labelsize = 26)
#            cbar = fig.colorbar(pos3, ax=axs[3])
#            cbar.ax.tick_params(labelsize=22)

            fig.tight_layout()
            plt.savefig(outDir+fileName+rf'_fft_{fignum}.png', bbox_inches='tight')
            plt.close()

            #PLOTTING (to be done at each time step)

#### we have now left the loop over frames

    Navg = 3
    Rei = eField_boxavg_z + alpha*vTe0
#    E_over_J_rolling = np.zeros(nFrames)
    for n in range(Navg,nFrames):
        E_over_J_rolling[n] = np.sum(eField_boxavg_z[n-Navg:n])/np.sum(J_boxavg_z[n-Navg:n])
    for n in range(Navg):
        E_over_J_rolling[n] =  E_over_J_rolling[Navg]  #bfill the first Navg values


    for n in range(Navg,nFrames):
        R_over_J_rolling[n] = np.sum(Rei[n-Navg:n])/np.sum(J_boxavg_z[n-Navg:n])
    for n in range(Navg):
        R_over_J_rolling[n] =  R_over_J_rolling[Navg]  #bfill the first Navg values
    nu_eff_rolling = -(R_over_J_rolling)

    nu_Sagdeev = 0.025*np.absolute(J_boxavg_z/np.sqrt(Telc_boxavg))*Telc_boxavg/Tion_boxavg

    #Energy plots
    fig, axs = plt.subplots(2,3,figsize=(30, 20), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace =.1)
    axs = axs.ravel()
    np.save('timesFile',times[0:fFrame-iFrame+1])
    np.save('eZboxAvgFile',eField_boxavg_z[0:fFrame-iFrame+1])
    np.save('jZboxAvgFile',J_boxavg_z[0:fFrame-iFrame+1])
    np.save('eSquaredFile',eField_fluct_squared[0:fFrame-iFrame+1])
    np.save('Tion_boxavgFile',Tion_boxavg[0:fFrame-iFrame+1])
    np.save('Telc_boxavgFile',Telc_boxavg[0:fFrame-iFrame+1])


    axs[0].plot(times[0:fFrame-iFrame+1],eField_fluct_squared,label=r'$ \left(\epsilon_0 \langle|\delta {E}|^2\rangle_{x,y}/2\right)/ T_{e0} $')
    axs[0].plot(times[0:fFrame-iFrame+1],eField_fluct_squared_byT,label=r'$ \left(\epsilon_0 \langle|\delta {E}|^2\rangle_{x,y}/2\right)/ T_{e} $')
    axs[0].set_xlabel(r'$t \ [\omega_{pe}^{-1}]$', fontsize=30)
    axs[0].set_ylabel(r'$ \left(\epsilon_0 \langle|\delta {E}|^2\rangle_{x,y}/2\right)/ T_{e0} $', fontsize=30)
    axs[0].tick_params(labelsize = 26)
    axs[0].legend(fontsize = 18)

    axs[1].plot(times[0:fFrame-iFrame+1],eField_boxavg_z)
    axs[1].set_xlabel(r'$t \ [\omega_{pe}^{-1}]$', fontsize=30)
    axs[1].set_ylabel(r'$\langle E_z \rangle$', fontsize=30)
    axs[1].tick_params(labelsize = 26)
#    axs[1].set_ylim(top = 0.00)

    axs[2].plot(times[0:fFrame-iFrame+1],E_over_J_rolling)
    axs[2].set_xlabel(r'$t \ [\omega_{pe}^{-1}]$', fontsize=30)
    axs[2].set_ylabel(r'$\langle E_z\rangle /\langle\, J_z\rangle \ [\nu_{\mathrm{eff}}/ \omega_{pe}]$', fontsize=30)
    axs[2].tick_params(labelsize = 26)
#    axs[2].set_ylim(0.0, 2*np.amax(E_over_J_rolling))

    axs[3].plot(times[0:fFrame-iFrame+1],Telc_boxavg)
    axs[3].set_xlabel(r'$t \ [\omega_{pe}^{-1}]$', fontsize=30)
    axs[3].set_ylabel(r'$T_e /T_{e0}$', fontsize=30)
    axs[3].tick_params(labelsize = 26)
    axs[3].set_ylim(0.0, 2*np.amax(Telc_boxavg) )

    axs[4].plot(times[0:fFrame-iFrame+1],Tion_boxavg, label = r'$T_i /T_{e0}$')
    axs[4].plot(times[0:fFrame-iFrame+1],Tion_boxavg/Telc_boxavg, label = r'$T_i /T_{e}$')
    axs[4].set_xlabel(r'$t \ [\omega_{pe}^{-1}]$', fontsize=30)
    axs[4].set_ylabel(r'$T_i /T_{e0}$', fontsize=30)
    axs[4].tick_params(labelsize = 26)
    axs[4].set_ylim(0.0, 0.1)
    axs[4].legend(fontsize = 18)
#    axs[4].set_ylim(0.0, 2*np.amax(Tion_boxavg) )

    axs[5].plot(times[0:fFrame-iFrame+1],nu_eff_rolling, label = 'nu_eff')
    axs[5].plot(times[0:fFrame-iFrame+1],nu_Sagdeev, label = 'nu_Sagdeev')
    axs[5].set_xlabel(r'$t \ [\omega_{pe}^{-1}]$', fontsize=30)
    axs[5].set_ylabel(r'$\langle R_z\rangle /\langle\, J_z\rangle \ [\nu_{\mathrm{eff}}/ \omega_{pe}]$', fontsize=30)
    axs[5].tick_params(labelsize = 26)
    axs[5].legend(fontsize = 18)

    fig.tight_layout()
    plt.savefig(outDir+fileName+rf'_energy_current.png', bbox_inches='tight')
    plt.close()


if (auto_loading):

    for n in range(len(directories)):
        os.chdir(directories[n])

        input_files = []
        #loading up all input files
        for filename in os.listdir(os.getcwd()):
            if filename.endswith("~"):
                    os.remove(filename)
            elif '.lua' in filename:
                input_files.append(filename)
        #looping over all input files; selecting the fFrame and then plotting
        for i in range(len(input_files)):
            fileroot = input_files[i].split('.lua',1)[0]
            files_output_num = np.zeros(1)
            #looping over the output files to find the largest output time
            for filename in os.listdir(os.getcwd()):
                if (fileroot+'_elc') in filename:
                    output_number = re.findall(r'\d+.bp', filename)
                    if len(output_number)!=0:
                        output_number = int(output_number[0].split('.bp',1)[0])
                        files_output_num = np.append(files_output_num,output_number)
            max_output = int(np.max(files_output_num))
            if max_output==0:
                continue
            main_plotting(fileroot,max_output)
        os.chdir('..')
else:

    os.chdir(directories[0])
    #looping over all input files; selecting the fFrame and then plotting
    for i in range(len(fileNames)):
        fileroot         = fileNames[i].split('.lua',1)[0]
        files_output_num = np.zeros(1)
        #looping over the output files to find the largest output time
        for filename in os.listdir(os.getcwd()):
            if (fileroot+'_elc') in filename:
                output_number = re.findall(r'\d+.bp', filename)
                if len(output_number)!=0:
                    output_number = int(output_number[0].split('.bp',1)[0])
                    files_output_num = np.append(files_output_num,output_number)
        max_output = int(np.max(files_output_num))
        if max_output==0:
            continue
        main_plotting(fileroot, max_output)
