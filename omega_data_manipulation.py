"""
OMEGA data manipulation.
Contained functions: omega_range_var_nsm, para_plot, addfile, omega_add_data, editfileAg
"""

#Import and reload all needed modules
import stellab
reload(stellab)
import omega
st = stellab.stellab()
import sygma
reload(sygma)
reload(omega)
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import math
import time as t_module


#Defining functions for name generation in plotting both coal time and power functions
def omega_range_var_nsm(c,e,n):
    """ omega_range_var_nsm(c,e,n). Input variable are to be in list format.
    c = Coalition time, e = Ejected mass, n = Neutron star mergers per solar mass
    Outpus as a string, must be converted to plot.
    """
    dummy=[]
    alpharange=[]
    crange=[]
    #Build up the naming convention used in the sim generation
    for i in range(len(c)):
        for j in range(len(e)):
            for k in range(len(n)):
                dummy='c'+str(c[i])+'e'+str(e[j])+'n'+str(n[k])
                alpharange.append(dummy)
    print alpharange
    return alpharange

#This is used when plotting the parameter plots
def para_plot(size,ej_plot,nsm_plot,nsm_def,nsm_ten,nsm_ms,title):
    """
    para_plot(size,ej_plot,nsm_plot,nsm_def,nsm_ten,nsm_ms,title).
    size       = size of plotted table
    ej_plot    = ejected mass plot range. List format
    nsm_plot   = neutron star merger mass density plot range. List format
    nsm_def    = NSM that is definitly within range. Embedded list format, expected length = nsm_plot
    nsm_ten    = NSM that is close to within range. Embedded list format, expected length = nsm_plot
    nsm_ms     = Inital nsm to be compared to the above nsm values, Likely same as nsm_plot, List format
    title      = Title to be given to plot
    
    """
    plt.figure(figsize=(5,5))
    for i in range(len(ej_plot)):#The ejected mass value range
        for j in range(len(nsm_plot)):#The nsm values
            if nsm_plot[j] in nsm_ten[i]:
                #Check to see if value is within Tentative bounds
                plt.scatter(x=ej_plot[i],y=nsm_plot[j],color='yellow',s=size, marker='s',label='Tentative')
            else:
                if nsm_plot[j] in nsm_def[i]:#Check to see if value is within a good fit's bounds
                    plt.scatter(x=ej_plot[i],y=nsm_plot[j],color='green',s=size, marker='s',label='Good')
                else:#Everything else is too far off
                    plt.scatter(x=ej_plot[i],y=nsm_plot[j],color='red',s=size, marker='s',label='Bad')
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.minorticks_on()
    plt.xlabel('Neutronstar mergers per solar mass (M$_{\odot}$$^{-1}$)'),plt.ylabel('Ejected mass (M$_{\odot}$)')
    plt.title('%s' %(title))
    plt.xlim((min(ej_plot)-0.005),(max(ej_plot)+0.005)), plt.ylim((min(nsm_ms)-0.5e-5),(max(nsm_ms)+0.5e-5))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.show()

def addfile(iniz,mgal,m_low,m_up,t,R,table):
    """
    addfile(iniz,mgal,m_low,m_up,t,R,table)
    iniz    = inital metallicity.
    mgal    = mass of Initial Gas Resevoir
    m_low   = Lower mass limit
    m_up    = Upper mass limit
    t       = dtd time values. list format
    R       = ratios for dtd at times t. list format. expect length = t
    table   = output table name.
    
    OUTPUTS:
    delayed_extra_dtd,delayed_extra_dtd_norm,delayed_extra_yields,delayed_extra_yields_norm
    """
    #These inputs are the Initial metalicity, std. mgal, Lower limit to mass range, Upper limit to mass range
    #t = ?, R=?, but likely the rate of nucleosynthsis, building upto a rate of 1 from 0 (0% to 100%) giving a build-up and cool-down
    #table is the input tanle you want to use
    print '############## Generating data from file %s ##############' %(table)
    # Run a SYGMA instance to access the Initial Mass Function (IMF)
    s_imf = sygma.sygma(iniZ=iniz, mgal=mgal)
    A = IMF = 1.0 / s_imf._imf(s_imf.imf_bdys[0], s_imf.imf_bdys[1], 2)
    nb_extra_star_per_m = A * s_imf._imf(m_low, m_up, 1)#SFR with IMF?
    # Create the DTD and yields information for the extra source
    # ==========================================================
    # Event rate [yr^-1] as a function of time [yr].
    # This assumes that all extra yields will be ejected
    # between 7.581E+06 and 9.588E+06 years (the lifetimes
    # of the 20 and 25 Msun models at Z = 0.02).    
    # Build the input DTD array
    dtd = []
    for i in range(0,len(t)):
        dtd.append([t[i], R[i]])

    # Add the DTD array in the delayed_extra_dtd array.
    delayed_extra_dtd = [[dtd]]
    
    # Define the total number of event per unit of Msun formed.  
    delayed_extra_dtd_norm = [[nb_extra_star_per_m]]

    # Define the total mass ejected by an extra source
    # Here, it would be best to find a correction factor 
    # to account for the different total mass ejected by
    # stars having different masses. For now, each star
    # in the mass range eject the same ejecta as the 20Msun
    # model.
    delayed_extra_yields_norm = [[1.0]]
    imf_yields_range = [m_low,m_up]
    # Define the yields path for the extra source
    extra_yields = [table]
    delayed_extra_yields = extra_yields
    print '####################################################################################' 
    return delayed_extra_dtd,delayed_extra_dtd_norm,delayed_extra_yields,delayed_extra_yields_norm #Outputs

#This is used as a quick way of generating omega plots for the above fn
def omega_add_data(delayed_extra_dtd,delayed_extra_dtd_norm,delayed_extra_yields,delayed_extra_yields_norm,mgal_in,coal,ej,nsm,sts):
    """
    omega_add_data(delayed_extra_dtd,delayed_extra_dtd_norm,delayed_extra_yields,delayed_extra_yields_norm,mgal_in,coal,ej,nsm,sts)
    delayed_extra_dtd                  = dtd generated by addfile
    delayed_extra_dtd_norm             = dtd generated by addfile
    delayed_extra_yields               = dtd generated by addfile
    delayed_extra_yields_norm          = dtd generated by addfile
    mgal_in                            = input mgal
    coal                               = coalition time
    ej                                 = ejected mass
    nsm                                = nsm input
    sts                                = special time steps    
    """
    run = omega.omega(galaxy='milky_way', mgal=mgal_in, ns_merger_on=True, special_timesteps=sts,\
                        t_nsm_coal=coal,m_ej_nsm=ej,nb_nsm_per_m=nsm,\
                        delayed_extra_dtd=delayed_extra_dtd,\
                        delayed_extra_dtd_norm=delayed_extra_dtd_norm, delayed_extra_yields=delayed_extra_yields,\
                        delayed_extra_yields_norm=delayed_extra_yields_norm, transitionmass=8.0)
    return run

#This takes data from ndw and strips all but the silver from it, it then changes it into Msun via ejected mass
def editfileAg(inputfile,outputfile,ej_m_mod,metalicity):
    """
    editfileAg(inputfile,outputfile,ej_m_mod,metalicity)
    inputfile  = Input file
    outputfile = name of outgoing file, will edit if there or create where needed
    ej_m_mod   = Change in the ejected mass
    metalicity = Metallicity used
    """
    fread = open(inputfile,'r')
    fprint = open(outputfile,'w')
    line = fread.readlines()
    del line[0]
    element=[]
    yields=[]
    for x in line:
        element.append(x.split('&')[1])
        yields.append(x.split('&')[2])
    fprint.write('H Neutrino-driven wind yields, from:%s  Using ejected mass:%s\n' %(inputfile,ej_m_mod))
    fprint.write('&Isotopes  &Z=%s\n' %metalicity)
    for i in range(len(element)):
        if element[i]=='Ag-107   ':
            fprint.write('&%s&%s\n' %(element[i],ej_m_mod*(float(yields[i]))))
        if element[i]=='Ag-109   ':
            fprint.write('&%s&%s\n' %(element[i],ej_m_mod*(float(yields[i]))))
    fread.close()
    fprint.close()
    print 'File created with name%s' %(outputfile)