from library_functions_FPZTF import *


#TODO: make this class ONLY a provider of correct photometry. let the epxlosion time and the peak magnitude be satellite functions or classes ? > YUP ON IT 
#TODO: adapt to the new filescheme what new file scheme??? 




###############################################

#               CLASS                         #

###############################################



class FP_ZTF( object ): 
    '''
    This class defines the object "photometric table" from the forced photometry table. It takes as input the forced
    photometry file and the filter in which you want to obtain the "final" product of the photometric table. 
    For more information, check the pdf file with al the information on forced photometry
    
    '''
    



    #--------#
    # SETTER #
    #--------#

    
    
    def __init__(self, forced_phot_filename , filtre  , redshift, RA, DEC, FD , DM = None  ):
        '''
        This function initisalises the class
        
        Parameters
        ----------
        self
        forced_phot_file [path+string] path + name of the file which we have to be treated
        filtre [string] : each object will be discriminated based on the ztf filter of the forced photometry
        name   [string] : object name
        plots  [bool]-optional- Whether you want to visualize some plots like the raw LC plot and the flux vs chi2 for the uncertainties. 
        other things can be added afterwards
        -optional- codes of processing status you want removed form the table
                                          ['255','64','56','61'] are the default status filtered out  
        
        Returns
        -------
        The actual things you're going to be working with goddamit
        
        '''
        self.cosmo = cosmology.Planck13 # We are using the Planck 13 cosmology for our calculations 
        
        self._filter_name   = {'ZTF_g':'darkolivegreen','ZTF_r':'orangered', 'ZTF_i':'gold'}
        
        self.set_metadata(redshift, RA, DEC, FD, DM)
        
        self.set_rawtable(forced_phot_filename)
        self.set_filtre(filtre)
        
        self.set_table()
    
        self._corrected_LC     = False 
        self._corrected_4_z    = False
        self._corrected_4_ext  = False
        self._converted_to_mag = False
        self._normalised_LC    = False
        self._corrected_unc    = False
        self._interpolatedLC   = False
        self._binned_LC        = False
        
        
    
        
        
    def set_rawtable(self, forced_phot_filename):
        '''
        
        '''
        self._rawtable = read_table_forced_phot(forced_phot_filename)
    
    
    
        
       
        
    def set_metadata(self, redshift, RA, DEC, FD , DM ):
       

        '''
        Sets the metadata such as redshift and coordinates which will be use afterwards 
        parameters
        ----------
        redshift    [float]  self explanatory
        RA          [float]  Right ascension (degrees) [0,180]
        DEC         [float]  Declinaisation (degrees)  [-90,90]
        FD          [float]  JD of first detection


        '''
   
        self._redshift     = redshift

        # Which cosmology should we assume? Plack 13.  Why? Cause 18 is fucked up
        if DM is not None: 
            self._DM = DM
        else:
            self._DM           = self.cosmo.distmod(z = self._redshift).value
        
        ### TODO: remmber in future iteration to add the error on the redshift, right now it's not necessary
        # self._e_DM         = (5*_table_infoobj['e_redshift'])/(np.log(10)*_table_infoobj['redshift']) # for small z!! 
        
        self._coordinates  = SkyCoord(np.array(RA), np.array(DEC),unit = 'deg', frame = 'fk5')
        
        self._jd_fd  = FD # This is given by the full table. 
      
              
    
    def set_filtre(self, filtre):
        '''
        sets the filter in which you want to obtain the photometry. the variable is only an entry variable. 
        parameter
        ---------
        filtre [string] name of the filter/ 'ZTF_{g,i,r}'
        
        '''
        self.filtre = filtre
        
        
    def set_table(self):
        '''
        
        This function prepares the rawtable to be worked with in a single band (filter). 
        It removes all the 'null' measurements. 

        #TODO: In a future version, give the procstatus on which to filter. 
        In this version, we implement quality checks on the seeing, bad pixels 
        
        Parameters 
        ----------
        self 
        
        Returns
        -------
        self.table
        '''
        # Selects only measurements in the relevant filters
        self.table = self._rawtable[self._rawtable['filter'] == self.filtre]
        
        # removing the "null" measurements : this is when the fit failed
        self.table = remove_null_measurements(self.table)


        #NOTE: sometimes the the measurements are not null, but the nearest reference source within 5" returns nothing so it still gives a null. 
        # Make sure to remove the null flux measurements before converting to nans and flaots
        if len(self.table['jd']) > 0 : 
            #removing the actual strings from the column which will be of no use for this object
            self.table.remove_columns(['filter']) # string removal, each light curbe is per filter

            self.table = convert_table_float(self.table)

            # Giving each observation a unique ID based on the CCD ID, the Quadrant and the Field IDs. 
            self.table = create_new_id(self.table)

            self.table['tfromFD'] = self.table['jd'] - self._jd_fd
            
            
        elif len(self.table['jd']) == 0 : 
            
            print(f'There does not seem to be any data {self.filtre}')
        

          
    #--------#
    # FILTER #
    #--------#

    def remove_big_errors(self, thres = 7):
        '''
        #IMPORTANT: this function removes a significant amount of data points? THis should be better investigated... 
       
        This function looks for the median of the error of the flux measurements and removes those with errors bigger than 3*median
        
        parameters
        ----------
        thres [float] the threshold of the MAD we remove 


        returns
        -------
        table without the faulty data points

        '''

        self.table = remove_big_err_point(self.table, threshold  = thres)
        # return remove_big_err_point(self.table, threshold  = thres)

    
    
    def reset_table(self):
        '''
        Useful function to reinitialise the table after applying cuts. 

        Warning
        -------
        Remember to correct the lightcurve and normalise it before using the fitter
        
        '''
       
        del self.table
        self.set_table()
        
        self._corrected_LC      = False 
        self._corrected_4_z     = False
        self._corrected_4_ext   = False
        self._normalised_LC     = False
        self._converted_to_mag  = False
        self._corrected_unc     = False
        self._binned_LC         = False
        self._clean_LC          = False
        self._added_jd_from_exp = False


        
    #--------#
    # GETTER #
    #--------#
    
    def get_baseline(self, binned = False):
        '''
        this functions aims at estimating the baseline of the non detections 
        Baseline correction is important vefore conversion to magnitude / upper limit. 
        Parameters
        ----------
        indiex [int]
        
        
        Returns
        -------
        Value of the median of the baseline and its standard deviation
        '''
        if self._corrected_LC  is True: 
            if binned is False:
                _temp = self.table['extcorrforcediffimflux'][self.table['tfromFD'] < -2.5]  
                return np.median(_temp ), np.std( _temp )

            elif binned is True:
                if self._binned_LC is True:
                    _temp = self.binned_lc['extcorrforcediffimflux'][self.binned_lc['tfromFD'] <= -2.5] 
                    return np.median( _temp ), np.std( _temp )
                else : 
                    print('You need to bin the LC first')
        else:
            print('You need to correct the LC first')

    
    
   
   
    
    #-----------#
    # Corrector #
    #-----------#
    
    def correct_lc(self, correct_unc = True , correct_z = False , correct_ext = True, correct_baseline = True, SNR_thres = 3 ,
                        correct_baseline_method = 'separated'):
        '''
        This function corrects the lightcurve for:
        0) remove known bad measurements
        1) recalibrating the uncertainties on the PSF flux measurement as suggested in the FP guide 
        2) correct for extinciton (by default)
        3) can correct for redshift (but this is turned off by default. ) #NOTE: this will be relevant when working with GP for LC parameters

        It can also "correct" the uncertainty which by default is true here. This only means multiplying the raw errors by the sqrt of the chi2
        Correcting the light curve is physical and is based on the measurements and on the redshift, zero point and etc 


        (il faut faire les quality cuts bien avant les corrections)
        
        Parameters
        ----------
        self
        correct_unc       [bool] whether you want to also correct the uncertainties. By default it will be true.  
        correct_z         [bool] whether to correct for redshift
        correct_ext       [bool] whether to correct for extinction
        correct_baseline  [bool] whether to correct for the baseline offset. By default True. 
        SNR_thres         [int] the S/N threshold to filter the SN LC from noisy points (not the baseline)
        
        
        returns
        -------
        corrected light curve (self table)
        
        
        '''
        ############# STEP 1: filter known bad measurement: 
        # Remove measurements with a seeing > 4" 
        self.table = filter_bad_seeing(self.table, seeing = 4)


        # remove high spatial nosie sigma per pixel , by default threshold is 25. Can be stricter.
        self.cuts_scisigpix()


        ############ STEP 2 : correct light curve of physical parameters

        # baseline correction (should be done before the zero point correction, no?)
        # important: should correct for the offset BEFORE cleaning the baseline... !!! 
        
        if correct_baseline is True:
            self.correct_baseline_offset(correction_method=correct_baseline_method)
            self.clean_baseline()
            
            

        # correction for zeropoint

        self.correct_zeropoint()


        # recalibration of uncertainties of diff flux
        
        if correct_unc is True:
            self.correct_uncertainties()
        self.compute_errors_on_final_diffflux()

        
        # Remove low S/N (3) in the Supernova LC (outside of the baseline)
    
        self.filter_on_SNR(SNR_thres=SNR_thres)


        # remove big big errors , there's a bug in this thing... 

        self.remove_big_errors()


        # correction for extinction and redshift

        self.correct_redshift(correct_z=correct_z)
        self.correct_extinction(correct_ext=correct_ext)

        # if add_jd_from_FD is True: 
        #     self.add_t_from_FD()



        # toggles this variable to True to let the user know that corrections were made.
        self._corrected_LC = True


    # def add_t_from_FD(self):
    #     '''
    #     If the table we provide has the explosion time already computed, this function adds a column of time from explosion

    #     parameters
    #     ----------

    #     returns
    #     -------

        
    #     '''

    #     if self._corrected_4_z is True: 
    #         self._jd_FD_zcor          = self.FD/(1+self._redshift)
    #         self.table['tfromFD_zc'] = self.table['jd_zcorr'] - self._jd_texp_zcor
    #         # print(self.table)
    #         self._added_jd_from_exp     = True
    #     else:
    #         self.table['tfromFD_znc'] = self.table['jd'] - self._jd_texp
    #         self._added_jd_from_exp = True
        


    ##################
    #    MAG LC      #
    ##################
    

    def add_magnitudes_detvndet_meas(self, return_magtable = False):
        '''
        This function computes the detections vs Non detections. 
        According to the Forced photometry guide, we choose SNT = 3 and SNU = 5. 
        If a the flux divided by the uncertainty is bigger/equal than 3, it is a dectection. If flux/unc < 3, it's a non-detection.
        Limiting magnitudeds are computed as 5 times sigma.

        parameters
        ----------
        return_magtable [bool] whether you want to return only the magnitude table. it will still be added to the main table

        returns
        -------

        note
        ----

        When computing the flux, we already account for the zero point so we don't need to take it into account again.

        '''

        if self._converted_to_mag == False:

            SNT_ = 3
            SNU_ = 5

            _magtable = Table(names = ('mag','emag','limmag', 'absmag', 'e_absmag'))

            for indy in range(len(self.table['jd'])):

                if self.table['extcorrforcediffimflux'][indy]/self.table['extcorrforcediffimfluxunc'][indy] >= SNT_: 

                    mag_          = to_mag(self.table['extcorrforcediffimflux'][indy])
                    absmag_       =  mag_ - self._DM
                    emag_         = 1.0857 * self.table['extcorrforcediffimfluxunc'][indy] / self.table['extcorrforcediffimflux'][indy]
                    # e_absmag_     = np.sqrt(emag_ **2 + self._e_DM**2)
                    e_absmag_     = np.sqrt(emag_ **2 ) ### VERSION WITHOUT THE ERROR ON THE REDSHIFT 
                    limmag_       = -2.5   * np.log10(SNU_*abs(self.table['extcorrforcediffimfluxunc'][indy]))
                else:
                    mag_          = 99.
                    absmag_       = 99.
                    emag_         = 99.
                    e_absmag_     = 99.
                    limmag_       = -2.5*np.log10(SNU_*abs(self.table['extcorrforcediffimfluxunc'][indy]))



                _magtable.add_row([  mag_ , emag_ , limmag_ , absmag_ , e_absmag_  ]  )

            self.table = hstack([self.table, _magtable ])      

            self._converted_to_mag = True
            

            if return_magtable is True:
                _magtable = hstack([self.table['tfromFD'], _magtable ])
                # if self._added_jd_from_exp is True:
                #     if self._corrected_4_z is True:
                #         _magtable = hstack([self.table['tfromexplo_zc'], _magtable ])    
                          
                #         return _magtable
                #     else:
                #         _magtable = hstack([self.table['tfromexplo_znc','jd'], _magtable ])      
                #         return _magtable

        else:
            print('The table was already converted to magnitudes. No magtable returned')
            return None




    def clean_baseline(self):
        '''
        this function computes the median value of the flux and the median absolute deviation in the baseline (i.e. <-2.5d from FD). 
        It removes any point more than 3 MAD away from the median flux

        NOTE: Clean the baseline BEFORE you correct for the baseline

        '''
        self.table = clean_baseline_phot(self.table)
    
    
    def correct_baseline_offset(self, correction_method = 'separated'):
        '''
        This function computes the median of the baseline and substracts this value to the full flux light curve.


        parameters
        ----------
        correction_method : choose 'separated' to correct the baseline according to the unique OBS ID (field/quadrant/ccd) to avoid any bias. Otherwise, choose 'all' to not distinguish between
        different field/quadrant/ccds. 


        '''
        if correction_method == 'separated':
            self.table = correct_4_baseline_offsets(self.table)

        elif correction_method == 'all':
            self.table = correct_4_baseline_offset_withoutsep(self.table)

        




    def normalise_flux(self):
        '''
        this functions normalises the flux to 1 in order to avoir very small numbers. 
        
        We normalise to the maximal flux point

        /!\ The flux which is being normalised is the extinction AND zero point corrected flux.
        parameters
        ----------
        photfluxe

        returns
        -------
        normalised photflux
        '''
        _normabite                                    = max(self.table['extcorrforcediffimflux'])
        self.table['extcorrforcediffimflux_norm']     = self.table['extcorrforcediffimflux']/_normabite
        self.table['extcorrforcediffimfluxunc_norm']  = self.table['extcorrforcediffimfluxunc']/_normabite

        self._normalised_LC = True
        


        
    def correct_uncertainties(self):
        '''
        This function re estimates the uncertainty on the forced flux measurement fits according to a chi_2 dependency
        As recommended by the forced photometry guide. 
        # Cprrect it before the zero point conversion? 
        
        parameters
        ----------
        self
        
        returns
        -------
        self.table['forcediffimfluxunc_chcorr'] corrected
        
        '''
        # for _ in range(len(self.table['forcediffimfluxunc'])):
        #     if self.table['forcediffimchisq'][_] >= 1.15: # completly arbitrary, let's be honest
        #         self.table['forcediffimfluxunc'][_] = self.table['forcediffimfluxunc'][_] * np.sqrt(self.table['forcediffimchisq'][_])
        
        self.table['forcediffimfluxunc'] = self.table['forcediffimfluxunc'] * np.sqrt(self.table['forcediffimchisq'])

        self._corrected_unc = True
    

    def compute_errors_on_final_diffflux(self):
        '''
      This function computes the errors on the computed Flux of the difference image (z_p corrected). Should also be 
      corrected for extinction, redshift and zeropoint
      Final means that these are the errors on the zp corrected flux... 
      
      parameters
      ----------
      the entire table, in the future this funciton will be in the class so it will be easier to handle
       
      returns
      -------
      
    
    
      '''
        _f_0 = 10**((-0.4) * self.table['zpdiff'])
        
        
        self.table['extcorrforcediffimfluxunc'] = np.sqrt( ( _f_0 *  self.table['forcediffimfluxunc'] )**2 + 
                                                            ( (-0.4) * np.log(10)* self.table['forcediffimflux'] * _f_0 * self.table['zpmaginpsciunc'] )**2  )  
            
    
    def filter_on_SNR(self,SNR_thres): 
        '''
        This function filters on the S/N ration of the diff flux

        parameters
        ----------
        S/R threshold
        '''

        _temp_SN       = self.table[self.table['tfromFD'] >= 0.] #to not kick out important first few points
        _temp_Baseline = self.table[self.table['tfromFD'] < 0.]

        _SN       = filter_StN_ratio(_temp_SN,SNR_thres)
        # _Baseline = filter_StN_ratio(_temp_Baseline,0.5)

        self.table = vstack([_temp_Baseline,_SN])

    
    def correct_zeropoint(self):
        '''
        This functions corrects for the zero point

        The zero point of an instrument, by definition, is the magnitude of an object that produces one count (or data number, DN) per second. 
        The magnitude of an arbitrary object producing DN counts in an observation of length EXPTIME is therefore:
        m = -2.5 x log10(DN / EXPTIME) + ZEROPOINT 

        For some obscur reason, the zero point converted to flux is already in DN ... So no need to divide by the exposure time. 
        
        parameters
        ----------
        table
    
        returns
        -------
        self.table[''] new column
    
        '''
    
        #return 10**((-0.4) * (-2.5*np.log10((table['forcediffimflux'])/table['exptime']) + table['zpdiff'] ) ) -> doesn't work???

        self.table['zpcorrfordifflux'] = (self.table['forcediffimflux']) * 10**((-0.4) * self.table['zpdiff'] ) 


    def correct_redshift(self, correct_z = False):
        '''
        This function creates a new column with the time corrected for redshift. We here do not consider correction on the flux because of redshift 
        
        parameters
        ----------
        table
        correct_z : bool. meaning we are considering only measured parameters and are not correcting for redshift. 
        
        NOTE: there's not real interest in correcting for redshift here? also not sustainable is no error on the redshift?? 
        
        returns
        -------
        
        '''
        
        
        # about correcting the flux with redsdhift, eran says we should, but ?
        if correct_z is True:
            # self.table['zp_zcorrfordifflux'] = self.table['zpcorrfordifflux'] * (1+self._redshift)
            self.table['zp_zcorrfordifflux'] = self.table['zpcorrfordifflux'] * (1+self._redshift) # correcting the flux for redshift effect ?

            self.table['jd_zcorr'] = self.table['jd']/ (1+self._redshift)

            self.table['tfromFD'] = (self.table['tfromFD']) / (1+self._redshift) 

            self._corrected_4_z = True
        else:
            self.table['zp_zcorrfordifflux'] = self.table['zpcorrfordifflux']
            self.table['jd_zcorr'] = self.table['jd']/ (1+self._redshift)
            

        
        # time from marshal first detection corrected for redshift
        
        
        
        
    def correct_extinction(self,correct_ext = True):
        '''
        ISM correction... 
        WARNING: make sure you have the dust maps and that they are stored in the right folder
        YOU NEED TO CORRECT FOR ZEROPOINT BEFORE USING THIS FUNCTION
        '''
        
        
        _p48   = { 'ZTF_r': 6339.6,'ZTF_g': 4722.7,'ZTF_i': 7886.1 }
        
        
        _gal_reddening = sfdmap.SFDMap('/Users/r.olivaw/Dropbox (Weizmann Institute)/ASTRO/code_library/python/sfddata-master/', scaling=0.86)
        _ebv           = _gal_reddening.ebv(self._coordinates)
        _alam          = extinction.ccm89(np.array([_p48.get(self.filtre)]) , _ebv * 3.1, 3.1, unit = 'aa')
        
        if correct_ext is True:
            self.table['extcorrforcediffimflux'] = self.table['zp_zcorrfordifflux']*(10**(0.4*_alam))
            self._corrected_4_ext = True

        else:
            self.table['extcorrforcediffimflux'] = self.table['zp_zcorrfordifflux']
        
        
        
        
        
    
    def cuts_scisigpix(self, thres_sci = None):
        '''
        We defined an abcolute value of 25 (Yuhan), twice what 

        Robust sigma per pixel in sci image [DN]

        parameters
        ---------
        thres_sci  [float] How many MAD away from median, by default = 5
        
        returns
        -------
        Does not return anything, just filters the table based on the scisigpix threshold 
        '''

        # faut checked en plus ou moins non? ...
        if type(thres_sci)==float:
            _thres_sci = np.median(self.table['scisigpix']) + thres_sci * median_abs_deviation(self.table['scisigpix']) #actually more strict than what yuhan does... 
            self.table = self.table[self.table['scisigpix'] <= _thres_sci ]
        else:
            self.table = self.table[self.table['scisigpix'] <= 25 ]
       



 
    #---------#
    # PLOTTER #
    #---------#
    
    def plot_lc(self, save_fig = False, path_save = None ,
                lc_corr = False, norma = False, add_plot = False, **kwargs):
        '''
        this function plots the light curve
        
        Parameters
        ----------
        add_plot [bool] if you want to plot several LC on the same figure, toggle to "True". Otherwise per default it opens a new figure per LC.
        
        Returns
        -------
        '''
        if add_plot is False:
            plt.figure()
        
        if lc_corr is False :
            plt.errorbar(self.table['jd'], self.table['forcediffimflux'], self.table['forcediffimfluxunc'], fmt='o',  ms = 2.5, color = self._filter_name.get(self.filtre),**kwargs)
            plt.xlabel('JD')
            plt.ylabel('#Counts')

        elif lc_corr is True and norma is False: 
            if self._corrected_LC is True :
                plt.errorbar(self.table['tfromFD'], self.table['extcorrforcediffimflux'], self.table['extcorrforcediffimfluxunc'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre),**kwargs)
                plt.xlabel('Time from first ASFD', size = 15)
                plt.ylabel('Flux [Jy]', size = 15)

            else :
                raise Exception('You need to correct the LC before plotting it')
        elif norma is True: 

            if self._corrected_LC is True and self._normalised_LC is True : 
                plt.errorbar(self.table['jd_zcorr'], self.table['extcorrforcediffimflux_norm'], self.table['extcorrforcediffimfluxunc_norm'], fmt='o',  ms = 2.5,color = self._filter_name.get(self.filtre),**kwargs)
                plt.xlabel('JD, z corrected')
                plt.ylabel('Jy, norm')
            else :
                raise Exception('You need to normalise the light curve first')
        


        
        
        if save_fig == True :
            if path_save != None:
                plt.savefig(path_save)
                
            else :
                print('You need to provide a path to save the figure')




    def plot_maglc(self, save_fig = False, path_save = None ,
                add_plot = False, **kwargs):
        '''
        this function plots the magnitude light curve
        
        Parameters
        ----------
        add_plot [bool] if you want to plot several LC on the same figure, toggle to "True". Otherwise per default it opens a new figure per LC.
        
        Returns
        -------
        '''
         ##### FUNCTION FOR THE SECONDARY AXIS 
        def to_abs_mag(x,dm = self._DM):
            '''
            converts to absolute magnitude knowing the distance modulus
            '''
            return x - dm

        def to_app_mag(x,dm = self._DM):
            '''
            converts to apparent magnitude knowing the distance modulus
            '''
            return x + dm



        if add_plot is False:
            plt.figure()
            ax = plt.subplot(111)
 
        if self._converted_to_mag==True:

            _det  = self.table[self.table['mag']!=99.]
            _ndet = self.table[self.table['mag']==99.]


            ax.errorbar(_det['tfromFD'], _det['mag'], _det['e_absmag'], fmt='o', alpha = 0.2 ,ms = 3.5, elinewidth=4, color = self._filter_name.get(self.filtre),**kwargs)
            ax.errorbar(_det['tfromFD'], _det['mag'], _det['emag'], fmt='o',  ms = 3.5, color = self._filter_name.get(self.filtre),**kwargs)
            
            ax.plot(_ndet['tfromFD'], _ndet['limmag'], lw=0 ,marker = 'v', ms = 3.5, color = self._filter_name.get(self.filtre),**kwargs)
            

            ax2 = ax.secondary_yaxis('right', functions=(to_abs_mag, to_app_mag))

        
            ax.set_ylabel('Apparent Magnitude', size = 15)
            ax2.set_ylabel('Absolute Magnitude', size = 15)
            plt.axvline(0,alpha = 0.3, ls = ':', lw = 0.75)
            plt.xlabel('Days from A.S. first detection', size = 15)
            # plt.xlim([-50,])
            plt.xlim(left=-3)
            # plt.ylim([16.5,21.5])

            plt.gca().invert_yaxis() 

            if save_fig == True :
                if path_save != None:
                    plt.savefig(path_save)
                else :
                    print('You need to provide a path to save the figure')

        else: 
            print('You need to add the magnitudes column before plotting it. See the function add_magnitudes_detvndet_meas. ')



    








