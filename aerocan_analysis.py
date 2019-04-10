from numpy import *
from scipy import signal
from datetime import datetime
from matplotlib import *
import pandas as pd
import matplotlib.pyplot as plt
from pysolar.solar import *
from dateutil import parser
import greg_jul
import statsmodels.api as sm

# some look-up constants etc
location='Egbert'
aerocan_date_start='20180718'
aerocan_date_finish='20180720'
aerocan_filepath='c:/Remote Sensing/Aerocan Data/'+aerocan_date_start+'_'+aerocan_date_finish+'_'+location+'.lev15'
aerocan_wavelength_list=['340','380','440','500']

# file path locator for Pandora data input
pandora_date='20180718'
pandora_serial='108'
pandora_filepath='c:/Remote Sensing/Pandora Data/Pandora'+pandora_serial+'s1_'+location+'_'+pandora_date+'_L1_smca1c2p1-5.txt'


for aerocan_wavelength in aerocan_wavelength_list:
    
    # file path locator for Pandora data output
    pandora_direct_sun_filepath='c:/Remote Sensing/Pandora Data/Pandora'+pandora_serial+'s1_'+location+'_'+pandora_date+'_L1_direct_sun_'+aerocan_wavelength+'nm.txt'
    pandora_direct_sun_aod_filepath='c:/Remote Sensing/Pandora Data/Pandora'+pandora_serial+'s1_'+location+'_'+pandora_date+'_L1_direct_sun_aod'+aerocan_wavelength+'nm.txt'
    figure_filepath='c:/Remote Sensing/Plots/'
        
    # a look-up table for photometer centre wavelengths and bandpasses
    if aerocan_wavelength=='340':
        centre_wavelength=339.7
        aerocan_bandpass=2.0
        et_spectrum=24.588
    elif aerocan_wavelength=='380':
        centre_wavelength=380.5
        aerocan_bandpass=4.0
        et_spectrum=24.785
    elif aerocan_wavelength=='440':
        centre_wavelength=439.9
        aerocan_bandpass=10.0
        et_spectrum=25.495
    elif aerocan_wavelength=='500':
        centre_wavelength=498.8
        aerocan_bandpass=10.0
        et_spectrum=25.569
        
    # do a pair of small conversions to use in Rayleigh calculations    
    centre_wavelength_um=centre_wavelength/1e3
    centre_wavelength_cm=centre_wavelength/1e7
    
    surface_count=2.546899e19 # N_s term from Bodhane, 1999, JAOT
    dry_refractive_index=1.00028702102066 # based on 360ppm CO2
    surface_pressure_cgs=9.94e5 # surface pressure in cgs units
    
    
    # import AOD from file into Pandas data frame
    aerocan_data=pd.read_csv(aerocan_filepath, header=6)#, delim_whitespace=True, header=2,skiprows=range(4,6))
    aerocan_data.columns=aerocan_data.columns.str.replace('(','_').str.replace(')','') # replace parentheses to remove read errors
    
        
    #aerocan_calibration_data=aerocan_data[(aerocan_data.day_of_year_fraction>=199.5) & (aerocan_data.day_of_year_fraction<=200.84)] # manual selection of time data for calibration
    aerocan_times=aerocan_data.Day_of_Year_Fraction
    aerocan_aod_column_header='AOD_'+aerocan_wavelength+'nm'
    aerocan_aod=aerocan_data.loc[:,aerocan_aod_column_header] # aerocan aod at selected wavelength
    aerocan_comparison_data=pd.concat([aerocan_times,aerocan_aod],axis=1)
    
    
    ###################### import full file pandora data from file by line to generate a pandas header index list #######################################
    with open(pandora_filepath) as f: # open the file 
        pandora_file_entries=f.readlines() # read all the lines in the file into an object
    
    column_headers_list=[] # initialize a column header list
    
    i=0 # initialize a line counter
    
    for entry in pandora_file_entries:
        
        lat_check=entry.find('Location latitude [deg]: ') # extract the latitude
        if lat_check>-1:
            lat=entry.split(': ')[1]
    
        lon_check=entry.find('Location longitude [deg]: ') # extract the longitude
        if lon_check>-1:
            lon=entry.split(': ')[1]
        
        altitude_check=entry.find('Location altitude [m]: ') # extract the altitude
        if altitude_check>-1:
            altitude=entry.split(': ')[1]
    
        
        wavelength_check=entry.find('Nominal wavelengths [nm]:') # fine the list of nominal wavelengths
        if wavelength_check>-1:
            wavelength_list=entry.split(': ')[1].split(' ') # populate a list of wavelengths
    
        column_check=entry.find('Column ') # look for the word column to start populating the column header list
        if column_check>-1:
            #print (entry.split(': ')[1].split(',')[0])
            column_headers_list.append(entry.split(': ')[1].split(',')[0]) # populate the coumn header list for first 80 columns
        data_check=entry.split(' ') # track the number of entries in a line
        if len(data_check)>3000: # look for the first line of data, which is more than 3000 entries wide
            break # do not bother reading any further
        skip_rows=i
        
        i+=1
        
    f.close() # close the file   

    print ('columns',len(column_headers_list))    
    
    # initialize and populate a set of header uncertainty lists
    L1_uncertainty_list=[] 
    L1_instrument_uncertainty_list=[]
    for wavelength in wavelength_list:
        L1_uncertainty_list.append(wavelength +'_nm_L1_Uncertainty')
        L1_instrument_uncertainty_list.append(wavelength +'_nm_L1_Instrument_Uncertainty')
    
    column_headers_list=column_headers_list+wavelength_list+L1_uncertainty_list+L1_instrument_uncertainty_list # concatenate all the lists to generate a master column header list
        
    
    # initialize an array of wavelengths matching the list of wavelengths headers
    wavelength_array=zeros(len(wavelength_list))
    i=0
    for wavelength in wavelength_list:
        wavelength_array[i]=float(wavelength)
        i+=1
    
    # locate the lower and upper Pandora wavelegnths based on the aerocan bandpass and aerocan centre wavelength
    lower_wavelength=centre_wavelength-aerocan_bandpass/2
    upper_wavelength=centre_wavelength+aerocan_bandpass/2
    #print (upper_wavelength)
    #print ('max wavelength',wavelength_list[searchsorted(wavelength_array,upper_wavelength)-1])
    
    
    # replace the array search floats with strings matching the wavelength headers
    lower_wavelength_header=wavelength_list[searchsorted(wavelength_array,lower_wavelength)].replace('.','_')
    upper_wavelength_header=wavelength_list[searchsorted(wavelength_array,upper_wavelength)-1].replace('.','_') # NOTE the -1 index adjustment
    
    
    # open the Level 1 data file as a Pandas data frame and extract the data centered around the wavelength of interest
    pandora_data=pd.read_csv(pandora_filepath,skiprows=range(skip_rows), names=column_headers_list, delim_whitespace=True) # open the file as a pandas dataframe
    pandora_data.columns=pandora_data.columns.str.replace(' ','_').str.replace('\n','').str.replace('.','_') # replace spaces and decimals with '_' and remove newline characters
    direct_sun_obs=pandora_data[pandora_data.Two_letter_code_of_measurement_routine=='SO'] # select direct sun measurements
    direct_sun_string_times=direct_sun_obs.loc[:,'UT_date_and_time_for_beginning_of_measurement']
    direct_sun_float_times=direct_sun_obs.loc[:,'Fractional_days_since_1-Jan-2000_midnight_for_beginning_of_measurement'] # select direct sun measurement float times
    direct_sun_integration_times=direct_sun_obs.loc[:,'Integration_time_[ms]'] # select direct sun measurement integration times
    direct_sun_counts=direct_sun_obs.loc[:,lower_wavelength_header:upper_wavelength_header] # select the 2um wavelength range centered at the photometer wavelength
    direct_sun_scale_factor=direct_sun_obs.loc[:,'Scale_factor_for_data']
    direct_sun_obs=pd.concat([direct_sun_float_times.to_frame(),direct_sun_integration_times,direct_sun_scale_factor],axis=1) # combine the times, integration times, scale factor into one dateframe
    
    # sum the counts across the wavelength range
    direct_sun_obs.loc[:,'mean_count_rate']=direct_sun_counts.mean(axis=1)
    direct_sun_obs.loc[:,'log_count_rate']=pd.Series(log(direct_sun_obs.mean_count_rate))
    
    ## create a normalized sum
    
    direct_sun_column_list=direct_sun_obs.columns.tolist()+['solar_zenith','rayleigh_air_mass_factor','aerosol_air_mass_factor']
    
    # create a dummy column to hold zenith angle and air mass factor
    direct_sun_obs=direct_sun_obs.reindex(columns=direct_sun_column_list)
    
    float_times=[] # initialize a blank float times list
    datetime_list=[]
    
    # calculate a solar elevation based on the UTC string time
    for index in direct_sun_obs.index:
        date_object=parser.parse(direct_sun_string_times[index]) # parse the string into a datetime object
        datetime_list.append(date_object)
        day_of_year=greg_jul.gregorian_to_julian(date_object.isoformat()[:10])[-3:] # use my simple date convertor to create a day of year string
        fraction_of_day=greg_jul.decimal_day(date_object.isoformat()[11:19]) # use my simple date convertor to create a fraction of day float
        float_time=float(day_of_year)+fraction_of_day # combine the day of year and fraction of day to create a time of year float
        float_times.append(str(format(float_time,'.6f'))) # append the float time to the float times list
        direct_sun_obs.solar_zenith[index]=90.-get_altitude(float(lat),float(lon),date_object) # calculate the solar zenith angle for the datetime object
    
    
    direct_sun_obs.loc[:,'day_of_year_fraction']=array(float_times) # create a dataframe column of float times
    direct_sun_obs.day_of_year_fraction=pd.to_numeric(direct_sun_obs.day_of_year_fraction) # convert datetime object to floats for merging
    direct_sun_obs.loc[:,'rayleigh_air_mass_factor']=pd.Series(1/(cos(arcsin((6371/(6371+6.2))*sin(radians(direct_sun_obs.loc[:,'solar_zenith'])))))) # create a dataframe column of rayleigh air mass factors from Kasten 1989
    direct_sun_obs.loc[:,'aerosol_air_mass_factor']=pd.Series(1/(cos(arcsin((6371./(6371+2.))*sin(radians(direct_sun_obs.loc[:,'solar_zenith'])))))) # create a dataframe column of aerosol air mass factors from Kasten 1989
#    direct_sun_obs.loc[:,'ozone_air_mass_factor']=pd.Series(1/(cos(arcsin((6371./(6371+20.4))*sin(radians(direct_sun_obs.loc[:,'solar_zenith'])))))) # create a dataframe column of aerosol air mass factors from Kasten 1989
        
    # calculate the unscaled log of the count rate
    direct_sun_obs.loc[:,'log_count_rate']=pd.Series(log(direct_sun_obs.loc[:,'mean_count_rate']/direct_sun_obs.loc[:,'Scale_factor_for_data']))
                                                                                                                           
#    # calculate the zenith equivalent count rate correcting for air_mass factor
#    direct_sun_obs.loc[:,'zenith_count_rate']=pd.Series(direct_sun_obs.loc[:,'mean_count_rate']*exp(direct_sun_obs.loc[:,'rayleigh_air_mass_factor']))
    
    # interpolate the aerocan aod onto pandora float times to match up pandora count rates with aerocan aod measurements
    interpolated_aerocan_aod=interp(float_times,aerocan_times,aerocan_aod)
    direct_sun_obs.loc[:,'coincident_aod']=interpolated_aerocan_aod
    direct_sun_obs.loc[:,'Date']=datetime_list
    direct_sun_obs.loc[:,'coincident_aod_with_air_mass_factor']=interpolated_aerocan_aod*direct_sun_obs.aerosol_air_mass_factor
    
    # calculate gravity for rayleigh optical depth calculation using List 1968
    geopotential_latitude=2*radians(float(lat)) # intemediate angle 
    latitude_cosine=cos(geopotential_latitude) # cosine term carried through
    g_0=980.616*(1-.0026373*latitude_cosine+5.9e-6*latitude_cosine**2) # sea-level gravity
    mass_weighted_altitude=0.73737*float(altitude)+5517.56
    g_z=g_0-(3.085462e-4+2.27e-7*latitude_cosine)*mass_weighted_altitude+(7.254e-11+1e-13*latitude_cosine)*mass_weighted_altitude**2-(1.517e-17+6e-20*latitude_cosine)*mass_weighted_altitude**3
    
    
    # calculate the scatter cross-section for dry air 
    nitrogen_depol=1.034+3.17e-4/centre_wavelength_um**2
    oxygen_depol=1.096+1.385e-3/centre_wavelength**2+1.448e-4/centre_wavelength_um**4
    king_factor=(78.084*nitrogen_depol+20.946*oxygen_depol+0.934+0.036*1.15)/(78.084+20.946+0.934+0.036)
    air_scatter=(24*math.pi**3*(dry_refractive_index**2-1)**2)/(centre_wavelength_cm**4*surface_count**2*(dry_refractive_index**2+2)**2)*king_factor
    
    # calculate the rayleigh optical depth from Bodhain et al. 1999
    rayleigh_optical_depth=6.0221367e23*air_scatter*surface_pressure_cgs/((28.9595+15.0556*3.6e-4)*g_z)
    
    print (rayleigh_optical_depth)
    
    # correct L1 count rate for aerosol optical depth and air mass factor
    direct_sun_obs.loc[:,'log_aerosol_corrected_count_rate']=pd.Series(direct_sun_obs.log_count_rate+direct_sun_obs.coincident_aod_with_air_mass_factor)
    
    
    ######### truncate the data to air mass factor <5
    truncated_direct_sun_obs=direct_sun_obs.loc[direct_sun_obs.rayleigh_air_mass_factor<5].copy()
    
   
    # determine mid-point of solar day
    
    # select only the morning data
    morning_direct_sun_obs=truncated_direct_sun_obs.loc[:truncated_direct_sun_obs.rayleigh_air_mass_factor.idxmin(),:].copy()
    
    print ('morning data',morning_direct_sun_obs.shape)
     
    # do a linear regression on the truncated date to start a first order filter
 #   linear_regression=sm.OLS(morning_direct_sun_obs.log_count_rate,sm.add_constant(morning_direct_sun_obs.rayleigh_air_mass_factor)).fit()
    linear_regression=sm.OLS(morning_direct_sun_obs.log_aerosol_corrected_count_rate,sm.add_constant(morning_direct_sun_obs.rayleigh_air_mass_factor)).fit()
    print ('lin regression params',linear_regression.params)
    et_count_rate=linear_regression.params['const'] # extract the LS intercept as the ET count rate
    air_mass_slope=linear_regression.params['rayleigh_air_mass_factor'] # extract the LS slope as the linear dependence of the count rate on air mass factor
    morning_direct_sun_obs.loc[:,'modeled_count_rate']=pd.Series(air_mass_slope*morning_direct_sun_obs.rayleigh_air_mass_factor+et_count_rate) # generate a linear model to emulate the count rate
    morning_direct_sun_obs.loc[:,'fractional_difference']=pd.Series(abs(linear_regression.resid)/morning_direct_sun_obs.modeled_count_rate*100) # compare the model to the actual count rate
    linear_regression_error=average(abs(linear_regression.resid))
    filtered_direct_sun_obs=morning_direct_sun_obs.loc[abs(linear_regression.resid)<linear_regression_error].copy() # filter the date with an arbitrary 2% cut-off
    
#    # filter the data recursively until the ET count rate converges
    print ('slope, intercept')
    print ('first run results ',air_mass_slope,et_count_rate)
    filtered_linear_regression=sm.OLS(filtered_direct_sun_obs.log_aerosol_corrected_count_rate,sm.add_constant(filtered_direct_sun_obs.rayleigh_air_mass_factor)).fit()
    filtered_et_count_rate=filtered_linear_regression.params['const']
    filtered_air_mass_slope=filtered_linear_regression.params['rayleigh_air_mass_factor']
    filtered_regression_error=average(abs(filtered_linear_regression.resid))
    et_count_rate_difference=abs(filtered_et_count_rate-et_count_rate)
    et_count_rate_uncertainty=filtered_linear_regression.bse['const']
    print ('filtered results ',filtered_air_mass_slope,filtered_et_count_rate)
    
    et_count_rate=filtered_et_count_rate
    filtered_direct_sun_obs=filtered_direct_sun_obs.loc[abs(filtered_linear_regression.resid)<filtered_regression_error] # filter the data using the residuals
    print('filtered data',filtered_direct_sun_obs.shape)

    average_aod=average(filtered_direct_sun_obs.loc[:,'coincident_aod'])

    print ('average_aod',average_aod)        
    
    ################calculate the aerosol optical depth from the Pandora data
    
    truncated_direct_sun_obs.loc[:,'pandora_aod']=pd.Series((et_spectrum-truncated_direct_sun_obs.log_count_rate-truncated_direct_sun_obs.rayleigh_air_mass_factor*rayleigh_optical_depth)/truncated_direct_sun_obs.aerosol_air_mass_factor)
    
    # calculate the difference
    
    truncated_direct_sun_obs.loc[:,'pandora_cimel_aod_difference']=pd.Series(truncated_direct_sun_obs.coincident_aod-truncated_direct_sun_obs.pandora_aod)
    difference_limit=average(truncated_direct_sun_obs.pandora_cimel_aod_difference)
    
    ######## output dataframe to file################
    truncated_direct_sun_obs.to_csv(pandora_direct_sun_filepath) # output the dataframe to csv file
    
    ############### plotting ####################
    
#    fig, ax = plt.subplots(figsize=(12,8 ))
#    fig.subplots_adjust(left=0.2, right=0.85)
#    
#    truncated_direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='pandora_aod', color='k', ylim=(0,1))#,logy=True)kind='scatter', 
#    truncated_direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='coincident_aod', secondary_y=False, color='r', ylim=(0,1))#,logy=True)kind='scatter', 
#    truncated_direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='pandora_cimel_aod_difference', secondary_y=True, color='b')
#    ax.set_xlabel('UTC Date & Hour')
#    ax.set_ylabel('AOD')
#    ax.right_ax.set_ylabel('Difference')
#    ax.right_ax.set_ylim(0,.1)
#    ax.legend([ax.get_lines()[0], ax.get_lines()[1],ax.right_ax.get_lines()[0]], ['Pandora AOD','Cimel AOD','Difference'], loc='lower left', frameon=False)
#    #ax.legend([ax.get_lines()[0], ax.get_lines()[1]], ['Pandora AOD','Cimel AOD'], loc='lower left', frameon=False)
#    plt.savefig(figure_filepath+'aod-comparison-'+aerocan_wavelength+'nm-'+pandora_date+'.png')
#    
#    plt.close(fig)
    
#    fig, ax = plt.subplots(figsize=(12,8 ))
#    fig.subplots_adjust(left=0.2, right=0.85)
#    
#    direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='log_count_rate', ylim=(0,14))#, ylim=(0,1))#,logy=True)kind='scatter', 
#    direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='coincident_aod', secondary_y=True, color='r')#,logy=True)kind='scatter', 
#    ax.set_xlabel('UTC Date & Hour')
#    ax.set_ylabel('ln(Count Rate (MHz))')
#    ax.right_ax.set_ylabel('Coincident AOD')
#    ax.legend([ax.get_lines()[0], ax.right_ax.get_lines()[0]], ['Rayleigh-Corrected Count Rate','Coincident AOD'], loc='lower right', frameon=False)
#    
#    plt.savefig(figure_filepath+'rayleigh-corrected-count-rate-'+aerocan_wavelength+'nm-'+pandora_date+'.png')
#    
#    plt.close(fig)
#    
#    
#    fig, ax = plt.subplots(figsize=(12,8 ))
#    fig.subplots_adjust(left=0.2, right=0.85)
#    
#    direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='log_rayleigh_corrected_zenith_count_rate', ylim=(0,14))#, ylim=(0,1))#,logy=True)kind='scatter', 
#    direct_sun_obs.plot(ax=ax, style=".", x='Date',  y='coincident_aod', secondary_y=True, color='r')#,logy=True)kind='scatter', 
#    ax.set_xlabel('UTC Date & Hour')
#    ax.set_ylabel('ln(Zenith Count Rate (MHz))')
#    ax.right_ax.set_ylabel('Coincident AOD')
#    ax.legend([ax.get_lines()[0], ax.right_ax.get_lines()[0]], ['Rayleigh-Corrected Zenith Count Rate','Coincident AOD'], loc='lower right', frameon=False)
#    
#    plt.savefig(figure_filepath+'rayleigh-corrected-zenith-count-rate-'+aerocan_wavelength+'nm-'+pandora_date+'.png')
#    
#    plt.close(fig)


    # plot the time series of the air mass factor
    
    fig, ax = plt.subplots(figsize=(12,8 ))
    fig.subplots_adjust(left=0.2, right=0.85)
    
    morning_direct_sun_obs.plot(ax=ax, style=".", x='rayleigh_air_mass_factor',  y='log_aerosol_corrected_count_rate', color='k')#,logy=True)kind='scatter', , secondary_y=True
    morning_direct_sun_obs.plot(ax=ax, style=".", x='rayleigh_air_mass_factor',  y='modeled_count_rate', color='b')#,logy=True)kind='scatter', 
    filtered_direct_sun_obs.plot(ax=ax, style=".", x='rayleigh_air_mass_factor',  y='log_count_rate', color='g')#,logy=True)kind='scatter', 
    morning_direct_sun_obs.plot(ax=ax, style=".", x='rayleigh_air_mass_factor',  y='fractional_difference', color='r', secondary_y=True)#,logy=True)kind='scatter', 
    ax.set_xlabel('Air Mass Factor')
    ax.set_ylabel('Log Count Rate')
    ax.right_ax.set_ylabel('% model difference', color='r')
    ax.legend([ax.get_lines()[0],ax.get_lines()[1],ax.get_lines()[2]], ['Count Rate', 'Modeled Count Rate','Filtered Count Rate'], loc='upper right', frameon=False)    
    ax.text(0.75,0.85,'ET Count Rate: '+str(filtered_et_count_rate)[:6],transform=ax.transAxes) 
    ax.text(0.75,0.8,'Slope (Total OD): '+str(-1.*filtered_air_mass_slope)[:6],transform=ax.transAxes)   
    ax.text(0.75,0.75,'Rayleigh OD: '+str(rayleigh_optical_depth)[:5],transform=ax.transAxes)
    ax.text(0.75,0.7,'Cimel AOD: '+str(average_aod)[:5],transform=ax.transAxes)
    ax.text(0.75,0.65,'Actual TOD: '+str(average_aod+rayleigh_optical_depth)[:5],transform=ax.transAxes)
    # alternate secondary y-axis method
    
    #ax.right_ax.set_ylabel('Modeled Count Rate')
    #ax.legend([ax.get_lines()[0], ax.right_ax.get_lines()[0]], ['OD Corrected Count Rate','Modeled Count Rate'], loc='lower left', frameon=False)
    
    plt.savefig(figure_filepath+'aod-corrected-count-rate-vs-air_mass-'+aerocan_wavelength+'nm-'+pandora_date+'.png')
    
    plt.close(fig)

#pandora_direct_sun_aod.to_csv(pandora_direct_sun_aod_filepath)

#pandora_data=pd.read_csv('c:/Remote Sensing/Pandora Data/Pandora108s1_Egbert_20180718_L1_smca1c2p1-5.txt',skiprows=range(108), names=column_headers_list, delim_whitespace=True, )
