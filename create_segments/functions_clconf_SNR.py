from __future__ import print_function
import numpy as np 
import sys 
import h5py 
from ConfigParser import ConfigParser
import os
import gwpy 
import scipy 
import math 
import random 
import time


from scipy import stats as sc #for d'agostino-pearson test

# Matched filter 
from pycbc.filter import matched_filter

from gwpy.timeseries import TimeSeries 
from gwpy.frequencyseries import FrequencySeries
from gwpy.segments import Segment, SegmentList

from pycbc.waveform import get_td_waveform
from pycbc.detector import Detector


def create_config(config_file):
    config = ConfigParser()
    config.read(config_file)
    return config 


#config = create_config('config_main_denoising_usnr.ini')
config = create_config(sys.argv[2])


# General settings 
detector   = config.get('conf', 'detector')
glitchfile = config.get('conf', 'glitchfile')
segmenfile = config.get('conf', 'segmenfile')
segmentdir = config.get('conf', 'segmentdir')

#check if a file with injection exists
injfile = 0
try:
        injfile    = config.get('conf', 'injfile')
        print('injction file',injfile)
except:
        print("ATTENTION: injection file not present in the configuration file")        

#check if an additional list of glitch exists        
glitchfile_gspy = 0
try:
        glitchfile_gspy    = config.get('conf', 'glitchfile_gspy')
        print('additional glitch file',glitchfile_gspy)
except:
        print("ATTENTION: additional glitch file not present in the configuration file")

#sigfile    = config.get('conf', 'sigfile')

# threshold on the snr of the glitches to be considered for the class "glitch"
# default value 0
snr_th = 0
try:
        snr_th = float(config.get('conf', 'snr_th'))
        print("Threshold in SNR",snr_th)
except:
        print("ATTENTION: you have not defined a threshold for the SNR of the glitches")

glitch_cat = None
try:
        glitch_cat = config.get('conf', 'glitch_cat')
        print("Selecting glitches of categlory",glitch_cat)
except:
        print("ATTENTION: you have not selected a glitch category")


noiseoutfile  = config.get('conf', 'noiseoutfile')
glitchoutfile = config.get('conf', 'glitchoutfile')
sigoutfile    = config.get('conf', 'sigoutfile')

fs         = int(config.get('conf', 'fs'))
crop       = float(config.get('conf', 'crop'))
duration   = float(config.get('conf', 'duration'))
highpass   = float(config.get('conf', 'highpass'))
fduration  = float(config.get('conf', 'fduration'))
fftlength  = float(config.get('conf', 'fftlength'))


# waveform specific config 
apx = config.get('wf', 'apx')
f_lower = config.get('wf', 'f_lower')

# limits on the snr of the signal 
snr_limit_min = float(config.get('wf', 'snr_limit_min')) 
snr_limit_max = float(config.get('wf', 'snr_limit_max')) 

# ranges of masses, spins and distance  
m1 = [float(x) for x in config.get('wf', 'm1range').split(",")]
m2 = [float(x) for x in config.get('wf', 'm2range').split(",")]
s1 = [float(x) for x in config.get('wf', 's1range').split(",")]
s2 = [float(x) for x in config.get('wf', 's2range').split(",")]
di = [float(x) for x in config.get('wf', 'distance').split(",")]



#critical value for the choice of the noise, default value 1e-1
pcrit = 1e-1
try:
        pcrit = float(config.get('conf', 'pcrit'))
        print("pcrit from config file:",pcrit)
except:
        print("default value of pcrit:",pcrit)
        

#mb !!! obsolete now !!! 
def random_signal_waveform(detector, gps, apx, m1, m2, s1, s2, di): 
    
    det = Detector(detector)

    #mb uniform distribution - OK for testing? 
    # random choice of parameters for masses, spins, inclination, 
    # phase at coalescence, distance    
    mass1 = np.random.uniform(m1[0], m1[1]) 
    mass2 = np.random.uniform(m2[0], m2[1])

    spin1z = np.random.uniform(s1[0], s1[1]) 
    spin2z = np.random.uniform(s2[0], s2[1])

    inclination = np.random.uniform(0, np.pi) 
    coa_phase = np.random.uniform(0, 2*np.pi)

    distance = np.random.uniform(di[0], di[1])

    #mb this generates the waveform 
    #!!! f_lower is set to 40 [Hz] for testing (config_main.ini,  
    # to make the time-domain waveform not too long) 
    hp, hc = get_td_waveform(approximant=apx,
            mass1=mass1, 
            mass2=mass2,
            spin1z=spin1z, 
            spin2z=spin2z,
            inclination=inclination, 
            coa_phase=coa_phase, 
            distance=distance, 
            delta_t=1.0/fs, 
            f_lower=f_lower)
   
    #mb selecting optimal orientation corresponding 
    # to a given gps time 

    right_ascension, declination = det.optimal_orientation(gps)
    polarization = np.random.uniform(0, 2*np.pi)

    signal   = det.project_wave(hp, hc, right_ascension, declination, polarization)
    metadata = [mass1, mass2, spin1z, spin2z, inclination, coa_phase, distance, right_ascension, declination, polarization]

    return signal.numpy(), metadata 


def random_waveform_parameters(detector, gps, m1, m2, s1, s2, di): 
    
    """
    Select random parameters for the waveform generation:
    mass1, mass2, spin1, spin2, distance from ranges defined 
    in the config file.

    Other parameters - inclination, coa_phase, polarization - 
    selected randomly from an uniform distribution 

    Sky position selected based on the gps time, optimally 
    for a given detector det 
    """

    det = Detector(detector)

    # random choice of parameters for masses, spins, 
    # inclination, phase at coalescence, distance    
    mass1 = np.random.uniform(m1[0], m1[1]) 
    mass2 = np.random.uniform(m2[0], m2[1])

    spin1z = np.random.uniform(s1[0], s1[1]) 
    spin2z = np.random.uniform(s2[0], s2[1])

    inclination = np.random.uniform(0, np.pi) 
    coa_phase = np.random.uniform(0, 2*np.pi)

    distance = np.random.uniform(di[0], di[1])

    right_ascension, declination = det.optimal_orientation(gps)

    polarization = np.random.uniform(0, 2*np.pi)

    metadata = [mass1, mass2, spin1z, spin2z, inclination, coa_phase, distance, right_ascension, declination, polarization]

    return metadata
 

def get_signal_waveform(detector, gps, apx, metadata): 
    
    det = Detector(detector)

    # metadata is [mass1, mass2, spin1z, spin2z, inclination, coa_phase, distance, right_ascension, declination, polarization]

    #mb this generates the waveform from the metadata 
    hp, hc = get_td_waveform(approximant=apx,
            mass1=metadata[0], 
            mass2=metadata[1],
            spin1z=metadata[2], 
            spin2z=metadata[3],
            inclination=metadata[4], 
            coa_phase=metadata[5], 
            distance=metadata[6], 
            delta_t=1.0/fs, 
            f_lower=f_lower)

    signal   = det.project_wave(hp, hc, metadata[7], metadata[8], metadata[9])

    return signal.numpy()


def signal_snr(signal, zoom, psd, gps): 

    # shifting the peak of signal to the middle of time series 
    sigshift = duration/2.0 - float(len(signal))/fs;           

    ##mb !!! random shift +- 0.1 -- 0.4 s 
    randshift = random.uniform(0.1, 0.4)  
    sigshift = sigshift + randshift 
    
    #sigshift has to be an exact multiple to zoom.dt otherwise the injection will not work after
    factor = int(sigshift/zoom.dt.value)
    sigshift = zoom.dt.value * factor


    sig = TimeSeries(signal, sample_rate=zoom.sample_rate, epoch=(gps + sigshift)).taper()
    sig.resize(len(zoom))

    # Signal-to-noise of waveform sig  
    snr_timed = matched_filter(sig.to_pycbc(), sig.to_pycbc(), 
                         psd=psd.to_pycbc(), low_frequency_cutoff=highpass)

    snr_timed_max = np.max(TimeSeries.from_pycbc(snr_timed).abs())

    return sig, snr_timed_max, sigshift 


# hdf5 dump info function 
def dump_info(name, obj):
    print ("{0} :".format(name))
    try:
        print ("   value: {0}".format(obj.value))
        for key in obj.attrs.keys():
            print ("     {0}:  {1}".format(key, obj.attrs[key]))
    except:
        pass


def my_whiten(data, fftlength=None, overlap=0, method='scipy-welch',
               window='hanning', detrend='constant', asd=None,
               fduration=2, highpass=None, **kwargs):
        """Whiten this `TimeSeries` using inverse spectrum truncation

        Parameters
        ----------
        fftlength : `float`, optional
            FFT integration length (in seconds) for ASD estimation,
            default: choose based on sample rate

        overlap : `float`, optional
            number of seconds of overlap between FFTs, defaults to the
            recommended overlap for the given window (if given), or 0

        method : `str`, optional
            FFT-averaging method, default: ``'scipy-welch'``,
            see *Notes* for more details

        window : `str`, `numpy.ndarray`, optional
            window function to apply to timeseries prior to FFT,
            default: ``'hanning'``
            see :func:`scipy.signal.get_window` for details on acceptable
            formats

        detrend : `str`, optional
            type of detrending to do before FFT (see `~TimeSeries.detrend`
            for more details), default: ``'constant'``

        asd : `~gwpy.frequencyseries.FrequencySeries`, optional
            the amplitude spectral density using which to whiten the data,
            overrides other ASD arguments, default: `None`

        fduration : `float`, optional
            duration (in seconds) of the time-domain FIR whitening filter,
            must be no longer than `fftlength`, default: 2 seconds

        highpass : `float`, optional
            highpass corner frequency (in Hz) of the FIR whitening filter,
            default: `None`

        **kwargs
            other keyword arguments are passed to the `TimeSeries.asd`
            method to estimate the amplitude spectral density
            `FrequencySeries` of this `TimeSeries`

        Returns
        -------
        out : `TimeSeries`
            a whitened version of the input data with zero mean and unit
            variance

        See Also
        --------
        TimeSeries.asd
            for details on the ASD calculation
        scipy.signal.fftconvolve
            for details on the convolution algorithm used here
        gwpy.signal.filter_design.fir_from_transfer
            for FIR filter design through spectrum truncation

        Notes
        -----
        The `window` argument is used in ASD estimation, FIR filter design,
        and in preventing spectral leakage in the output.

        Due to filter settle-in, a segment of length `0.5*fduration` will be
        corrupted at the beginning and end of the output.

        The input is detrended to give the whitened `TimeSeries` zero mean,
        and the output is normalised to give it unit variance.

        For more on inverse spectrum truncation, see arXiv:gr-qc/0509116.
        """
        
        # compute the ASD
        fftlength = fftlength if fftlength else _fft_length_default(data.dt)
        if asd is None:
            asd = data.asd(fftlength, overlap=overlap,
                           method=method, window=window, **kwargs)
        asd = asd.interpolate(1./data.duration.decompose().value)
        
        # design whitening filter, with highpass if requested
        ncorner = int(highpass / asd.df.decompose().value) if highpass else 0
        ntaps = int((fduration * data.sample_rate).decompose().value)
        #print ncorner, ntaps
        tdw = gwpy.signal.filter_design.fir_from_transfer(1/asd.value, ntaps=ntaps,
                                              window=window, ncorner=ncorner)
        
        # condition the input data and apply the whitening filter
        in_ = data.copy().detrend(detrend)
        pad = int(np.ceil(tdw.size/2))
        window = scipy.signal.get_window(window, tdw.size)
        in_[:pad] *= window[:pad]
        in_[-pad:] *= window[-pad:]
        out = type(data)(scipy.signal.fftconvolve(in_.value, tdw, mode='same'))
        out.__array_finalize__(data)
        
        return math.sqrt(2 * in_.dt.value) * out ## No normalization!


def read_segment_info():

    seggps_sta, seggps_end, seggps_dur, = np.loadtxt(segmenfile, usecols=(0,1,2), unpack=True)
    
    # sort the list w.r.t the first column  
    seggps_sta, seggps_end, seggps_dur = (np.asarray(t) for t in zip(*sorted(zip(seggps_sta, seggps_end, seggps_dur))))    
    return seggps_sta, seggps_end, seggps_dur


def read_injections():
    injs = []
    inje = []
    if (injfile):
        f=open(injfile,"r")
        lines=f.readlines()
        for x in lines:
            injs.append(float(x.split('\t')[0]))
            inje.append(float(x.split('\t')[1]))
    
    return injs,inje

def read_segment_and_injection():
        
        #read the files with segments of available data and build segments out of it
        data = SegmentList()
        f=open(segmenfile,"r")
        lines=f.readlines()
        for x in lines:
                segment = Segment(float(x.split(' ')[0]),float(x.split(' ')[1]))
                data.append(segment)
        f.close()
        
        #Sort the elements of a list into ascending order, and merge continuous segments into single segments.
        data.coalesce() 

        #read the files with injection (if it exists) and build segments out of it
        injections = SegmentList()
        if (injfile):
                f=open(injfile,"r")
                lines=f.readlines()
                for x in lines:
                        segment = Segment(float(x.split('\t')[0]),float(x.split('\t')[1]))
                        injections.append(segment)
                        f.close()
        
                #Sort the elements of a list into ascending order, and merge continuous segments into single segments.        
                injections.coalesce() 
                
        #remove from the data segments, the segments with injections
        new_segments = data - injections
        #Sort the elements of a list into ascending order, and merge continuous segments into single segments.
        new_segments.coalesce()
        
        return new_segments 
                            

def read_glitch_files():

    glgps, gldur, glsnr = np.loadtxt(glitchfile, usecols=(0,1,2), unpack=True)

    if (np.mean(gldur)>5):
            print("something is wrong in the duration of the glitches, mean value",np.mean(gldur))
            print("the list of glitches should contain 3 columns with")
            print("GPS time     duration        SRN(of the detector) ")
            exit()
              
    #if another file is present concatenate the two ndarray
    if(glitchfile_gspy): 
            glgps_gspy, gldur_gspy, glsnr_gspy = np.loadtxt(glitchfile_gspy, usecols=(0,1,2), unpack=True)
            glgps = np.concatenate((glgps,glgps_gspy), axis=None)
            glsnr = np.concatenate((glsnr,glsnr_gspy), axis=None)   
            gldur = np.concatenate((gldur,gldur_gspy), axis=None)    

    
    return glgps, gldur, glsnr

def read_glitch_files_gspy(sel=None):
    #glgps, gldur, glsnr = np.loadtxt(glitchfile_gspy, usecols=(0,1,2), unpack=True)
    glgps, gldur, glsnr, label = np.loadtxt(glitchfile_gspy, usecols=(0,1,2,3), unpack=True, dtype={'names': ('a','b','c','d'),'formats': (np.float64, np.float64,np.float64,'|S20')})
    
    print("dentro read_glitch_files_gspy") 
    if(sel):
        print("dovrebbe entrare qui",sel)
        #select only one type of glicth
        #index_sel = [index for index, value in enumerate(label) if value == "Blip"]
        index_sel = [index for index, value in enumerate(label) if value == sel]
        glgps_sel = [glgps[i] for i in index_sel]
        gldur_sel = [gldur[i] for i in index_sel]
        glsnr_sel = [glsnr[i] for i in index_sel]
        return glgps_sel, gldur_sel, glsnr_sel
    
    else:            
        return glgps, gldur, glsnr
            
            
    
def get_hoft_glitch_sample(hoft, gpsstart, t, duration): 

    # gps start time of the sample, corrected for the crop 
    gpsstart += crop 

    # index corresponding to the gps time t
    index = int((t - gpsstart)*fs) 

    # slice the [duration] s part of the array 
    d = int(duration*fs/2.)
    return hoft[index - d:index + d]



def available_glitches_excluding_injections():

    """ 
    Checks for glitches lying in between of at least [duration s] 
    to adjascent glitches, start or end of segment taking into account
    also the presence of injections (by excluding segments with injections)
        
    In this version it is assumed that the txt files containing the glitches
    have at leat 3 columns containing: GPS duration SNR        
    """

#    glgps_all, glsnr = np.loadtxt(glitchfile, usecols=(0,2), unpack=True)

#    print('glgps_all',len(glgps_all))
#    print('glsnr',len(glsnr))    
    #if another file is present concatenate the two ndarray
#    if(glitchfile_gspy): 
#            glgps_gspy, glsnr_gspy = np.loadtxt(glitchfile_gspy, usecols=(0,2), unpack=True)
#            glgps_all = np.concatenate((glgps_all,glgps_gspy), axis=None)
#            glsnr = np.concatenate((glsnr,glsnr_gspy), axis=None)   
    

#    glgps_all, gldur, glsnr = read_glitch_files()
    glgps_all, gldur, glsnr = read_glitch_files_gspy(glitch_cat)

    #find elements for which SNR>snr_th
    index_srn = [index for index, value in enumerate(glsnr) if value > snr_th]
    #find elements for which dur<1
    index_dur =  [index for index, value in enumerate(gldur) if value <1]
    #intersection
    index_inter = list(set(index_srn) & set(index_dur))
    glgps = [glgps_all[i] for i in index_inter]
    selsnr = [glsnr[i] for i in index_inter] 

    glgps_snr = zip(glgps,selsnr)
    #list of segments of data excluding injections
    segments_list = read_segment_and_injection()

    # sorting the glitch gps times list
    #gs = np.sort(glgps)
    glgps_snr_sorted = sorted(glgps_snr, key = lambda x: x[0])
    gs, gsnr = (zip(*glgps_snr_sorted))
    gs = np.asarray(gs)
    gsnr = np.asarray(gsnr)

    result = [] 
    for s in segments_list:
   
        # sorted gps times in a given segment  
        a = gs[(gs>=s[0]) & (gs<=s[1])]
        snr = gsnr[(gs>=s[0]) & (gs<=s[1])]
    
        # adding gps start and end to array  
        a = np.insert(a, 0, s[0]) 
        a = np.append(a, s[1]) 

        gsel = [[g, snr[ind], a[ind], a[ind+2]] for ind, g in enumerate(a[1:-1]) if a[ind+2]-a[ind] >= duration]

        # gps time of glitch, times of preceeding and following glitches 
        result.extend(gsel)   

    print(result[:14])###------------cancella
    return result

def choose_glitch_samples(glitches, howmany): 

    # from the list of glitches, select randomly glitch gps time 
    l = [item[0] for item in glitches]
    snr = [item[1] for item in glitches]

    #extract randomly the index
    index = np.linspace(0,len(l)-1,len(l),dtype=int)
    indsel = np.random.choice(index, howmany, replace=False)

    g = [l[i] for i in indsel]
    snrs = [snr[i] for i in indsel]

    g = np.asarray(g)
    snrs = np.asarray(snrs) 

    
    seggps_sta, seggps_end, seggps_dur = read_segment_info()

    nind = np.searchsorted(seggps_sta, g) - 1

    segment_start = [] 
    filename = [] 

    segment_start = [seggps_sta.item(ind) for ind in nind] 
    filename = [detector + '_' + str(int(seggps_sta.item(ind))) + '_' + str(int(seggps_dur.item(ind))) + '.hdf5' for ind in nind] 

    fn, gp, ss, snrp = zip(*sorted(zip(filename, g, segment_start,snrs)))

    return gp, ss, fn, snrp


def get_hoft_noise_sample(hoft, gpsstart, t, duration): 

    # gps start time of the sample, corrected for the crop 
    gpsstart += crop 

    # index corresponding to the gps time t
    index = int((t - gpsstart)*fs) 
   
    # slice the [duration] s part of the array
    return hoft[index:index + int(duration*fs)] 


def available_noise_stretches_glitch_duration_excluding_injections():

    """
    Finds streches of noise of at least [duration s] 
    between glitches (taking into account their duration), 
    start or end of the segment
    
    The segments with known injections are excluded 
        
    In this version it is assumed that the txt files containing the glitches
    have at leat 3 columns containing: GPS duration SNR        
        
    """

    # glitches gps times and glitches duration from glitchfile   
    glgps, gldur, glsnr = read_glitch_files()
    
    # Half-duration 
    gldur_half = [x/2. for x in gldur]

    # sorting the glitch gps times list (and respective half duration)
    # gs, gd = zip(*sorted(zip(glgps, gldur_half)))

    # beginning and end gps times of glitches
    gbeg = [x - y for x, y in zip(glgps, gldur_half)]
    gend = [x + y for x, y in zip(glgps, gldur_half)]
    
    
    # I add the injections to the list of glitches since these are times that have to be excluded
    injs,inje = read_injections()
    
    #I extend the lists of glicthes start and stop with the injections
    gbeg.extend(injs)
    gend.extend(inje)

    be = sorted(zip(gbeg, gend))

    #list of segments of data
    seggps_sta, seggps_end, seggps_dur = read_segment_info() 

    result = []
    glitches = []
    begends = []

    for s0,s1 in sorted(zip(seggps_sta, seggps_end)):

        # beginning and ends of glitches within a segment     
        be_seg = [(x, y) for x, y in be if x > s0 and y < s1]

        """
        Some glitches may be placed on top of each other, 
        for example:  

        1. b1----b2--e2-----e1 
            or 
        2. b1-----------b2-----e1---e2

        What is needed is instead  

        1. b1---------------e1
        2. b1-----------------------e2      
        """

        # iterating from the end of the tuple 
        # to be able to remove elements 
        for i in xrange(len(be_seg) - 1, -1, -1):

            # case 1. (glitch i is within the duration of glitch i-1) 
            # b_{i} > b_{i-1} and e_{i} < e_{i-1}
            if be_seg[i][0] > be_seg[i-1][0] and be_seg[i][1] < be_seg[i-1][1]:
                del be_seg[i]

            # case 2. (glitch i is partially within the duration of glitch i-1) 
            # b_{i} < e_{i-1} and e_{i} > e_{i-1} 
            elif be_seg[i][0] < be_seg[i-1][1] and be_seg[i][1] > be_seg[i-1][1]:
                # extend the end of the glitch period 
                be_seg[i-1] = (be_seg[i-1][0], be_seg[i][1])
                del be_seg[i]

        gs_seg = [x for x in glgps if x > s0 and x < s1]


        """
        Starting from the full segment length (s0, s1) 
        divide it by removing parts corresponding to glitches: 

        1. s0--------------------------------s1
        2. glitch b1--e1 
        3. s0-----b1  e1---------------------s1 
        4. glitch           b2---e2 
        5. s0-----b1  e1----b2   e2----------s1 
        6. etc. 
        """

        # full length of the segment 
        r = np.array([(s0, s1)])


        # iterative scheme, working on last tuple 
        # of r (glitch beginings and ends are sorted) 
        for b, e in be_seg:

            r = np.append(r, [(r[-1][0], b)], axis=0)
            r = np.append(r, [(e, r[-2][1])], axis=0)
            # remove the undivided part 
            r = np.delete(r, -3, axis=0)

        # list of gps time start of noise sample, its duration, 
        # gps time of segment's start and its duration 
        gaps = [[x, y - x, s0, s1-s0] for x, y in r if y - x > duration]


        result.extend(gaps)

    return result

def choose_noise_samples(noise_stretches, howmany): 

    """
    Having a set of selected stretches of at least `duration` s - noise_stretches -  
    select `howmany` random gps start times that allow for `duration` length, 
    with corresponding segment data (gps time segment_start and filename) 
    """

    # maximal duration of each noise stretch to accomodate the preset duration 
    gsegl = [] 
    for d in noise_stretches: 
        gsegl.append(d[1] - duration) 

    # cummulative sum of the lengths (adding 0 at the begining) 
    cumsum = np.insert(np.cumsum(gsegl), 0, 0) 
    
    # selecting random numbers from the range 0, cumsum[-1] 
    choices = np.random.uniform(0, cumsum[-1], howmany) 

    seggps_sta, seggps_end, seggps_dur = read_segment_info()

    gps = [] 
    segment_start = []  
    filename = [] 

    # checking in which time-series sample the randomly chosen 
    # numbers belong; returning gps times, segment start times, segment file names     
    for c in choices:
        
        ind = np.searchsorted(cumsum, c)

        g = noise_stretches[ind-1][0] + c - cumsum[ind-1] 
        gps.append(g) 

        nind = np.searchsorted(seggps_sta, g) - 1
        segment_start.append(seggps_sta[nind]) 

        filename.append(detector + '_' + str(int(seggps_sta[nind])) + '_' + str(int(seggps_dur[nind])) + '.hdf5') 

    fn, gp, ss = zip(*sorted(zip(filename, gps, segment_start)))

    return list(gp), list(ss), list(fn) #list is needed to convert tuples in lists
 
#    return gps, segment_start, filename 
   

#def get_random_batch(sigfile, batchsize, selection=None):                   
#
#    with h5py.File(sigfile, "r") as f:
#
#        groups = f.keys()
#
#        if not selection:
#            import random  
#            selection = random.sample(groups, batchsize)   
#
#        signal = [] 
#        metadata = [] 
#
#        for s in selection:
#        
#            grp = f[s]
#        
#            signal.append(grp['signal'][:]) 
#            metadata.append(grp['metadata'][:])
# 
#    f.close() 
# 
#    # subset of groups with indices selection
#    return signal, metadata, selection


def output_noise_selection(howmany, add_signal=False): 

#    noise_stretches = available_noise_stretches()
    
    noise_stretches = available_noise_stretches_glitch_duration_excluding_injections()

    print("Found {} noise stretches of at least {} s length in available segment '{}' and glitch '{}' lists".format(len(noise_stretches), duration, segmenfile, glitchfile))

    print("Selecting {} noise samples...".format(howmany)) 

    gps, segment_start, filename = choose_noise_samples(noise_stretches, howmany)

    # chose the file according if it's the signal or the noise case
    if(add_signal):
            outfile = sigoutfile
    else:
            outfile = noiseoutfile

    # write samples to file             
    with h5py.File(outfile, "w") as outf:

        prevfilename = None 
        i=0 #index of the list of selected gps times
        samples_selected=0 #counter of the number of selected samples
        while (samples_selected<howmany):
        #for i in range(howmany):
            try:

                    print("Selecting noise sample #{} at gps {} from segment {} ...".format(samples_selected, gps[i], filename[i])) 
    
                    # Assuming sorted filename from choose_noise_samples() 
                    if prevfilename == filename[i]:
                            pass 

                    else:  
                            # load new data 
                            f = h5py.File(segmentdir + filename[i], 'r')
                            hoft = f['Strain']['Whitened_Strain'][()]
                            
                            df   = f['PSD']['PSD'].attrs['Delta_f']
                            psd  = FrequencySeries(f['PSD']['PSD'][()], df=df)
                            asd  = FrequencySeries(np.sqrt(f['PSD']['PSD'][()]), df=df)
                            f.close()
 
                            prevfilename = filename[i] 

                    h = get_hoft_noise_sample(hoft, segment_start[i], gps[i], duration) 

                    # normality test (D'Agostino-Pearson) #-------------------------------------
                    #Test whether a sample differs from a normal distribution.
                    statistic, pvalue = sc.normaltest(h)
                    #the second argument is the 2-sided chi squared probability for the hypothesis test
                    #print("pvalue:",pvalue,i)
                    if (pvalue<pcrit or h.size<int(duration*fs)):
                            
                            raise Exception("pvalue<pcrit")
                            
                    #case in which a signal nedds to be added        
                    if (add_signal):
                            write_signal(h, fs, gps[i],segment_start[i],psd,asd,outf)
                    
                    #only noise case        
                    if not add_signal:
                            # Create a hdf5 group for each noise instance, 
                            # containing the array and its possible meta parameters 
                            grp = outf.create_group(str(gps[i]))
   
                            ts  = grp.create_dataset("time series", data=h)
                            
                            psd_file  = grp.create_dataset("psd", data=psd)
                            # meta data is: gps segment_start, gps at the beginning of ts, duration [s] 
                            met = grp.create_dataset("metadata", data=[segment_start[i], gps[i], duration, fs, df])


                    samples_selected = samples_selected+1            
                    i = i+1 #update the index of the list

            except:
                    print("\nD'Agostino-Pearson not passed (or size less than expected), extracting a new noise sample")
                    #extract a new random gps time and add it to the existing list
                    gps_new, segment_start_new, filename_new = choose_noise_samples(noise_stretches, 1)
                    #use method extend because the output of the previous method consist in lists
                    gps.extend(gps_new)
                    segment_start.extend(segment_start_new)
                    filename.extend(filename_new)
                    i=i+1 #update the index of the list 
      
    outf.close() 


def output_glitch_selection(howmany): 

    glitches = available_glitches_excluding_injections()

    print("Found {} glitches in gaps of at least {} s length in available segment '{}' and glitch '{}' lists".format(len(glitches), duration, segmenfile, glitchfile))

    print("Selecting {} glitches:".format(howmany)) 

    gps, segment_start, filename, snr = choose_glitch_samples(glitches, howmany)

    # write samples to file 
    with h5py.File(glitchoutfile, "w") as outf:

        prevfilename = None

        for i in range(howmany):

            print("Selecting glitch #{} at gps {} from segment {} ...".format(i, gps[i], filename[i])) 

            # Assuming sorted filename from choose_noise_samples() 
            if prevfilename == filename[i]:
                pass

            else:
                # load new data
                f = h5py.File(segmentdir + filename[i], 'r')
                hoft = f['Strain']['Whitened_Strain'][()]
                f.close() 

                prevfilename = filename[i]

            h = get_hoft_glitch_sample(hoft, segment_start[i], gps[i], duration) 
            

#            # adding some additional glitchiness 
#            fc = 50*np.random.random_sample() + 50 
#            amp = 2 + np.random.random_sample()
                   
#            t = np.linspace(-1, 1, duration*fs, endpoint=False)
#            r = amp*scipy.signal.gausspulse(t, fc=fc, retquad=False, retenv=False)

#            h = h + r  

            # Create a hdf5 group for each noise instance, 
            # containing the array and its possible meta parameters 
            grp = outf.create_group(str(gps[i]))
   
            ts  = grp.create_dataset("time series", data=h)
            # metadata is: gps segment_start, gps at the beginning of ts, duration [s] 
            met = grp.create_dataset("metadata", data=[segment_start[i], gps[i], duration,snr[i]])

    outf.close() 


def output_signal_selection(howmany): 

#    noise_stretches = available_noise_stretches()

    noise_stretches = available_noise_stretches_glitch_duration_excluding_injections()

    print("Found {} noise stretches of at least {} s length in available segment '{}' and glitch '{}' lists".format(len(noise_stretches), duration, segmenfile, glitchfile))

    print("Selecting {} noise samples...".format(howmany)) 

    gps, segment_start, filename = choose_noise_samples(noise_stretches, howmany)

    # write samples to file 
    with h5py.File(sigoutfile, "w") as outf:

        prevfilename = None 

        for i in range(howmany):

            print("Selecting noise sample #{} at gps {} from segment {} ...".format(i, gps[i], filename[i])) 

            # Assuming sorted filename from choose_noise_samples() 
            if prevfilename == filename[i]:
                pass

            else:

                # load new data 
                f = h5py.File(segmentdir + filename[i], 'r')

                hoft = f['Strain']['Whitened_Strain'][()]


                df   = f['PSD']['PSD'].attrs['Delta_f']
                psd  = FrequencySeries(f['PSD']['PSD'][()], df=df)
                asd  = FrequencySeries(np.sqrt(f['PSD']['PSD'][()]), df=df)
         
                f.close() 
                prevfilename = filename[i] 
     
            h = get_hoft_noise_sample(hoft, segment_start[i], gps[i], duration) 
            zoom = TimeSeries(h, sample_rate=fs, epoch=gps[i])

            snr_random_selection = np.random.uniform(snr_limit_min, snr_limit_max, 1)

            metadata = random_waveform_parameters(detector, gps[i], m1, m2, s1, s2, di)

            # generating a waveform using the metadata values  
            signal = get_signal_waveform(detector, gps[i], apx, metadata)
            sig, snr_timed_max, _ = signal_snr(signal, zoom, psd, gps[i]) 

            #print("Have: ", snr_timed_max, " selected: ", snr_random_selection, metadata[6], metadata[0], metadata[1])
            # correctin the distance to get randomly selected snr 
            metadata[6] = metadata[6]*(snr_timed_max/snr_random_selection)

            # generating a waveform using the metadata values  
            signal = get_signal_waveform(detector, gps[i], apx, metadata)
            sig, snr_timed_max, sigshift = signal_snr(signal, zoom, psd, gps[i]) 

            #print("Modified: ", snr_timed_max, " selected: ", snr_random_selection, metadata[6], metadata[0], metadata[1])

            sig_whiten = my_whiten(sig, fftlength=fftlength, asd=asd, fduration=fduration, highpass=highpass)

            # Inject (add) whiten signal to the background
            bck_plus_sig = zoom.inject(sig_whiten) 

            # Create a hdf5 group for each instance, 
            # containing the time series, signal 
            # and their meta parameters 
            grp = outf.create_group(str(gps[i]))
  
            ts = grp.create_dataset("time series", data=bck_plus_sig.value)

#-------------------------
            #Timeseries containing only the signal
            #I inject the whitened signal in a timeseries made all by zeros with the epoch=gps[i]
            #in this way the signal will be in the correct position
            zeros_ts = TimeSeries(np.zeros(len(h)), sample_rate=fs, epoch=gps[i])
            just_signal = zeros_ts.inject(sig_whiten)
#-------------------------
            #just_signal = np.insert(sig_whiten.value[:int((duration-sigshift)*fs)], 0, np.zeros(int(sigshift*fs)+1), axis=0)
            
            si = grp.create_dataset("signal", data=just_signal.value)            
            
            # metadata is: gps segment_start, gps at the beginning of ts, duration [s]
            # signal parameters metadata[i] and the SNR
            m = np.append(np.array([segment_start[i], gps[i], duration/4]), metadata, 0) 
            m = np.append(m, np.array([snr_timed_max]), 0) 

            m = np.array(m.flatten(), dtype=np.float)  

            met = grp.create_dataset("metadata", data=m)

    outf.close() 


def testing_plot(filename, plotname): 

    # a few random instances 
    xnum = 2 
    ynum = 2 
    howmany = xnum*ynum 

    import matplotlib
    matplotlib.use('agg')
 
    import matplotlib.pyplot as plt
    matplotlib.use('Agg')
    plt.rc('legend', fontsize=6) 
    plt.rc('xtick', labelsize=8) 
    plt.rc('ytick', labelsize=8) 

    fig, axs = plt.subplots(xnum, ynum)
    fig.tight_layout()

    f = h5py.File(filename, 'r')

    #import random 
    groups = f.keys()

    random_selection_of_keys = random.sample(groups, howmany) 

    import itertools 
    a = list(itertools.product(range(xnum), range(ynum)))

    for i, key in enumerate(random_selection_of_keys):

        grp = f[key]

        h = grp['time series'][:]

        gpsstart = grp['metadata'][1]
        duration = grp['metadata'][2]
        snr      = grp['metadata'][13] 
 
        t = np.arange(gpsstart, gpsstart + duration, 1./fs)

        axs[a[i][0], a[i][1]].plot(t, h, label="{}, SNR: {:.3f}".format(gpsstart, snr))
        axs[a[i][0], a[i][1]].legend(loc="best")

        try: 
            s = grp['signal'][:]
            axs[a[i][0], a[i][1]].plot(t, s, label="{}, SNR: {:.3f}".format(gpsstart, snr))
        except: 
            pass 

        print("Metadata:\n")  
        print(grp['metadata'][:]) 

    print("Test plot at {}".format(testplotwww + plotname))
    testplot = testplotdir + plotname

    plt.savefig(testplot, format="pdf")




def write_signal(h, fs, gps, segment_start, psd,asd,outf):
        """
        method to add the signal to the noise and write it to the output file 
        """
        
        zoom = TimeSeries(h, sample_rate=fs, epoch=gps)

        snr_random_selection = np.random.uniform(snr_limit_min, snr_limit_max, 1)

        metadata = random_waveform_parameters(detector, gps, m1, m2, s1, s2, di)

        # generating a waveform using the metadata values  
        signal = get_signal_waveform(detector, gps, apx, metadata)
        sig, snr_timed_max, _ = signal_snr(signal, zoom, psd, gps) 

        #print("Have: ", snr_timed_max, " selected: ", snr_random_selection, metadata[6], metadata[0], metadata[1])
        # correctin the distance to get randomly selected snr 
        metadata[6] = metadata[6]*(snr_timed_max/snr_random_selection)

        # generating a waveform using the metadata values  
        signal = get_signal_waveform(detector, gps, apx, metadata)
        sig, snr_timed_max, sigshift = signal_snr(signal, zoom, psd, gps) 

        #print("Modified: ", snr_timed_max, " selected: ", snr_random_selection, metadata[6], metadata[0], metadata[1])

        sig_whiten = my_whiten(sig, fftlength=fftlength, asd=asd, fduration=fduration, highpass=highpass)

        # Inject (add) whiten signal to the background
        bck_plus_sig = zoom.inject(sig_whiten) 
        print("---> adding signal with SNR",snr_timed_max)

        # Create a hdf5 group for each instance, 
        # containing the time series, signal 
        # and their meta parameters 
        grp = outf.create_group(str(gps))
  
        ts = grp.create_dataset("time series", data=bck_plus_sig.value)

        #Timeseries containing only the signal
        #I inject the whitened signal in a timeseries made all by zeros with the epoch=gps
        #in this way the signal will be in the correct position
        zeros_ts = TimeSeries(np.zeros(len(h)), sample_rate=fs, epoch=gps)
        just_signal = zeros_ts.inject(sig_whiten)
            
        si = grp.create_dataset("signal", data=just_signal.value)            
                    
        # metadata is: gps segment_start, gps at the beginning of ts, duration [s]
        # signal parameters metadata and the SNR
        m = np.append(np.array([segment_start, gps, duration/4]), metadata, 0) 
        m = np.append(m, np.array([snr_timed_max]), 0) 

        m = np.array(m.flatten(), dtype=np.float)  

        met = grp.create_dataset("metadata", data=m)






def main(): 


    start_time = time.time()

    howmany = int(sys.argv[1]) 

    np.set_printoptions(precision=3, suppress=True, formatter={'float_kind':'{:3f}'.format})

    try:
        output = sys.argv[3]
    except:
        print('!!!default case glitch!!!')
        print('add one of the following options:  noise glitch signal')
        output = 'glitch'

    if (output=='noise'):
        print('case:',output)
        output_noise_selection(howmany)

    elif(output=='glitch'):
        print('case:',output)
        output_glitch_selection(howmany)

    elif(output=='signal'):
        print('case:',output)
        output_noise_selection(howmany,add_signal=True) 

    else:
        print('add one of the following options: noise glitch signal')

    print("--- %s seconds ---" % (time.time() - start_time))
   
#    output_noise_selection(howmany)     #------------------------------------- for noise
#    testing_plot(noiseoutfile, 'noise.pdf') 

#    # load data 
#    f = h5py.File(noiseoutfile, "r")
#    f.visititems(dump_info)    
#    f.close() 

    #output_glitch_selection(howmany) #------------------------------------- for glitches     
#    testing_plot(glitchoutfile, 'loud_glitches.pdf')

#    # load data 
#    f = h5py.File(glitchoutfile, "r")
##    f.visititems(dump_info)    
#    f.close() 

#    output_signal_selection(howmany)     #------------------------------------- for signal
#    output_noise_selection(howmany,add_signal=True)     #------------------------------------- for signal with D'agostino-pearson test

#    testing_plot(sigoutfile, 'signals_test.pdf')

#    # load data 
#    f = h5py.File(sigoutfile, "r")
##    f.visititems(dump_info)    
#    f.close() 
#


if __name__ == "__main__": 
    main() 



