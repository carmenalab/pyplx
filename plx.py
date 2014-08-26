import plxread 
import scipy.io, os 
from numpy import *
from scipy.io import savemat

STROBED_CHANNEL = 257

def convert_plx_file_to_mat(plx_fname, save_fname, AD_channels=range(33,41),
        save_wf=False):
    """ Import a plexon file and save it in matlab format

    Parameters
    ----------
    plx_fname : plexon source file name
    save_fname : destination file name 

    Returns
    -------
    None
    """
    if not os.path.exists(plx_fname):
        print "File does not exist: %s" % plx_fname
        raise Exception

    if not isinstance(save_wf, bool) or isinstance(save_wf, basestring):
        raise Exception('Improper dtype: save_wf')
    store_wf = int(bool(save_wf))
    
    data = plxread.import_file(plx_fname, AD_channels=AD_channels, 
        store_wf=save_wf)

    # Empty dict for matlab/numpy formatted data 
    mdata = {}

    # pull spike timestamps out of data
    unit_lut = {0:'i', 1:'a', 2:'b', 3:'c', 4:'d' }
    spike_keys = filter( lambda x: isinstance(x, tuple), data.keys() )
    spike_ts_keys = filter( lambda x : len(x) == 2, spike_keys)
    spike_wf_keys = filter( lambda x : len(x) == 3, spike_keys)

    for key in spike_ts_keys:
        spike_name = 'sig%03d%s' % (key[0], unit_lut[key[1]])
        mdata[spike_name] = mat( data[key] ).T

    for key in spike_wf_keys:
        spike_name = 'sig%03d%s_%s' % (key[0], unit_lut[key[1]], key[2])
        mdata[spike_name] = mat( data[key] ).T
        

    strobed = data[STROBED_CHANNEL] # plxread.get_strobed(plx_fname)
    strobe_times = mat([ s[0] for s in strobed]).T
    strobe_codes = mat([ s[1] for s in strobed]).T
    strobed = hstack((strobe_times, strobe_codes))
    mdata['Strobed'] = strobed

    #channels = plxread.get_AD_channels(plx_fname, channels=range(33,41) )
    for ch in AD_channels:
        mdata[ "AD%02d" % ch ] = mat(data[ch]).T
    mdata['AD33_ts'] = data['t_start']
    mdata['t_start'] = data['t_start']

    scipy.io.savemat(save_fname, mdata)


def get_strobed(plx_fname, save_fname):
    mdata = {}
    strobed = plxread.get_strobed(plx_fname)
    strobe_times = mat([ s[0] for s in strobed]).T
    strobe_codes = mat([ s[1] for s in strobed]).T
    strobed = hstack((strobe_times, strobe_codes))
    mdata['Strobed'] = strobed

    savemat(save_fname, mdata)
