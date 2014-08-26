/////////////////////////////////////////////////////////////////////
// plxread.cpp - sample functionality for reading PLX files.
//
// Derivative work, based on 


#include "plxread.h"

void 
load_header_dict(PyObject* header_dict)
{
    //
    // add metadata
    //
    PyObject* metadata_key = PyString_FromString("metadata");
    PyObject* metadata_dict = PyDict_New();

    PyObject* date = PyTuple_Pack(6, PyInt_FromLong(fileHeader.Year), PyInt_FromLong(fileHeader.Month), PyInt_FromLong(fileHeader.Day), PyInt_FromLong(fileHeader.Hour), PyInt_FromLong(fileHeader.Minute), PyInt_FromLong(fileHeader.Second ));
    PyDict_SetItem( metadata_dict, PyString_FromString("Version"), PyInt_FromLong(fileHeader.Version));
    PyDict_SetItem( metadata_dict, PyString_FromString("Date"), date );
    PyDict_SetItem( metadata_dict, PyString_FromString("LastTimestamp"), PyFloat_FromDouble(fileHeader.LastTimestamp));
    PyDict_SetItem( metadata_dict, PyString_FromString("NumPointsPreThr"), PyInt_FromLong(fileHeader.NumPointsPreThr));
    PyDict_SetItem( metadata_dict, PyString_FromString("NumPointsWave"), PyInt_FromLong(fileHeader.NumPointsWave));
    PyDict_SetItem( metadata_dict, PyString_FromString("WaveformFreq"), PyInt_FromLong(fileHeader.WaveformFreq));
    PyDict_SetItem( metadata_dict, PyString_FromString("ADFrequency"), PyInt_FromLong(fileHeader.ADFrequency));

    if (fileHeader.Version >= 103)
    {
        PyDict_SetItem( metadata_dict, PyString_FromString("SlowMaxMagnitudeMV"), PyInt_FromLong(fileHeader.SlowMaxMagnitudeMV));
        PyDict_SetItem( metadata_dict, PyString_FromString("SpikeMaxMagnitudeMV"), PyInt_FromLong(fileHeader.SpikeMaxMagnitudeMV));
        PyDict_SetItem( metadata_dict, PyString_FromString("NumSlowChannels"), PyInt_FromLong(fileHeader.NumSlowChannels));
        PyDict_SetItem( metadata_dict, PyString_FromString("NumEventChannels"), PyInt_FromLong(fileHeader.NumEventChannels));
        PyDict_SetItem( metadata_dict, PyString_FromString("NumDSPChannels"), PyInt_FromLong(fileHeader.NumDSPChannels));
        PyDict_SetItem( metadata_dict, PyString_FromString("BitsPerSpikeSample"), PyInt_FromLong(fileHeader.BitsPerSpikeSample));
        PyDict_SetItem( metadata_dict, PyString_FromString("DataTrodalness"), PyInt_FromLong(fileHeader.DataTrodalness));
        PyDict_SetItem( metadata_dict, PyString_FromString("Trodalness"), PyInt_FromLong(fileHeader.Trodalness));
    }
    PyDict_SetItem( header_dict, metadata_key, metadata_dict );
    
    //
    // add spike_counts data
    //
    PyObject* spike_counts_key = PyString_FromString("spike_counts");
    PyObject* spike_counts_dict = PyDict_New();
    // Dump the spike channel counts
    int iChannel, iUnit;
    for (iChannel = 0 ; iChannel < 130 ; iChannel++)
    {
        for (iUnit = 0 ; iUnit < 5 ; iUnit++)
        {
            int count = fileHeader.TSCounts[iChannel][iUnit] ;
            if (count > 0)
            {
                int countWF = fileHeader.WFCounts [iChannel][iUnit] ;
                printf("  Channel %d, Unit %d: Time Stamp %d, Waveform %d\n", iChannel, iUnit, count, countWF) ;
                PyObject* sig = PyTuple_Pack(2, PyInt_FromLong(iChannel), PyInt_FromLong(iUnit));
                PyObject* sig_info = PyTuple_Pack(2, PyInt_FromLong(count), PyInt_FromLong(countWF));
                PyDict_SetItem( spike_counts_dict, sig, sig_info);
            }
        }
    }
    PyDict_SetItem( header_dict, spike_counts_key, spike_counts_dict );
    

    //
    // Event Counts
    //
    PyObject* event_counts_key = PyString_FromString("event_counts");
    PyObject* event_counts_dict = PyDict_New();

    printf ("\nEvent Counts\n") ;
    for (iChannel = 0 ; iChannel < 300 ; iChannel++)
    {
        int count = fileHeader.EVCounts[iChannel] ;
        if (count > 0)
        {
            printf ("  Index %d Event Channel %d: Count %d\n", iChannel, iChannel, count) ;
            PyDict_SetItem( event_counts_dict, PyInt_FromLong(iChannel), PyInt_FromLong(count));
        }
    }
    PyDict_SetItem( header_dict, event_counts_key, event_counts_dict );


    //
    // Slow Channel Counts
    //

    PyObject* slow_ch_counts_key = PyString_FromString("slow_ch_counts");
    PyObject* slow_ch_counts_dict = PyDict_New();
    for (iChannel = 300 ; iChannel < 512 ; iChannel++)
    {
        int count = fileHeader.EVCounts[iChannel] ;
        if (count > 0)
        {
            printf ("  Index %d Slow Channel %d: Count %d\n", iChannel, iChannel-300+1, count) ;
            PyDict_SetItem( slow_ch_counts_dict, PyInt_FromLong(iChannel), PyInt_FromLong(count));
        }
    }
    PyDict_SetItem( header_dict, slow_ch_counts_key, slow_ch_counts_dict );

    return;
}

int
read_plx_headers(char* fname)
{
    // Open the specified PLX file.
    fp = fopen(fname, "rb");
    if(fp == 0){
        cout << "Cannot open PLX file: " << fname << endl;
        return -1;
    }

    // Read the file header
    fread(&fileHeader, sizeof(fileHeader), 1, fp);

    if (fileHeader.Version >= 103)
        max_spike_magnitude = fileHeader.SpikeMaxMagnitudeMV; 
    if (fileHeader.Version >= 105) {
        spike_preamp_gain = fileHeader.SpikePreAmpGain; 
        signed_adc_levels = 0.5*pow(2.0, fileHeader.BitsPerSpikeSample);
    }

    // Read the spike channel headers
    if(fileHeader.NumDSPChannels > 0)
        fread(spikeChannels, fileHeader.NumDSPChannels*sizeof(PL_ChanHeader), 1, fp);

    // Read the event channel headers
    if(fileHeader.NumEventChannels> 0)
        fread(eventChannels, fileHeader.NumEventChannels*sizeof(PL_EventHeader), 1, fp);

    // Read the slow A/D channel headers
    if(fileHeader.NumSlowChannels)
        fread(slowChannels, fileHeader.NumSlowChannels*sizeof(PL_SlowChannelHeader), 1, fp);

    // save the position in the PLX file where data block begin
    data_start = sizeof(fileHeader) + fileHeader.NumDSPChannels*sizeof(PL_ChanHeader)
                        + fileHeader.NumEventChannels*sizeof(PL_EventHeader)
                        + fileHeader.NumSlowChannels*sizeof(PL_SlowChannelHeader);
    return data_start;
}

PyObject *
py_read_plx_headers(PyObject *self, PyObject *args)
{
    char* fname;
    if (!PyArg_ParseTuple(args, "s", &fname))
        return Py_None;

    // get the header data into the struct
    data_start = read_plx_headers( fname );

    // create new python dictionary
    PyObject* header_dict = PyDict_New();
    load_header_dict(header_dict);

    return header_dict;
}

PyObject* import_file(PyObject* self, PyObject* args, PyObject* kwargs)
{
    bool verbose = 0;
    // Get data filename
    char* fname; 
    PyArg_ParseTuple( args, "s", &fname) ;

    int data_start = read_plx_headers(fname);  
    if (data_start == -1) 
        return Py_None;

    bool AD_channels_to_import[fileHeader.NumSlowChannels+1];
    for (int i = 0; i < fileHeader.NumSlowChannels+1; i += 1) {
        AD_channels_to_import[i] = 0;
    }

    bool store_wf = 0; 
    
    PyObject* AD_channels;
    // fill kwarg structs
    if (kwargs != NULL) {
        // if ( PyDict_Contains(kwargs, PyString_FromString("sig_names")) ) 
        //     sig_names = PyDict_GetItem(kwargs, PyString_FromString("sig_names"));
        if ( PyDict_Contains(kwargs, PyString_FromString("AD_channels")) ) 
            AD_channels = PyDict_GetItem(kwargs, PyString_FromString("AD_channels"));
            for (int i=0; i < PyList_Size(AD_channels); i += 1) {
                AD_channels_to_import[PyInt_AsLong(PyList_GetItem(AD_channels, i))] = 1;
            }
        if ( PyDict_Contains(kwargs, PyString_FromString("store_wf")) ) 
            store_wf = (bool) PyInt_AS_LONG(PyDict_GetItem(kwargs, PyString_FromString("store_wf")));
        if ( PyDict_Contains(kwargs, PyString_FromString("verbose")) ) 
            verbose = (bool) PyInt_AS_LONG(PyDict_GetItem(kwargs, PyString_FromString("verbose")));
    }

    PyObject* data = PyDict_New();

    int n_samples;
    char* ch = (char*) malloc(sizeof(char)*6);

    float** AD_inds = (float**) malloc(sizeof(float*) * fileHeader.NumSlowChannels);
    for (int i=0; i < fileHeader.NumSlowChannels; i+=1 ){
        AD_inds[i] = 0;
    }

    int k;
    for (k=0; k < fileHeader.NumSlowChannels; k+=1 )
    {
        n_samples = fileHeader.EVCounts[k+300];
        if (n_samples > 0 && AD_channels_to_import[k+1]) 
        {
            int arr_size[] = {n_samples};
            PyObject* ls = PyArray_FromDims(1, arr_size, NPY_FLOAT32); 
            AD_inds[k] = (float*) ((PyArrayObject*) ls)->data;
            snprintf(ch, 6, "AD%d", k+1); 
            PyDict_SetItemString( data, ch, ls );
            Py_DECREF(ls);
        }
    }

    int Channel, Unit, n_spikes;
    char unit_lut[] = "iabcd\0";
    char* sig_name = (char*) malloc(sizeof(char)*8);
    double** spike_ts_inds = (double**) malloc(sizeof(double*)*(MAX_SPIKE_CHANNELS+1)*MAX_SPIKES_PER_ELECTRODE);
    for (int i = 0; i < MAX_SPIKE_CHANNELS*MAX_SPIKES_PER_ELECTRODE; i+=1) {
        spike_ts_inds[i] = 0;
    }

    if (verbose) {
        cout << "creating dict of sig names" << endl;
    }

    for (Channel = 0 ; Channel < 130 ; Channel++) 
    {
        for (Unit = 1 ; Unit < MAX_SPIKES_PER_ELECTRODE ; Unit++) 
        {
            n_spikes = fileHeader.TSCounts[Channel][Unit] ;
            if (n_spikes > 0) 
            {
                snprintf(sig_name, 8, "sig%03d%c", Channel, (char) unit_lut[Unit]); 
                int arr_size[] = {n_spikes}; 
                if (verbose)
                    cout << sig_name << endl;
                PyObject* ls = PyArray_FromDims(1, arr_size, NPY_DOUBLE); 
                spike_ts_inds[MAX_SPIKES_PER_ELECTRODE*Channel+Unit] = 
                    (double*) ((PyArrayObject*) ls)->data;
                PyDict_SetItemString(data, sig_name, ls);
                Py_DECREF(ls);
            }
        }
    }

    int wf_arr_size[] = {0, 32};
    char* wf_name = (char*) malloc(sizeof(char)*12);
    short** spike_wf_inds = (short**) malloc(sizeof(short*)*(MAX_SPIKE_CHANNELS+1)*MAX_SPIKES_PER_ELECTRODE); 
    for (int i = 0; i < MAX_SPIKE_CHANNELS*MAX_SPIKES_PER_ELECTRODE; i+=1) {
        spike_wf_inds[i] = 0;
    }
    if (store_wf) {
        if (verbose)
            cout << "Adding _wf elements to dict" << endl;
        for (Channel = 0 ; Channel < 130 ; Channel++) 
        {
            for (Unit = 1 ; Unit < MAX_SPIKES_PER_ELECTRODE ; Unit++) 
            {
                n_spikes = fileHeader.TSCounts[Channel][Unit];
                //cout << n_spikes << endl;
                if (n_spikes > 0)  
                {
                    snprintf(wf_name, 11, "sig%03d%c_wf", Channel, unit_lut[Unit]); 
                    //cout << wf_name << endl;
                    wf_arr_size[0] = n_spikes;
                    PyObject* ls = PyArray_FromDims(2, wf_arr_size, NPY_INT16); 
                    spike_wf_inds[MAX_SPIKES_PER_ELECTRODE*Channel+Unit] = 
                        (short*) ((PyArrayObject*) ls)->data;
                    PyDict_SetItemString(data, wf_name, ls);
                    Py_DECREF(ls);
                }
            }
        }
    }

    if (verbose)
        cout << "Allocating array for wf_gain" << endl;
    int wf_gain_arr_size[] = {130};
    PyObject* wf_gain_arr = PyArray_FromDims(1, wf_gain_arr_size, NPY_DOUBLE); 
    PyDict_SetItemString(data, "wf_gain", wf_gain_arr);
    Py_DECREF(wf_gain_arr);
    double* wf_gain = (double*) ((PyArrayObject*) wf_gain_arr)->data;
    
    // log events besides Strobed? 
    double** evt_data_inds = (double**) malloc(sizeof(double*)*MAX_EVENT_CHANNELS);
    int strobed_size[] = { fileHeader.EVCounts[STROBED_CHANNEL], 2 };
    PyObject* strobed_arr = PyArray_FromDims(2, strobed_size, NPY_DOUBLE);
    PyDict_SetItemString(data, "Strobed", strobed_arr); 
    evt_data_inds[STROBED_CHANNEL] = (double*) ((PyArrayObject*) strobed_arr)->data;
    Py_DECREF(strobed_arr);

    //
    // Begin file parsing
    // 
    if (verbose)
        cout << "Beginning file parsing" << endl;
    PL_DataBlockHeader dataBlock;
    PL_ChanHeader* pSpikeChannelHeader;
    PL_SlowChannelHeader* pSlowChannel;
    short buf[MAX_SAMPLES_PER_WAVEFORM];

    // Seek to the beginning of the data blocks in the PLX file
    fseek(fp, data_start, SEEK_SET);

    int nbuf ;
    double t_start = -1;

    // Rip through the rest of the file
    for (int iBlock = 0 ; ; iBlock++) {
        // Read the next data block header.
        if (fread(&dataBlock, sizeof(dataBlock), 1, fp) != 1) break ;
        
        if (iBlock % 1000000 == 0 && verbose) {
            cout << iBlock << endl;
            //verbose = 0;
        }

        // Read the waveform samples if present.
        nbuf = 0;
        if(dataBlock.NumberOfWaveforms > 0) {
            nbuf = dataBlock.NumberOfWaveforms*dataBlock.NumberOfWordsInWaveform;
            if (fread(buf, nbuf*2, 1, fp) != 1) break ;
        }

        // Convert the timestamp to seconds
        LONGLONG ts = ((static_cast<LONGLONG>(dataBlock.UpperByteOf5ByteTimestamp)<<32) + static_cast<LONGLONG>(dataBlock.TimeStamp)) ;
        double seconds = (double) ts / (double) fileHeader.ADFrequency ;

        
        int unit_idx = 0;
        if (dataBlock.Type == PL_SingleWFType) { 
            n_spikes = fileHeader.TSCounts[dataBlock.Channel][dataBlock.Unit] ;
            if (n_spikes > 0 && dataBlock.Unit >= 1) 
            {
                unit_idx = MAX_SPIKES_PER_ELECTRODE*dataBlock.Channel + dataBlock.Unit;
                if (verbose) {
                    if (1) {
                    //if (dataBlock.Channel == 128 && dataBlock.Unit == 2) {
                        printf("SingleWFType: Channel: %d, Unit: %d\n", dataBlock.Channel, dataBlock.Unit); 
                        cout << unit_idx << endl;
                        cout << spike_ts_inds[unit_idx] << endl;
                    }
                }

                *(spike_ts_inds[unit_idx]) = seconds; 
                spike_ts_inds[unit_idx]++;

                // get waveform
                if (store_wf) {
                    pSpikeChannelHeader = &spikeChannels[dataBlock.Channel-1];
                    wf_gain[dataBlock.Channel] = max_spike_magnitude/(signed_adc_levels*spike_preamp_gain*pSpikeChannelHeader->Gain);
                    if (dataBlock.Unit >= 1) {
                        for(int i=0; i<dataBlock.NumberOfWordsInWaveform; i++)
                        {
                            if (verbose)
                                cout << "getting wf" << endl;
                            *spike_wf_inds[unit_idx] = buf[i];
                            spike_wf_inds[unit_idx]++;
                        }
                    }
                }
            }
        } else if(dataBlock.Type == PL_ExtEventType) { 
            if (dataBlock.Channel == STROBED_CHANNEL) {
                if (verbose) {
                    printf("strobed\n"); 
                    cout << evt_data_inds[STROBED_CHANNEL] << endl; 
                }

                *(evt_data_inds[STROBED_CHANNEL]) = seconds; // set timestamp
                evt_data_inds[STROBED_CHANNEL]++; 
                *(evt_data_inds[STROBED_CHANNEL]) = (double) dataBlock.Unit; // set event 
                evt_data_inds[STROBED_CHANNEL]++; 
            }
        } else if(dataBlock.Type == PL_ADDataType) { 
            if (t_start == -1) 
                t_start = seconds;
            if ( AD_channels_to_import[dataBlock.Channel+1] ) {
                if (verbose)
                    printf("AD Channel--Channel: %d\n", dataBlock.Channel); 
                pSlowChannel = &slowChannels[dataBlock.Channel];
                int gain = pSlowChannel->Gain;
                int preamp_gain = pSlowChannel->PreAmpGain;
                for(int i=0; i < dataBlock.NumberOfWordsInWaveform; i++)
                {
                    n_samples = fileHeader.EVCounts[dataBlock.Channel+300];
                    if (n_samples > 0) {
                        *(AD_inds[dataBlock.Channel]) = (buf[i]*fileHeader.SlowMaxMagnitudeMV)/(0.5*pow(2,fileHeader.BitsPerSlowSample)*gain*preamp_gain);
                        AD_inds[dataBlock.Channel]++;
                    }
                }
            }
        }
    }

    PyDict_SetItemString(data, "t_start", PyFloat_FromDouble(t_start));
    if (verbose)
        cout << "finished processing file" << endl;

    free(ch);
    free(sig_name);
    free(AD_inds);
    free(wf_name);
    free(evt_data_inds);
    free(spike_wf_inds);
    free(spike_ts_inds);
    return data;
}


// 
//  Distutils stuff
//
static PyMethodDef methods[] = {
    {"read_plx_headers", py_read_plx_headers, METH_VARARGS, 
        "Creates a dictionary of the plx headers"},
    {"import_file", (PyCFunction) import_file, METH_VARARGS|METH_KEYWORDS, 
        "get all the components of the *_raw.mat"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC 
initplxread(void)
{
    (void) Py_InitModule("plxread", methods);   
    import_array();
}
