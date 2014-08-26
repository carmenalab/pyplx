// Modification instructions from original library
//------------------------------------------------
//
// remove #include <windows.h>
// add #include <stdlib.h>
// add #include <string.h>
// #include "..\Plexon.h" -> #include "../Plexon.h"
// change return type of main function from "void" to "int"

#ifndef _PLXREAD_H
#define _PLXREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Plexon.h"
#include <iostream>
#include <stdint.h>
#include "Python.h"
#include "arrayobject.h"

#define LONGLONG int64_t 
#define MAX_SPIKE_CHANNELS   (128)
#define MAX_EVENT_CHANNELS   (512)
#define MAX_SLOW_CHANNELS    (256)

#define MAX_SAMPLES_PER_WAVEFORM (256)

#define STROBED_CHANNEL 257

#define MAX_SPIKES_PER_ELECTRODE 5

// constants for spike waveform conversion (avoid version number issues)
double max_spike_magnitude = 3000;
double signed_adc_levels = 2048;
double spike_preamp_gain = 1000;

// Open PLX file
FILE* fp = 0 ;

// PLX File header structure
PL_FileHeader fileHeader;

// PLX Spike Channel headers
PL_ChanHeader spikeChannels[MAX_SPIKE_CHANNELS];

// PLX Event Channel headers
PL_EventHeader eventChannels[MAX_EVENT_CHANNELS];

// PLX Slow A/D Channel headers
PL_SlowChannelHeader slowChannels[MAX_SLOW_CHANNELS];

//
// position in file where data begins
int data_start = 0 ;

// Dump header information
bool bDumpHeader = false ;
// Dump spike data blocks
bool bDumpSpike = false ;
// Dump event data blocks
bool bDumpEvent = false ;
// Dump slow A/D data blocks
bool bDumpSlow  = false ;

// Maximum number of data blocks to dump
#define MAX_DATA_BLOCKS (5000)

PyObject* py_get_raw_data(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject* py_get_strobed(PyObject* self, PyObject* args);
PyObject* py_read_plx_headers(PyObject *self, PyObject *args);
PyObject* py_get_spike_ts(PyObject *self, PyObject *args, PyObject *keywds);

void get_spike_ts ( PyObject* spike_times );
int read_plx_headers(char* fname);
void get_events (PyObject* event_times);
void datablock_to_spike_ts(PyObject*, PyObject*, PL_DataBlockHeader, short*, bool); 
void datablock_to_events(PyObject*, PyObject*, PL_DataBlockHeader);
void datablock_to_AD_samples(PyObject*, PyObject*, PL_DataBlockHeader, short int*);
void datablock_to_AD_samples2(float*, PL_DataBlockHeader, short int*); 

void get_raw_data( PyObject* spike_times, PyObject* ad_chan_times, PyObject* event_times, bool ); 
PyArrayObject* pyarray_alloc(int data_len);
using namespace std;


#endif
