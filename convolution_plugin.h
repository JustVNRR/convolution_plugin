/********************************************************************
 *                                                                  *
 * THE SOURCE CODE IS (C) COPYRIGHT                                 *
 * by the University of Liege http://www.montefiore.ulg.ac.be       *
 *                                                                  *
 ********************************************************************

 function:

 last mod: $Id: convolution.h 463 2006-07-10 22:46:30Z NicolasW $

*********************************************************************/


#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_

#include "fft_plugin.h"
#include "glib.h"



#ifdef BUILD_BLOCK_CONVO
    #define BLOCK_CONVO_IE __declspec(dllexport)
#else
    #define BLOCK_CONVO_IE __declspec(dllimport)
#endif

/**@defgroup convolution_plugin Module de convolution
 * \brief Fast convolution functions. Overlap-add with subdivision of impulse responses into blocks
 * \section depconv Dépendances
 * - @link fft plugin FFT@endlink
 * - @link wav plugin de gestion des fichiers ".wav"@endlink
 *  @{
 */

/*!\brief Impulse response (mono) */
typedef struct ImpRespInfo
{  
	int blockSize;

	int RIRblockAmount;

	int convoBlockAmount;

	float *irFreq;

	int lengthIr;

	int ch;

	void* fftLook;

} ImpRespInfo;

/*!\brief Common info for all multi IR and for all anechoic sources */
typedef struct MultiConvoGeneralInfo
{
	int blockSize;

	int ch;

	void* fftLook;

} MultiConvoGeneralInfo;

/*!\brief Info for one multichannel IR, for all anechoic sources */
typedef struct MultiIR
{
	MultiConvoGeneralInfo *mcgi;

	ImpRespInfo **iri;

	char irFileName[256];

} MultiIR;

/*!\brief Convolution block info for one anechoic source (mono) */
typedef struct ConvoBlockInfo
{
	ImpRespInfo *iri;/*pointer on ir mono*/

	float *buffer_freq;

	float *temporary_buffer;

	float **ptr_buffer_freq;

	float *buffer_tmp;

	int normalisation;

	int nin;

	float *anechoFreq;

	float *outpcm;

	GThread *convoThread;

	int chID;

	int maxOutputChannels;

	GMutex* mutex;

	GMutex *mutexProcessing;

	int processed;
  
} ConvoBlockInfo;

/*!\brief Convolution block info for one source and one multichannel IR*/
typedef struct MultiConvoBlockInfo
{
	MultiConvoGeneralInfo *mcgi;

	ConvoBlockInfo **cbi;

	float *anechoFreq;

	int maxBlock;

	int* listchIR;

	int maxOutputChannels;

	int ncbi;

} MultiConvoBlockInfo;

BLOCK_CONVO_IE MultiIR *initMultiIr(int blockSize, int ch, float **ir, int lengthIR);

BLOCK_CONVO_IE MultiIR *initMultiIr2(int blockSize, int ch, float **ir, int lengthIR);

BLOCK_CONVO_IE int multiConvolution(float *pcmAnecho, int nin, float *outpcm, MultiConvoBlockInfo *mcbi, int maxBlockAmount, GMutex** mutex);

BLOCK_CONVO_IE int noConvolution(float *pcmAnecho, int nin, float *outpcm, MultiConvoBlockInfo *mcbi, GMutex** mutex);

BLOCK_CONVO_IE int changeMultiIr(void *MultiConvoBlockInfo, MultiIR * newIR);

BLOCK_CONVO_IE void *initMultiConvo(int blockSize, int ch, int* listOfChannels, int maxOutPutChannels);

BLOCK_CONVO_IE int clearMultiConvo (void *Info);

BLOCK_CONVO_IE int changeMultiConvo(MultiConvoBlockInfo *mcbi, int blockSize, int ch, int* chList, float **ir, int lengthIR);

BLOCK_CONVO_IE MultiIR* initMultiIrFile(char *irFileName);

BLOCK_CONVO_IE MultiIR* initMultiIrFile2(void *mcgi2, char *irFileName);

BLOCK_CONVO_IE MultiConvoGeneralInfo *initMultiConvoGeneral(int blockSize, int ch);

BLOCK_CONVO_IE int convolFreqOvAdd(float *blockInFreq, float *hFr, float *resultFreq, int n);

BLOCK_CONVO_IE void *blockConvol(void *Info);

BLOCK_CONVO_IE int getConvoBlockAmount(ConvoBlockInfo *cbi, int maxBlockAmount);

BLOCK_CONVO_IE int initConvolutionSession();

BLOCK_CONVO_IE int closeConvolutionSession();

typedef struct freqBlocConvoFunctions
{
	void* (*malloc)(int len);

	void (*free)(void*);

	void (*fftInit)(void**, int n);

	void (*fftClear)(void*);

	void (*fftForward)(void*,float *data);

	void (*fftBackward)(void*,float *data);

	int(*mulFreqAdd)(float *src1, float *src2,float *dst, float *temp, int len);

} freqBlocConvoFunctions;

BLOCK_CONVO_IE freqBlocConvoFunctions fbc_ogg;

BLOCK_CONVO_IE freqBlocConvoFunctions fbc_ipp;

/** @} */

#endif
