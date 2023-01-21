
/********************************************************************
 *                                                                  *
 * THE SOURCE CODE IS (C) COPYRIGHT                                 *
 * by the University of Liege http://www.montefiore.ulg.ac.be       *
 *                                                                  *
 ********************************************************************

*********************************************************************/

#include "convolution_plugin.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ipps.h"
#include "ippcore.h"
#include <stdio.h>
#include "wave_plugin.h"
#include <omp.h>

GThreadPool* pool;

static freqBlocConvoFunctions *fbc = &fbc_ipp;

int getConvoBlockAmount(ConvoBlockInfo *cbi, int maxBlockAmount)
{
	int maxbc;

	int RIRblockAmount = cbi->iri->RIRblockAmount;

	maxbc = ((maxBlockAmount < RIRblockAmount) ?  maxBlockAmount : RIRblockAmount);

	cbi->iri->convoBlockAmount = maxbc;

	return 0;
}


/*decoupage de ir en blocs, et calcul de la fft de
  chaque bloc. Cette  fonction est appelee par initImpResp et changeIR*/
static int initBlockIr(ImpRespInfo *iri, float *ir, int lengthIr)
{
  int i;

  int blockSize = iri->blockSize;

  int RIRblockAmount = 0;

  iri->lengthIr = lengthIr;

  RIRblockAmount = (int) ceil((float) iri->lengthIr / (float) iri->blockSize);

  iri->RIRblockAmount = RIRblockAmount;

  if (RIRblockAmount == 0)
  {
    fprintf(stdout, "fatal error in initBlockIr\n");

    exit(EXIT_FAILURE);
  }

  /*allocation de la memoire pour la TF de la RI*/
  iri->irFreq = (float *) fbc->malloc(iri->blockSize * 2 * RIRblockAmount);

  /* calcul de la FFT de chaque bloc de la RI */
  for (i = 0; i < (iri->RIRblockAmount - 1); i++)
  {
    memcpy(iri->irFreq + i * 2 * blockSize, ir + i * blockSize, sizeof(float) * blockSize);

    fbc->fftForward(iri->fftLook, iri->irFreq + (i * 2 * blockSize));

  }

  /* pour le dernier bloc (qui est incopmplet).
     ... mais ce n'est peut etre pas necessaire... 
     voir specification fonction memcpy*/

  i = (iri->RIRblockAmount - 1);

  memcpy(iri->irFreq + i * 2 * blockSize, ir + i * blockSize, sizeof(float) * (lengthIr % blockSize));

  fbc->fftForward(iri->fftLook, iri->irFreq + (i * 2 * blockSize));

  return 0;
}

/*Init impulse response info (mono)*/
int initImpResp(ImpRespInfo *iri, void *fftLook, float *ir, int lengthIr, int blockSize, int ch)
{
  iri->fftLook = fftLook;

  iri->blockSize = blockSize;

  iri->ch = ch;

  initBlockIr(iri, ir, lengthIr);

  return 0;
}

/*Clear impulse response info (mono)*/
//int clearImpResp(ImpRespInfo *iri)
//{
//  if(iri->irFreq != NULL) g_free(iri->irFreq);
//
//  return 0;
//}

/* Init convolution block info for one anechoic source (mono), given one impulse response (mono) */
int initBlockConvol(ConvoBlockInfo *cbi, ImpRespInfo *iri)
{
  int i;

  cbi->iri = iri;

  cbi->normalisation = 1; /** @todo ceci doit etre un parametre !!! ou le calculer!!!!!!*/

  /*allocation de buffer pour le resultat*/
  cbi->buffer_freq = (float *) fbc->malloc(iri->RIRblockAmount * iri->blockSize * 2);

  cbi->ptr_buffer_freq = (float **) g_malloc(iri->RIRblockAmount * sizeof(float *));

  cbi->temporary_buffer = (float *) fbc->malloc(iri->blockSize * 2);
 
  for (i = 0; i < iri->RIRblockAmount; i++)
  {
    cbi->ptr_buffer_freq[i] = cbi->buffer_freq + i * iri->blockSize * 2;
  }

  cbi->buffer_tmp = (float *) fbc->malloc(iri->blockSize * 2);

  cbi->iri->convoBlockAmount = iri->RIRblockAmount;

  cbi->anechoFreq = (float *) fbc->malloc(2 * iri->blockSize);

  cbi->mutexProcessing = g_mutex_new();

  return 0;
}

/* Clear convolution block info for one anechoic source (mono) */
int clearBlockConvol(ConvoBlockInfo *cbi)
{
	g_mutex_lock(cbi->mutexProcessing); 
	g_mutex_unlock(cbi->mutexProcessing); 
	g_mutex_free(cbi->mutexProcessing);

	fbc->free(cbi->buffer_freq);

	fbc->free(cbi->temporary_buffer);

	fbc->free(cbi->buffer_tmp);

	g_free(cbi->ptr_buffer_freq);

  return 0;
}

/* Update convolution block info for one anechoic source (mono), given one NEW impulse response (mono) */
int changeBlockConvol(ConvoBlockInfo *cbi, ImpRespInfo *newIri)
{
  ImpRespInfo *oldIri = cbi->iri;    

  float *oldbufferFreq;

  float **oldptr_buffer_freq;

  int i, j;

  float *newBuffer_freq;

  float **newPtr_buffer_freq;

  int minRIRblockAmount;

  if (oldIri->blockSize != newIri->blockSize)
  {
    fprintf(stdout, "ERROR: block size of new and old ir must be the same!\n");

    exit(EXIT_FAILURE);
  }

  if(newIri->RIRblockAmount == oldIri->RIRblockAmount)
  {
	  g_mutex_lock(cbi->mutexProcessing); 

	  cbi->iri = newIri; 

	  g_mutex_unlock(cbi->mutexProcessing);
  }
  else
  {
	   oldbufferFreq = cbi->buffer_freq;

	   oldptr_buffer_freq = cbi->ptr_buffer_freq;

	   /*allocation de buffer pour le résultat*/
	  newBuffer_freq = (float *) fbc->malloc(newIri->RIRblockAmount * newIri->blockSize * 2);

	  newPtr_buffer_freq = (float **) g_malloc(newIri->RIRblockAmount * sizeof(float *));

	  minRIRblockAmount = (newIri->RIRblockAmount < oldIri->RIRblockAmount) ? newIri->RIRblockAmount : oldIri->RIRblockAmount;

	  for (i = 0; i < newIri->RIRblockAmount; i++)
	  {
		newPtr_buffer_freq[i] = newBuffer_freq + i * newIri->blockSize * 2;
	  }

	  for (i = 0; i < minRIRblockAmount; i++)
	  {
		for (j = 0; j < 2 * oldIri->blockSize; j++)
		{
		  newPtr_buffer_freq[i][j] += cbi->ptr_buffer_freq[i][j];
		}
	  }

	  g_mutex_lock(cbi->mutexProcessing); 

	  cbi->buffer_freq = newBuffer_freq;  

	  cbi->ptr_buffer_freq = newPtr_buffer_freq;   

	  cbi->iri = newIri; 

	  cbi->iri->convoBlockAmount = cbi->iri->RIRblockAmount;

	  g_mutex_unlock(cbi->mutexProcessing);

	  fbc->free(oldbufferFreq);

	  g_free(oldptr_buffer_freq);
  }

  if(oldIri->irFreq != NULL)
  {
	  fbc->free(oldIri->irFreq);

	  fbc->fftClear(oldIri->fftLook);
  }
 
  return 0;
}

/*Multiply the complex numbers xy1 and xy2 and add the result to the complex number r*/
int multComplexOvAdd(float *r, float *xy1, float *xy2)
{
  r[0] += (xy1[0] * xy2[0]) - (xy1[1] * xy2[1]);

  r[1] += (xy1[0] * xy2[1]) + (xy1[1] * xy2[0]);

  return 0;
}

/*Multiply the 2 complex vectors (blockInFreq and hFr) of n/2 complex numbers, 
  and add the result to the vector of complex numbers (resultFreq). 
  WARNING: VECTORS HAVE SPECIAL ORDER, SEE smallft.c
  (convolution is a multiplication in the frequency domain)
*/
int convolFreqOvAdd(float *blockInFreq, float *hFr, float *resultFreq , float *temp, int n)
{
  int i, imax;

  imax = n/2 - 1;

  /* WARNING: we must be take into acount the order of frq cmpt (see "smallft.c")*/

  resultFreq[0] += blockInFreq[0] * hFr[0];

  resultFreq[n-1] += blockInFreq[n-1] * hFr[n-1];

  /*#pragma omp parallel
  {

	   #pragma omp for nowait private(i)*/

  for (i = 0; i < imax; i++)
  {
    //multComplexOvAdd(resultFreq + i*2 + 1, blockInFreq + i*2 + 1, hFr + i*2 + 1);
  
	resultFreq[i*2 +1] += (blockInFreq[i*2 + 1] * hFr[i*2 + 1]) - (blockInFreq[i*2 + 2] * hFr[i*2 + 2]);

    resultFreq[i*2 + 2] += (blockInFreq[i*2 + 1] * hFr[i*2 + 2]) + (blockInFreq[i*2 + 2] * hFr[i*2 + 1]);
  
  }
  /*}*/

  return 0;
}

void *blockConvol(void *Info)
{
	ConvoBlockInfo *cbi = (ConvoBlockInfo *) Info;

	int i;

	int blockSize = cbi->iri->blockSize;

	int RIRblockAmount = cbi->iri->RIRblockAmount;

	int norm = cbi->normalisation * 2 * blockSize; /*WARNING : this is not very good!!*/

	float *ptrtmp;

	if(cbi->nin > 0)
	{
		for (i = 0; i < cbi->iri->convoBlockAmount; i++)
		{
			fbc->mulFreqAdd(cbi->anechoFreq, cbi->iri->irFreq + (2 * blockSize * i),cbi->ptr_buffer_freq[i], cbi->temporary_buffer, 2 * blockSize);
		}
	}

	fbc->fftBackward(cbi->iri->fftLook, cbi->ptr_buffer_freq[0]);

	ippsDivC_32f_I(norm, cbi->ptr_buffer_freq[0], blockSize * 2 - 1);

	ippsAdd_32f_I(cbi->ptr_buffer_freq[0], cbi->buffer_tmp, blockSize * 2 - 1);

	g_mutex_lock(cbi->mutex);

	for (i = 0; i < blockSize; i++)
	{
		cbi->outpcm[cbi->maxOutputChannels*i +cbi->chID] += cbi->buffer_tmp[i];
	}

	g_mutex_unlock(cbi->mutex);

	ippsMove_32f(cbi->buffer_tmp + blockSize, cbi->buffer_tmp, blockSize);

	ippsZero_32f(cbi->buffer_tmp + blockSize, blockSize);

	ptrtmp = cbi->ptr_buffer_freq[0];

	ippsMove_32f(cbi->ptr_buffer_freq + 1, cbi->ptr_buffer_freq, RIRblockAmount - 1);

	cbi->ptr_buffer_freq[(RIRblockAmount - 1)] = ptrtmp;

	ippsZero_32f(cbi->ptr_buffer_freq[(RIRblockAmount - 1)], 2 * blockSize);

	cbi->processed = 1;

	g_mutex_unlock(cbi->mutexProcessing); 

	return 0;
}

/*Init common info for all multi IR and for all anechoic sources*/
MultiConvoGeneralInfo *initMultiConvoGeneral(int blockSize, int ch)
{
	MultiConvoGeneralInfo *mcgi = (MultiConvoGeneralInfo *) g_malloc(sizeof(MultiConvoGeneralInfo)); 

	mcgi->blockSize = blockSize;

	mcgi->ch = ch;

	fbc->fftInit(&mcgi->fftLook, blockSize * 2);

  return mcgi;
}

/*Clear common info for all multi IR and for all anechoic sources*/
int clearMultiConvoGeneral(MultiConvoGeneralInfo *mcgi)
{
  fbc->fftClear(mcgi->fftLook); /* clear fft*/

  return 0;
}

/* Init info for one multichannel IR, for all anechoic sources
   ir is the multichannel impulse response buffer, with lengthIr samples*/
MultiIR *initMultiIr(int blockSize, int ch, float **ir, int lengthIR)
{
	int i;

	MultiConvoGeneralInfo *mcgi;

	MultiIR *mir = (MultiIR *) g_malloc(sizeof(MultiIR));

	mir->mcgi = initMultiConvoGeneral( blockSize, ch); 

	mir->iri = (ImpRespInfo **) g_malloc0(ch* sizeof(ImpRespInfo *));

	if(lengthIR == 1)
	{			
		for (i = 0; i < ch; i++)
		{
			ir[i] = (float *) g_malloc(sizeof(float)); 
			
			ir[i][0] = 1;
		}
		lengthIR = 1;
	}

	for (i = 0; i < ch; i++)
	{
		mir->iri[i] = (ImpRespInfo *) g_malloc0(sizeof(ImpRespInfo));

		fbc->fftInit(&mir->iri[i]->fftLook, blockSize * 2);

		initImpResp(mir->iri[i], mir->iri[i]->fftLook, ir[i], lengthIR, blockSize, i);
	}

	return mir;
}


MultiIR *initMultiIr2(int blockSize, int ch, float **ir, int lengthIR)
{
	int i;

	MultiConvoGeneralInfo *mcgi;

	MultiIR *mir = (MultiIR *) g_malloc(sizeof(MultiIR));

	mir->iri = (ImpRespInfo **) g_malloc0(ch* sizeof(ImpRespInfo *));

	for (i = 0; i < ch; i++)
	{
		mir->iri[i] = (ImpRespInfo *) g_malloc0(sizeof(ImpRespInfo));

		fbc->fftInit(&mir->iri[i]->fftLook, blockSize * 2);

		initImpResp(mir->iri[i], mir->iri[i]->fftLook, ir[i], lengthIR, blockSize, i);
	}

	return mir;
}

/*Update convolution block info for one anechoic source, given one NEW multichannel impulse response */
int changeMultiIr(void *Info, MultiIR *newIr)
{
	int i, ch;

	MultiConvoBlockInfo *mcbi = (MultiConvoBlockInfo *) Info;

	ch = newIr->mcgi->ch;

	//g_mutex_lock(mcbi->mutexProcessing); 

	for (i = 0; i < ch; i++)
	{
		ImpRespInfo *iri = newIr->iri[i];

		ConvoBlockInfo *cbi = mcbi->cbi[i];

		changeBlockConvol(cbi, iri);
	}
	//g_mutex_unlock(mcbi->mutexProcessing); 

	return 0;
}

/*Init convolution block info for one source, given one multichannel IR*/
void *initMultiConvo(int blockSize, int ch, int* listchIR, int maxOutputChannels)
{
	MultiConvoBlockInfo *mcbi =  (MultiConvoBlockInfo *) g_malloc0(1* sizeof(MultiConvoBlockInfo));

	int lengthIR = 1;

	float **ir = (float **) g_malloc(sizeof(float *) * ch);

	MultiIR *mir = initMultiIr(blockSize, ch, ir, lengthIR);

	int i;

	mcbi->ncbi = 0;

	mcbi->mcgi = mir->mcgi;

	mcbi->cbi = (ConvoBlockInfo **) g_malloc0(ch* sizeof(ConvoBlockInfo *));

	mcbi->anechoFreq = (float *) fbc->malloc(2 * blockSize);

	ippsZero_32f(mcbi->anechoFreq, 2 * blockSize);

	mcbi->maxBlock = 1000;

	mcbi->listchIR = listchIR;

	mcbi->maxOutputChannels =maxOutputChannels;

	for (i = 0; i < ch; i++)
	{ 
		mcbi->cbi[i] = (ConvoBlockInfo *) g_malloc0(1* sizeof(ConvoBlockInfo));

		initBlockConvol(mcbi->cbi[i], mir->iri[i]); /* allocate memory */

		mcbi->ncbi += 1;
	}

	return mcbi;
}

int changeMultiConvo(MultiConvoBlockInfo *mcbi, int blockSize, int ch, int* listchIR, float **ir, int lengthIR)
{
	int i, oldch;

	MultiIR *newIr;

	ImpRespInfo *iri;

	ConvoBlockInfo *cbi;

	oldch = mcbi->mcgi->ch;
	

	if(oldch< ch)
	{
		newIr = initMultiIr(blockSize, ch, ir, lengthIR);

		mcbi->mcgi = newIr->mcgi;

		mcbi->listchIR = listchIR;

		for (i = 0; i < oldch; i++)
		{
			iri = newIr->iri[i];

			cbi = mcbi->cbi[i];

			changeBlockConvol(cbi, iri);
		}

		mcbi->cbi = (ConvoBlockInfo **) g_realloc(mcbi->cbi, ch* sizeof(ConvoBlockInfo *));

		for (i = oldch; i < ch; i++)
		{
			iri = newIr->iri[i];

			cbi = mcbi->cbi[i];

			mcbi->cbi[i] = (ConvoBlockInfo *) g_malloc0(sizeof(ConvoBlockInfo));

			initBlockConvol(mcbi->cbi[i], iri);

			mcbi->ncbi += 1;
		}

	}
	else if(oldch > ch)
	{
		newIr = initMultiIr(blockSize, ch, ir, lengthIR);

		mcbi->mcgi = newIr->mcgi;

		mcbi->listchIR = listchIR;

		for (i = 0; i < ch; i++)
		{
			iri = newIr->iri[i];

			cbi = mcbi->cbi[i];

			changeBlockConvol(cbi, iri);
		}

		for (i = ch; i < oldch; i++)
		{
			g_mutex_lock(cbi->mutexProcessing); 
			mcbi->ncbi -= 1;
			g_mutex_unlock(cbi->mutexProcessing); 

			clearBlockConvol(mcbi->cbi[i]);
		}

		mcbi->cbi = (ConvoBlockInfo **) g_realloc(mcbi->cbi, ch* sizeof(ConvoBlockInfo *));
	}
	else 
	{
		newIr = initMultiIr2(blockSize, ch, ir, lengthIR);

		for (i = 0; i < ch; i++)
		{
			iri = newIr->iri[i];

			cbi = mcbi->cbi[i];

			changeBlockConvol(cbi, iri);
		}

		g_free(newIr);
	}

	return 0;
}
/*Clear convolution block info for one source and one multichannel IR*/
int clearMultiConvo (void *Info)
{
	int i, ch;

	MultiConvoBlockInfo *mcbi;

	if(Info == NULL) return 1;

	mcbi = (MultiConvoBlockInfo *) Info;

	ch = mcbi->mcgi->ch;

	fbc->free(mcbi->anechoFreq);

	//g_free(mcbi->anechoFreq);

	//g_free(mcbi->listchIR);

	for (i = 0; i < ch; i++)
	{
		if (mcbi->cbi + i != NULL) clearBlockConvol(mcbi->cbi[i]);
	}

	if (mcbi->cbi != NULL) g_free(mcbi->cbi);

	clearMultiConvoGeneral(mcbi->mcgi);

	return 0;
}

int multiConvolution(float *pcmAnecho, int nin, float *outpcm, MultiConvoBlockInfo *mcbi, int maxBlockAmount, GMutex** mutex)
{
	int i, j, ch, blockSize, ncbi;

	ch = mcbi->mcgi->ch;

	blockSize = mcbi->mcgi->blockSize;

	if (nin > blockSize)
	{
		fprintf(stdout, "fatal error multiConvolution: nin > blockSize (%d > %d)\n", nin, blockSize);

		exit(EXIT_FAILURE);
	}

	ippsZero_32f(mcbi->anechoFreq, 2 * blockSize);

	ippsCopy_32f(pcmAnecho, mcbi->anechoFreq, nin);

	fbc_ipp.fftForward(mcbi->mcgi->fftLook, mcbi->anechoFreq);

	ncbi = mcbi->ncbi;

	for (i = 0; i < ncbi; i++)
	{
		ConvoBlockInfo *cbi = mcbi->cbi[i];

		g_mutex_lock(cbi->mutexProcessing); 

		ippsCopy_32f(mcbi->anechoFreq, cbi->anechoFreq, 2 * blockSize);

		cbi->chID = mcbi->listchIR[i];

		cbi->mutex = mutex[cbi->chID];

		cbi->outpcm = outpcm;

		cbi->nin = nin;

		cbi->maxOutputChannels = mcbi->maxOutputChannels;

		getConvoBlockAmount(cbi, mcbi->maxBlock);

		cbi->processed = 0;

		g_thread_pool_push(pool, cbi, NULL);
	}

	for (i = 0; i < ncbi; i++)
	{
		while(mcbi->cbi[i]->processed == 0)
		{
			g_usleep((gulong) (1000));
		}
	}

	return 0;
}



int noConvolution(float *pcmAnecho, int nin, float *outpcm, MultiConvoBlockInfo *mcbi, GMutex** mutex)
{
	int i, j;

	for(i = 0; i < mcbi->mcgi->ch; i++)
	{
		g_mutex_lock(mutex[i]);
		for(j = 0; j < nin; j++)
		{
			outpcm[mcbi->maxOutputChannels*j +i] += pcmAnecho[j];
		}
		g_mutex_unlock(mutex[i]);
	}

	return 0;
}

int initConvolutionSession()
{
	int numberOfCores, numberOfThreads = 0;

	ippGetNumThreads(&numberOfCores);

	switch (numberOfCores)
	{

		case 1:

			numberOfThreads = 1;

			break;

		case 2:

			numberOfThreads = 1;

			break;

		case 3:

			numberOfThreads = 2;

			break;

		case 4:

			numberOfThreads = 3;

			break;

		case 8:

			numberOfThreads = 5;

			break;

		case 16:

			numberOfThreads = 12;

			break;

		default :
				
			break;

	}

	if(numberOfThreads == 0) numberOfThreads = numberOfCores - 1;

	pool = g_thread_pool_new (blockConvol, NULL, numberOfThreads, TRUE, NULL);

	return 0;
}

int closeConvolutionSession()
{

	g_thread_pool_free (pool, TRUE, FALSE);

	return 0;
}