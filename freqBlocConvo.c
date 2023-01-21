#include "convolution_plugin.h"
#include "stdio.h"
#include <math.h>
#include "ipps.h"
#include "ippcore.h"
typedef struct ipp_drft_lookup
{
	int n;

	IppsFFTSpec_R_32f* pFFTSpec;

	Ipp8u *pBuffer;

} ipp_drft_lookup;

void* mall_ogg(int len)
{
	float* data;

	data = g_malloc0(len*sizeof(float));

	return data;
}

void Init(void** infos, int n)
{
	*infos = (drft_lookup*) g_malloc(sizeof(drft_lookup));

	FFT_OggVorbis.init(*infos, n);
}

void Clear(void* infos)
{
	infos = (drft_lookup*) infos;

	FFT_OggVorbis.clear(infos);
}

void Forward(void* infos,float *data)
{
	infos = (drft_lookup*) infos;

	FFT_OggVorbis.forward(infos, data);
}
void Backward(void* infos,float *data)
{
	infos = (drft_lookup*) infos;

	FFT_OggVorbis.backward(infos, data);
}



void* mall(int len)
{
	Ipp32f* data = ippsMalloc_32f(len);

	ippsZero_32f(data, len);

	return data;
}


void (fre)(void* data)
{
	ippsFree(data);
}

void ippFFTInit(void** infos2, int n)
{
	IppStatus status;

	int BufferSize, order = 0, poww = 0;

	ipp_drft_lookup *infos = (ipp_drft_lookup*) ippMalloc(sizeof(ipp_drft_lookup));

	while(poww < n)
	{
		poww = pow(2, order);

		order++;
	}

	order -=1;

	if(n != poww)
	{
		fprintf(stdout,"\nError in audio settings :\n ->audio buffer block size must a power of 2\n");
	}

	status = ippsFFTInitAlloc_R_32f(&(infos->pFFTSpec),
		order,
		IPP_FFT_NODIV_BY_ANY /*IPP_FFT_DIV_FWD_BY_N IPP_FFT_NODIV_BY_ANY IPP_FFT_DIV_INV_BY_N*/,
		ippAlgHintFast /*,ippAlgHintNone */);

	if( status != ippStsNoErr)
	{
		fprintf(stdout, "\nError in ftt initialisation\n");
	}

	status = ippsFFTGetBufSize_R_32f(infos->pFFTSpec, &BufferSize);

	if( status != ippStsNoErr)
	{
		fprintf(stdout, "\nError in ftt buffer size\n");
	}

	infos->pBuffer = ippsMalloc_8u(BufferSize);

	infos->n = n;

	*infos2 = infos;

}

void ippFFTClear(void* infos2)
{
	IppStatus status;

	ipp_drft_lookup* infos = (ipp_drft_lookup*) infos2;

	status = ippsFFTFree_R_32f(infos->pFFTSpec);

	if( status != ippStsNoErr)
	{
		fprintf(stdout, "\nError in ftt termination\n");
	}

	ippsFree(infos->pBuffer);

	ippFree(infos);
}

void ippFFTForward(void* infos2,float *data)
{
	IppStatus status;

	ipp_drft_lookup* infos = (ipp_drft_lookup*) infos2;

	status = ippsFFTFwd_RToPack_32f_I(data, infos->pFFTSpec, infos->pBuffer);

	if( status != ippStsNoErr)
	{
		fprintf(stdout, "\nError in ftt forward\n");
	}
}

void ippFFTBackward(void* infos2,float *data)
{
	IppStatus status;

	ipp_drft_lookup* infos = (ipp_drft_lookup*) infos2;

	status = ippsFFTInv_PackToR_32f_I(data, infos->pFFTSpec, infos->pBuffer);

	if( status != ippStsNoErr)
	{
		fprintf(stdout, "\nError in ftt backward\n");
	}
}



int addVectors2(float *acc, float *v2, int n)
{
  int i;

  for (i = 0; i < n; i++) {acc[i] += v2[i];}

  return 0;
}

int ippFFTmulfreqadd(float *src1, float *src2,float *dst, float *temp, int len)
{
	IppStatus status;

	status = ippsMulPack_32f(src1, src2, temp, len);

	if( status != ippStsNoErr)
	{
		fprintf(stdout, "\nError in mulfreqadd function\n");

		return 1;
	}

	status = ippsAdd_32f_I(temp, dst, len);

	return 0;
}


freqBlocConvoFunctions fbc_ipp =
{
	&mall, 

	&fre,

	&ippFFTInit,  

	&ippFFTClear, 

	&ippFFTForward,

	&ippFFTBackward,

	&ippFFTmulfreqadd,

};

freqBlocConvoFunctions fbc_ogg =
{
	&mall_ogg,

	&g_free,

	&Init,

	&Clear,

	&Forward,

	&Backward,

	&convolFreqOvAdd,

};
