
/***************************************************************************************
*** pcmx  - by Amish S. Dave (8/23/96)
***
*** TODO:
***
*** HISTORY:
***	  8/21/96
***		- actually, I wrote this program a while ago, and have been updating it
***			relatively frequently.  However, it has reached a very happy point now.
***		- pcmx now supports processing sampled data before writing it out.  It
***			builds a chain of various processes (I'll call them filters), which are
***			applied in the order given to each input entry.
***		- pcmx also allows the input entry and output entry to be specified.  Before
***			I was using a massive kludge with the -outentrystart parameter.  Now,
***			you give -entry <num> after the input or output file.  If after the input
***			filespec, it means to process only the specified entry.  If after the output
***			filespec, it means to start with the specified entry.
***		- I've implemented the "SCALE", "LOPASS", and "HIPASS" filters.  I ported
***			the code for filtering from Dan's VMS Fortran program UNITS.
***
***	  8/23/96
***		- I've implemented the "INVERT", "CUT" and "REVERSE" filters.  Now, pcmx can
***			be used to generate BOS and RBOS stimuli...  Pcmx now replaces many
***			little utility programs (like negpolarity, etc.) that I had lying
***			around.  The procedure for creating a stimulus:
***
***			   UX: aplot zy_hv12s_19960823a.pcm_seq2 		// to find good stimulus.  Determine <strt> and <stop>
***			   UX: pcmx -cut <strt> <stop> -hipass 500 -scale 96 zy_hv12s_19960823a.pcm_seq2 -entry 1 hv12.pcm_fix
***			   UX: rcp hv12.pcm_fix fred:
***			   UX: rlogin fred
***			   FR: signals
***             S> packin hv12
***             S> pkey				// answer question about saving time, etc. with 'Y'
***             S> outcreate zy_hv12
***             S> output
***             S> outclose
***             S> <press CTRL-Z>
***			   FR: copy zy_hv12.exp_idx $disk_fr1:[pcmlib.zf_sng]*
***			   FR: delete zy_hv12.exp_idx;
***
***			as can be seen, we still need signals because I can't write an .exp_idx file,
***			which is the only file format UNITS is really happy with.  The solution, I think,
***			is to modify units to accept .pcm_fix files.
***	
***	  8/29/96
***		- pcmx now copies the timestamp from the input entries, allowing, eg:
***			a repacking of many small .pcm_seq2 files into fewer .pcm_seq2 files
***			without losing the time of creation information, etc.
***		- there is a new option: '-noout', which forces pcmx to just read the input
***			files without creating an output file...  This is useful, for example, if
***			you just want to read the datestamp information, or learn the # of entries,
***			etc.
***
***	  9/2/96
***		- pcmx now copies the samplerate from the input entries.  I had forgotten to
***			do that before (and noticed it after adding the ability to write .wav files)
***		- minor cleanup to the -nooutput implementation.
***
***	  10/1/96
***		- pcmx now supports the option: '-pad <ms>', which allows taking a small segment
***			of data, and enlarging it by adding zeros to the end.  Mostly this is useful
***			for trying to plot a stimulus signal aligned with a PSTH or raster display...
***
***	  2/19/98
***		- pcmx now supports the option: '-rmsscale <dB>', which is like '-scale', but
***			sets the root-mean-square value to the desired value.  It warns of clipped
***			samples.
***		- '-scale' and '-rmsscale' now provide better formatted output, including RMS values.
***
***	  5/24/98
***		- pcmx now supports the '-dc' switch, which calculates and subtracts the DC offset
***		- pcmx now supports the '-rectify' switch, which rectifies the signal
***
***	  6/21/98
***		- pcmx now supports the '-sum' and '-avg' switches.
***
***	  7/6/98
***		- pcmx now supports resampling (using code adapted from rateconv.c,
***			Copyright (c) 1992, 1995 by Markus Mummert).
***
***   10/16/99
***		- pcmx now supports copying and pasting (cut doesn't work as expected - I need
***			to rename it to crop and add a real cut).  This is useful to, for example,
***			replace a bird song syllable with silence from elsewhere in the sound.
***
***   10/20/99
***		- pcmx now supports ramping at the beginning and/or end of a sound.
***			use 'linramp' or 'expramp' to ramp linearly or exponentially (linearly
***			in dB) respectively.  These commands require at least a single parameter
***			which is the range, in milliseconds, over which to ramp up.  If two
***			parameters are given, the first is the range of the onset while the
***			second is the offset range.  If only 1 is given, it will be applied
***			to both the onset and offset.  Give a zero to restrict the operation.
***			This is useful for removing clicks caused by a sudden transition at
***			the onset or offset of a stimulus.
***		- The simple addition of this feature turns pcmx into one of the better
***			sound players out there (well, okay - not really.  There's no support for
***			anything but mono sounds, so stereo wav files will only have one
***			channel played...)  But I can now use pcmx more easily from within
***			the festival text to speech convertor. (to amplify the annoyingly 
***			faint speech that program produces)
***
***   2/23/00
***		- minor change: scale, rmsscale, and agcscale now accept real values
***			for the dB amplitude
***************************************************************************************/


/******************************************************************************
** Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include "pcmio.h"
#include "version.h"


/******************************************************************************
** Definitions
*/
#define SUCCESS 1
#define ERROR 2

#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef PI
#define PI M_PI
#endif
#ifndef SQR
#define SQR(a)  ((a) * (a))
#endif

#define PROCESS_SCALE 1
#define PROCESS_LOPASS 2
#define PROCESS_HIPASS 3
#define PROCESS_CUT 4
#define PROCESS_INVERT 5
#define PROCESS_REVERSE 6
#define PROCESS_PAD 7
#define PROCESS_PREPEND 8
#define PROCESS_RMSSCALE 9
#define PROCESS_SMOOTH 10
#define PROCESS_AGC 11
#define PROCESS_STRIPDC 12
#define PROCESS_RECTIFY 13
#define PROCESS_SUM 14
#define PROCESS_AVG 15
#define PROCESS_RESAMPLE 16
#define PROCESS_SEXCONV 17
#define PROCESS_COPY 18
#define PROCESS_PASTE 19
#define PROCESS_PASTEOVER 20
#define PROCESS_LINRAMP 21
#define PROCESS_EXPRAMP 22

#define MAXENTRIES 10000

/******************************************************************************
** Structure declarations
*/
typedef struct procchain_t
{
	struct procchain_t *next;
	int type;
	int intarg1;
	int intarg2;
	double fltarg1;
	double fltarg2;
} procchain_t;


/******************************************************************************
** Global variables (low and hi-pass filter coefficients
*/
static double c1[6][3], c2[6][3];
procchain_t *proclist = NULL;

static short *bufsamples = NULL; 		/* For copying and pasting */
static int bufnsamples = 0; 			/* For copying and pasting */


/******************************************************************************
** Prototypes
*/
static void usage(void);
static procchain_t *add_proc(int type);
static void invert_samples(short *samples, int nsamples);
static void smooth_samples(short *samples, int nsamples, int nsmooth);
static void prepend_samples(short **samples_p, int *nsamples_p, int samplerate, double prependms);
static void pad_samples(short **samples_p, int *nsamples_p, int samplerate, double desiredms);
static void resample_samples(short **samples_p, int *nsamples_p, int *samplerate_p, int up, int down);
static void copy_samples(short **samples_p, int *nsamples_p, int samplerate, double start, double stop);
static void paste_samples(short **samples_p, int *nsamples_p, int samplerate, double start);
static void pasteover_samples(short **samples_p, int *nsamples_p, int samplerate, double start, double stop);
static void cut_samples(short **samples_p, int *nsamples_p, int samplerate, double start, double stop);
static void reverse_samples(short *samples, int nsamples);
static void rectify_samples(short *samples, int nsamples);
static void sum_samples(short *samples, int nsamples, int lasttime);
static void avg_samples(short *samples, int nsamples, int lasttime);
static void stripDC_samples(short *samples, int nsamples);
static void sexconv_samples(short *samples, int nsamples);
static void scale_samples(short *samples, int nsamples, double scaledB);
static void rmsscale_samples(short *samples, int nsamples, double scaledB);
static void linramp_samples(short *samples, int nsamples, int samplerate, double onramp, double offramp);
static void expramp_samples(short *samples, int nsamples, int samplerate, double onramp, double offramp);
static void agcscale_samples(short *samples, int nsamples, double scaledB, int winsamples);
static void lpfilt_samples(short *samples, int nsamples, int freq);
static void hpfilt_samples(short *samples, int nsamples, int freq);
static void lpcoef(double fc, double *tg_p);
static void hpcoef(double fc, double *tg_p);
static void hphelp(double ccc1, double ccc2, int isgn, double alp, double *cc1_p, double *cc2_p);
extern int rateconv(double fsi, double p_up, double p_down, short *src, int n_in, short *dst);


static procchain_t *add_proc(int type)
{
	procchain_t *nextproc;

	nextproc = proclist;
	if (nextproc != NULL) while (nextproc->next != NULL) nextproc = nextproc->next;
	if (nextproc == NULL) { proclist = (procchain_t *)malloc(sizeof(procchain_t)); nextproc = proclist; }
	else { nextproc->next = (procchain_t *)malloc(sizeof(procchain_t)); nextproc = nextproc->next; }
	nextproc->next = NULL;
	nextproc->type = type;
	return nextproc;
}

static void usage(void)
{
	fprintf(stdout, "\nUSAGE: pcmx accepts the following arguments:\n");
	fprintf(stdout, "  pcmx [-malloc] [-scale dB] [-hipass Hz] [-lopass Hz] [-invert] [-cut ms ms]\n");
	fprintf(stdout, "       [-reverse] <in-file-name> [-entry #] <out-file-name> [-entry #] [-noout]\n");
	fprintf(stdout, "       [-h|-help] [-v|-version]\n\n");
	fprintf(stdout, "    [-h|-help]     = prints this help message and then exits\n");
	fprintf(stdout, "    [-v|-version]  = prints the version and then exits\n");
	fprintf(stdout, "    [-malloc]      = optional, might speed things up at cost of memory\n");
	fprintf(stdout, "    [-scale dB]    = optional, can scale signal to, eg. 96dB (peak)\n");
	fprintf(stdout, "    [-rmsscale dB] = optional, can scale signal to, eg. 96dB (rms)\n");
	fprintf(stdout, "    [-linramp ms ms]= optional, ramps signal linearly over first ms samples\n");
	fprintf(stdout, "                     (note: 1st ms is on-ramp, 2nd is off-ramp.  If only 1st\n");
	fprintf(stdout, "                     number is specified, then it is used for both.)\n");
	fprintf(stdout, "    [-expramp ms ms]= optional, like linramp, but exponential (linear in dB)\n");
	fprintf(stdout, "    [-hipass Hz]   = optional, 4th order elliptical hipass filter\n");
	fprintf(stdout, "    [-lopass Hz]   = optional, 4th order Butterworth lopass filter\n");
	fprintf(stdout, "    [-invert]      = optional, inverts the polarity of signal (mults by -1)\n");
	fprintf(stdout, "    [-sexconv]     = optional, flips the byte sex of the samples\n");
	fprintf(stdout, "    [-prepend ms]  = optional, adds ms zeros to start of output\n");
	fprintf(stdout, "    [-pad ms]      = optional, adds zeros to end to get desired length\n");
	fprintf(stdout, "    [-smooth n]    = optional, does n-point smoothing\n");
	fprintf(stdout, "    [-cut ms ms]   = optional, crops the samples to within a time range\n");
	fprintf(stdout, "    [-copy ms ms]  = optional, copies samples within time range to buffer\n");
	fprintf(stdout, "    [-paste ms]    = optional, pastes samples from copy buffer\n");
	fprintf(stdout, "    [-pasteover ms]= optional, pastes samples from copy buffer (overwriting)\n");
	fprintf(stdout, "    [-resample u d]= optional, resamples signal by u/d (must be integers)\n");
	fprintf(stdout, "                     (eg. -resample 1 2 will halve the sampling rate)\n");
	fprintf(stdout, "    [-reverse]     = optional, reverses the time order of the samples\n");
	fprintf(stdout, "    [-dc]          = optional, subtracts the DC component\n");
	fprintf(stdout, "    [-rectify]     = optional, rectifies the samples\n");
	fprintf(stdout, "    [-sum]         = optional, sums all the input, producing one output entry\n");
	fprintf(stdout, "                     Note that inputs must be same length (can use -cut)\n");
	fprintf(stdout, "    [-avg]         = optional, like -sum, but normalizes by number of inputs\n\n");
	fprintf(stdout, "    <infilename>   = required, describes input file[s].  Many file formats\n");
	fprintf(stdout, "                     are supported.  If %%d is part of filename, then multiple\n");
	fprintf(stdout, "                     input files are assumed\n");
	fprintf(stdout, "    [-entry #]     = optional, start with input file entry %%d (default == 1)\n");
	fprintf(stdout, "    <outfilename>  = required, describes output file[s].  Many file formats\n");
	fprintf(stdout, "                     are supported.  If %%d is part of filename, then multiple\n");
	fprintf(stdout, "                     output files are assumed\n");
	fprintf(stdout, "    [-entry #]     = optional, start with output file entry %%d (default == 1)\n");
	fprintf(stdout, "    [-noout]       = optional, instead of outfilename; create no output file\n\n");
	fprintf(stdout, "  For example, to extract pcm entries from a pcm_seq2 file:\n");
	fprintf(stdout, "    pcmx zy_hv12n*.pcm_seq2 hv12n_%%03d.pcm\n");
	return;
}


int main(int argc, char *argv[])
{
	PCMFILE *infp = NULL, *outfp = NULL;
	int status = SUCCESS;
	char outfname[132], *endptr;
	short *samples, *samples_copy, *samples_strt;
	int nsamples, samples_cnt, i, inentry, outentry, arg, output_is_sequential = 0, scaledB, samplerate, oneoutput;
	long ltime;
	char **filelist = NULL;
	int *entrylist = NULL;
	int filecnt = 0;
	int forcesr = 0;
	struct stat statbuf;

	/* Arguments */
	int domalloc, nooutput;
	char *infname, *rootoutname = NULL;
	procchain_t *nextproc;

	/*
	** Validate arguments
	*/
	if (argc == 1)
	{
		usage();
		return -1;
	}

	/*
	** The user can specify if he wants us to use malloc or mmap in reading files
	*/
	scaledB = -1;
	domalloc = 0;
	nooutput = 0;
	outentry = 1;
	oneoutput = 0;
	proclist = NULL;
	filelist = (char **)malloc(MAXENTRIES * sizeof(char *));
	entrylist = (int *)malloc(MAXENTRIES * sizeof(int));
	filecnt = 0;
	forcesr = 0;
	arg = 1;
	while (arg < argc)
	{
		if (*argv[arg] != '-')
		{
			filecnt++;
			if (filecnt == MAXENTRIES)
			{
				fprintf(stderr, "ERROR: too many inputs!\n");
				return -1;
			}
			filelist[filecnt] = strdup(argv[arg]);
			entrylist[filecnt] = 0;
		}
		else if (!strcmp(argv[arg], "-noout"))
			nooutput = 1;
		else if (!strcmp(argv[arg], "-malloc"))
			domalloc = 1;
		else if (!strcmp(argv[arg], "-resample") && (arg + 2 < argc))
		{
			nextproc = add_proc(PROCESS_RESAMPLE);
			nextproc->intarg1 = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -resample\n");
				return ERROR;
			}
			nextproc->intarg2 = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -resample\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-paste") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_PASTE);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -paste\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-pasteover") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_PASTEOVER);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -pasteover\n");
				return ERROR;
			}
			nextproc->fltarg2 = -1.0;
			if ((arg + 1 < argc) && (isdigit(argv[arg + 1][0])))
			{
				nextproc->fltarg2 = strtod(argv[++arg], &endptr);
				if (*endptr != '\0')
				{
					fprintf(stderr, "Incorrect usage: -pasteover\n");
					return ERROR;
				}
			}
		}
		else if (!strcmp(argv[arg], "-copy") && (arg + 2 < argc))
		{
			nextproc = add_proc(PROCESS_COPY);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -copy\n");
				return ERROR;
			}
			nextproc->fltarg2 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -copy\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-cut") && (arg + 2 < argc))
		{
			nextproc = add_proc(PROCESS_CUT);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -cut\n");
				return ERROR;
			}
			nextproc->fltarg2 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -cut\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-sexconv"))
		{
			nextproc = add_proc(PROCESS_SEXCONV);
		}
		else if (!strcmp(argv[arg], "-invert"))
		{
			nextproc = add_proc(PROCESS_INVERT);
		}
		else if (!strcmp(argv[arg], "-reverse"))
		{
			nextproc = add_proc(PROCESS_REVERSE);
		}
		else if (!strcmp(argv[arg], "-dc"))
		{
			nextproc = add_proc(PROCESS_STRIPDC);
		}
		else if (!strcmp(argv[arg], "-rectify"))
		{
			nextproc = add_proc(PROCESS_RECTIFY);
		}
		else if (!strcmp(argv[arg], "-sum"))
		{
			nextproc = add_proc(PROCESS_SUM);
			oneoutput = 1;
		}
		else if (!strcmp(argv[arg], "-avg"))
		{
			nextproc = add_proc(PROCESS_AVG);
			oneoutput = 1;
		}
		else if (!strcmp(argv[arg], "-prepend") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_PREPEND);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -prepend\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-forcesr") && (arg + 1 < argc))
		{
			forcesr = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -forcesr\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-smooth") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_SMOOTH);
			nextproc->intarg1 = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -smooth\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-pad") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_PAD);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -pad\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-scale") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_SCALE);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -scale\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-rmsscale") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_RMSSCALE);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -rmsscale\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-linramp") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_LINRAMP);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -linramp\n");
				return ERROR;
			}
			nextproc->fltarg2 = nextproc->fltarg1;
			if ((arg + 1 < argc) && (isdigit(argv[arg + 1][0])))
			{
				nextproc->fltarg2 = strtod(argv[++arg], &endptr);
				if (*endptr != '\0')
				{
					fprintf(stderr, "Incorrect usage: -linramp\n");
					return ERROR;
				}
			}
		}
		else if (!strcmp(argv[arg], "-expramp") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_EXPRAMP);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -expramp\n");
				return ERROR;
			}
			nextproc->fltarg2 = nextproc->fltarg1;
			if ((arg + 1 < argc) && (isdigit(argv[arg + 1][0])))
			{
				nextproc->fltarg2 = strtod(argv[++arg], &endptr);
				if (*endptr != '\0')
				{
					fprintf(stderr, "Incorrect usage: -expramp\n");
					return ERROR;
				}
			}
		}
		else if (!strcmp(argv[arg], "-agcscale") && (arg + 2 < argc))
		{
			nextproc = add_proc(PROCESS_AGC);
			nextproc->fltarg1 = strtod(argv[++arg], &endptr);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -agcscale\n");
				return ERROR;
			}
			nextproc->intarg1 = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -agcscale\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-lopass") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_LOPASS);
			nextproc->intarg1 = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -lopass\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-hipass") && (arg + 1 < argc))
		{
			nextproc = add_proc(PROCESS_HIPASS);
			nextproc->intarg1 = strtol(argv[++arg], &endptr, 10);
			if (*endptr != '\0')
			{
				fprintf(stderr, "Incorrect usage: -hipass\n");
				return ERROR;
			}
		}
		else if (!strcmp(argv[arg], "-entry") && (arg + 1 < argc))
		{
			int entry, fentry, lentry;

			arg++;
			fentry = strtol(argv[arg], &endptr, 10);
			if (((*endptr != '\0') && (*endptr != ':')) || (fentry < 1))
			{
				fprintf(stderr, "Incorrect usage: -entry\n");
				return ERROR;
			}
			if (*endptr == ':')
			{
				lentry = strtol(endptr + 1, &endptr, 10);
				if ((*endptr != '\0') || (lentry < fentry))
				{
					fprintf(stderr, "Incorrect usage: -entry\n");
					return ERROR;
				}
			}
			else lentry = fentry;
			entrylist[filecnt] = fentry;
			for (entry = fentry+1; entry <= lentry; entry++)
			{
				filecnt++;
				if (filecnt == MAXENTRIES)
				{
					fprintf(stderr, "ERROR: too many inputs!\n");
					return -1;
				}
				filelist[filecnt] = strdup(filelist[filecnt - 1]);
				entrylist[filecnt] = entry;
			}
		}
		else if ((!strcmp(argv[arg], "-version")) || (!strncmp(argv[arg], "-v", 2)))
		{
		  //			fprintf(stderr, "pcmx build %d version %d\nWritten by Amish S. Dave.\n", PCMX_BUILD, PCMX_VERSION);
			return 0;
		}
		else if ((!strcmp(argv[arg], "-help")) || (!strncmp(argv[arg], "-h", 2)))
		{
			usage();
			return 0;
		}
		else
		{
			fprintf(stderr, "Incorrect usage\n");
			return ERROR;
		}
		arg++;
	}

	/*
	** Determine whether we are to write 1 or many output output files
	** If the rootoutname doesn't allow for entry number substition, then we
	** write multiple entries to the single output file.
	** Otherwise, we create multiple output files, each with a single entry.
	*/
	if (nooutput == 0)
	{
		rootoutname = filelist[filecnt];
		if (strchr(rootoutname, '%'))
		{
			output_is_sequential = 0;
		}
		else
		{
			output_is_sequential = (oneoutput == 0) ? 1 : 0;
			sprintf(outfname, rootoutname);
			if (stat(outfname, &statbuf) != -1)
			{
				fprintf(stderr, "ERROR: output file '%s' already exists! ABORTING.\n", outfname);
				return ERROR;
			}
			fprintf(stderr, "fucking outfname='%s'\n", outfname);
			if ((outfp = pcm_open(outfname, "w")) == NULL)
			{
			  perror("bullshit");
				fprintf(stderr, "ERROR: can't create output file '%s'\n", outfname);
				return ERROR;
			}
		}
		outentry = (entrylist[filecnt] == 0) ? (1) : (entrylist[filecnt]);
	}
	else
	{
		filecnt++;
	}

	/*
	** Loop through each input file
	*/
	infp = NULL;
	infname = NULL;
	for (i = 1; i < filecnt; i++)
	{
		if ((infp != NULL) && (infname != NULL) && (strcmp(infname, filelist[i]) != 0))
		{
			pcm_close(infp);
			infp = NULL;
			infname = NULL;
		}
		infname = filelist[i];
		if (infp == NULL)
			infp = pcm_open(infname, "r");
		if (infp != NULL)
		{
			pcm_ctl(infp, (domalloc == 1) ? PCMIOMALLOC : PCMIOMMAP, NULL);
			inentry = (entrylist[i] == 0) ? (1) : (entrylist[i]);
			for (;;)
			{
				if (pcm_seek(infp, inentry) != -1)
				{
					if (pcm_read(infp, &samples, &nsamples) != -1)
					{
						ltime = -1;
						if (pcm_ctl(infp, PCMIOGETTIME, &ltime) != -1)
							fprintf(stderr, "DATETIME: %s\n", (char *)ctime(&ltime));
						if (pcm_ctl(infp, PCMIOGETSR, &samplerate) != -1)
							fprintf(stderr, "SAMPLERATE: %d\n", samplerate);
						if (forcesr != 0)
							samplerate = forcesr;

						/*
						** Open the output file, if we need to; otherwise, position ourselves at
						** the next entry.
						*/
						if (nooutput == 0)
						{
							if (output_is_sequential == 0)
							{
								if (oneoutput == 0)
								{
									sprintf(outfname, rootoutname, outentry);
									if (stat(outfname, &statbuf) != -1)
									{
										fprintf(stderr, "ERROR: output file '%s' already exists! ABORTING.\n", outfname);
										return ERROR;
									}
									if ((outfp = pcm_open(outfname, "w")) == NULL)
									{
										fprintf(stderr, "ERROR: can't create output file '%s'\n", outfname);
										pcm_close(infp);
										return ERROR;
									}
									pcm_seek(outfp, 1);
								}
							}
							else
							{
								if (pcm_seek(outfp, outentry) == -1)
									fprintf(stderr, "WARNING: seek failed on entry %d of output file - APPENDING\n", outentry);
							}
							if (pcm_ctl(outfp, PCMIOSETSR, &samplerate) == -1)
								fprintf(stderr, "WARNING: couldn't set samplerate to %d\n", samplerate);
							if (ltime != -1)
								if (pcm_ctl(outfp, PCMIOSETTIME, &ltime) == -1)
									fprintf(stderr, "WARNING: couldn't set timestamp to %s\n", (char *)ctime(&ltime));
						}

						/*
						** Process the signal, if so requested
						*/
						samples_copy = NULL;
						samples_strt = samples;
						samples_cnt = nsamples;
						if (proclist != NULL)
						{
							samples_copy = (short *)malloc(nsamples * sizeof(short));
							memcpy(samples_copy, samples, nsamples * sizeof(short));
							samples_strt = samples_copy;
							samples_cnt = nsamples;
							for (nextproc = proclist; nextproc != NULL; nextproc = nextproc->next)
							{
								if (nextproc->type == PROCESS_SCALE)
									scale_samples(samples_strt, samples_cnt, nextproc->fltarg1);
								if (nextproc->type == PROCESS_RMSSCALE)
									rmsscale_samples(samples_strt, samples_cnt, nextproc->fltarg1);
								if (nextproc->type == PROCESS_LINRAMP)
									linramp_samples(samples_strt, samples_cnt, samplerate, nextproc->fltarg1, nextproc->fltarg2);
								if (nextproc->type == PROCESS_EXPRAMP)
									expramp_samples(samples_strt, samples_cnt, samplerate, nextproc->fltarg1, nextproc->fltarg2);
								if (nextproc->type == PROCESS_LOPASS)
									lpfilt_samples(samples_strt, samples_cnt, nextproc->intarg1);
								if (nextproc->type == PROCESS_HIPASS)
									hpfilt_samples(samples_strt, samples_cnt, nextproc->intarg1);
								if (nextproc->type == PROCESS_AGC)
									agcscale_samples(samples_strt, samples_cnt, nextproc->fltarg1, nextproc->intarg1 * samplerate);
								if (nextproc->type == PROCESS_CUT)
									cut_samples(&samples_strt, &samples_cnt, samplerate, nextproc->fltarg1, nextproc->fltarg2);
								if (nextproc->type == PROCESS_COPY)
									copy_samples(&samples_strt, &samples_cnt, samplerate, nextproc->fltarg1, nextproc->fltarg2);
								if (nextproc->type == PROCESS_INVERT)
									invert_samples(samples_strt, samples_cnt);
								if (nextproc->type == PROCESS_REVERSE)
									reverse_samples(samples_strt, samples_cnt);
								if (nextproc->type == PROCESS_SEXCONV)
									sexconv_samples(samples_strt, samples_cnt);
								if (nextproc->type == PROCESS_STRIPDC)
									stripDC_samples(samples_strt, samples_cnt);
								if (nextproc->type == PROCESS_RECTIFY)
									rectify_samples(samples_strt, samples_cnt);
								if (nextproc->type == PROCESS_SUM)
									sum_samples(samples_strt, samples_cnt, (i == (filecnt - 1)));
								if (nextproc->type == PROCESS_AVG)
									avg_samples(samples_strt, samples_cnt, (i == (filecnt - 1)));
								if (nextproc->type == PROCESS_SMOOTH)
									smooth_samples(samples_strt, samples_cnt, nextproc->intarg1);
								if (nextproc->type == PROCESS_PASTE)
									paste_samples(&samples_strt, &samples_cnt, samplerate, nextproc->fltarg1);
								if (nextproc->type == PROCESS_PASTEOVER)
									pasteover_samples(&samples_strt, &samples_cnt, samplerate, nextproc->fltarg1, nextproc->fltarg2);
								if (nextproc->type == PROCESS_PREPEND)
									prepend_samples(&samples_strt, &samples_cnt, samplerate, nextproc->fltarg1);
								if (nextproc->type == PROCESS_RESAMPLE)
									resample_samples(&samples_strt, &samples_cnt, &samplerate, nextproc->intarg1, nextproc->intarg2);
								if (nextproc->type == PROCESS_PAD)
									pad_samples(&samples_strt, &samples_cnt, samplerate, nextproc->fltarg1);

								/*
								** Do the following in case samples_strt has changed - i.e. if one of the functions
								** freed the samples it was given and then reallocated (at possibly different size)
								** a new samples_strt.  samples_copy is what we'll free() below.
								*/
								samples_copy = samples_strt;
							}
						}

						/*
						** Do the actual writing
						*/
						if (nooutput == 0)
						{
							if ((oneoutput == 0) || (i == filecnt - 1))
							{
								if (pcm_write(outfp, samples_strt, samples_cnt) == -1)
								{
									fprintf(stderr, "ERROR: can't write to file '%s'\n", outfname);
									pcm_close(outfp);
									pcm_close(infp);
									return ERROR;
								}
								fprintf(stderr, "INFO: wrote %d samples to entry %d, file '%s'\n", samples_cnt, outentry, outfname);
								outentry++;
							}
						}

						/*
						** Free malloced() memory
						*/
						if (samples_copy != NULL)
							free(samples_copy);
						samples_copy = NULL;
						samples_cnt = 0;
						samples_strt = NULL;

						/*
						** Close the output file if we need to
						*/
						if (nooutput == 0)
						{
							if ((output_is_sequential == 0) && (oneoutput == 0))
							{
								pcm_close(outfp);
								outfp = NULL;
							}
						}
					}
					else 
					{
						fprintf(stderr, "INFO: can't read entry %d of '%s' - MOVING TO NEXT FILE\n", inentry, infname);
						break;
					}
				}
				else
				{
					fprintf(stderr, "INFO: can't seek entry %d of '%s' - MOVING TO NEXT FILE\n", inentry, infname);
					break;
				}
				if (entrylist[i] != 0)
					break;
				inentry++;
			}
		}
		else
		{
			fprintf(stderr, "ERROR: couldn't open input file '%s'\n", infname);
			if (nooutput == 0)
			{
				if (((output_is_sequential == 1) || (oneoutput == 1)) && (outfp != NULL))
					pcm_close(outfp);
			}
			return ERROR;
		}
	}

	if (infp)
	{
		pcm_close(infp);
		infp = NULL;
	}

	/*
	** Close the output file if it was a sequential (ie: pcm_seq2) file
	*/
	if (nooutput == 0)
	{
		if (output_is_sequential == 1)
			pcm_close(outfp);
	}

	/*
	** If there's anything in the copy/paste buffer, free it
	*/
	if (bufsamples != NULL) {free(bufsamples); bufsamples = NULL; bufnsamples = 0;}

	return status;
}

static void invert_samples(short *samples, int nsamples)
{
	int i;

	for (i=0; i < nsamples; i++)
		samples[i] *= -1;
}

static void smooth_samples(short *samples, int nsamples, int nsmooth)
{
	int i, j, num;
	short *sstore;
	float tot;

	sstore = (short *)calloc(nsamples, sizeof(short));

	for (i=0; i < nsamples; i++)
	{
		tot = 0.0;
		num = 0;
		for (j = -1 * (nsmooth/2); j < (nsmooth/2); j++)
		{
			if ((i + j >= 0) && (i + j < nsamples))
			{
				tot += samples[i + j];
				num++;
			}
		}
		sstore[i] = tot/num;
	}
	for (i=0; i < nsamples; i++)
		samples[i] = sstore[i];

	free(sstore);
}

static void pad_samples(short **samples_p, int *nsamples_p, int samplerate, double desiredms)
{
	short *samples;
	int nsamples = *nsamples_p, sampduration;

	sampduration = desiredms * (samplerate / 1000.0);
	samples = (short *)calloc(sampduration, sizeof(short));
	memcpy(samples, *samples_p, MIN(nsamples, sampduration) * sizeof(short));
	free(*samples_p);
	*samples_p = samples;
	*nsamples_p = sampduration;
}

static void resample_samples(short **samples_p, int *nsamples_p, int *samplerate_p, int up, int down)
{
	int n, new_nsamples, old_nsamples;
	short *new_samples, *old_samples;
	int samplerate;

	samplerate = *samplerate_p;
	old_nsamples = *nsamples_p;
	old_samples = *samples_p;
	new_nsamples = old_nsamples * up / down;
	new_nsamples += ((up / down) * 200);
	new_samples = (short *)calloc(new_nsamples + 330, sizeof(short));
	n = rateconv(samplerate, up, down, old_samples, old_nsamples, new_samples);
	if (n > new_nsamples)
	{
		short *new_new_samples;
		int ediff;

	    fprintf(stderr, "WARNING: calculation of returned nsamples is wrong\n");
	    fprintf(stderr, "n = %d, new_nsamples = %d, old_nsamples = %d\n", n, new_nsamples, old_nsamples);
		ediff = n - new_nsamples;
		new_new_samples = (short *)malloc(n * sizeof(short));
		memcpy(new_new_samples, new_samples + (ediff / 2), new_nsamples * sizeof(short));
		free(new_samples);
		new_samples = new_new_samples;
		n = new_nsamples;
	}
	free(old_samples);
	*samples_p = new_samples;
	*nsamples_p = n;
	*samplerate_p = samplerate * (float)(up) / (float)(down);
}

static void prepend_samples(short **samples_p, int *nsamples_p, int samplerate, double prependms)
{
	short *samples;
	int nsamples = *nsamples_p, sampduration, prependsamples;

	prependsamples = (prependms * (samplerate / 1000.0));
	sampduration = *nsamples_p + prependsamples;
	samples = (short *)calloc(sampduration, sizeof(short));
	memcpy(samples + prependsamples, *samples_p, MIN(nsamples, sampduration) * sizeof(short));
	free(*samples_p);
	*samples_p = samples;
	*nsamples_p = sampduration;
}

static void cut_samples(short **samples_p, int *nsamples_p, int samplerate, double start, double stop)
{
	short *samples = *samples_p;
	int nsamples = *nsamples_p, sampstart, sampstop;

	sampstart = start * (samplerate / 1000.0);
	if (stop > start)
		sampstop = stop * (samplerate / 1000.0);
	else
		sampstop = nsamples;
	if (sampstop > sampstart)
	{
		samples += sampstart;
		nsamples = sampstop - sampstart;
	}
	else fprintf(stderr, "ERROR: improper cut parameters.  Skipping cut operation.\n");
	*samples_p = samples;
	*nsamples_p = nsamples;
}

static void copy_samples(short **samples_p, int *nsamples_p, int samplerate, double start, double stop)
{
	short *samples = *samples_p;
	int nsamples = *nsamples_p, sampstart, sampstop;

	sampstart = start * (samplerate / 1000.0);
	if (stop > start)
		sampstop = stop * (samplerate / 1000.0);
	else
		sampstop = nsamples;
	if (sampstop <= sampstart)
	{
		fprintf(stderr, "ERROR: improper copy parameters.  Skipping copy operation.\n");
		return;
	}

	if (bufsamples != NULL) {free(bufsamples); bufsamples = NULL; bufnsamples = 0;}
	bufsamples = (short *)calloc(sampstop - sampstart, sizeof(short));
	bufnsamples = sampstop - sampstart;
	memcpy(bufsamples, samples + sampstart, bufnsamples * sizeof(short));
}

static void paste_samples(short **samples_p, int *nsamples_p, int samplerate, double start)
{
	short *samples = *samples_p;
	short *newsamples = NULL;
	int newnsamples = 0;
	int nsamples = *nsamples_p, sampstart;

	sampstart = start * (samplerate / 1000.0);
	if (sampstart > nsamples)
	{
		fprintf(stderr, "ERROR: improper paste parameter.  Skipping copy operation.\n");
		return;
	}
	if (bufsamples == NULL)
	{
		fprintf(stderr, "Nothing to paste! Aborting...\n");
		exit(-1);
	}

	newnsamples = nsamples + bufnsamples;
	newsamples = (short *)calloc(newnsamples, sizeof(short));
	memcpy(newsamples, samples, sampstart * sizeof(short));
	memcpy(newsamples + sampstart, bufsamples, bufnsamples * sizeof(short));
	memcpy(newsamples + sampstart + bufnsamples, samples + sampstart, (nsamples - sampstart) * sizeof(short));
	free(samples);
	*samples_p = newsamples;
	*nsamples_p = newnsamples;
}

static void pasteover_samples(short **samples_p, int *nsamples_p, int samplerate, double start, double stop)
{
	short *samples = *samples_p;
	short *newsamples = NULL;
	int newnsamples = 0;
	int nsamples = *nsamples_p, sampstart, sampstop = 0;
	int towrite;

	sampstart = start * (samplerate / 1000.0);
	if (stop > 0.0)
	{
		sampstop = stop * (samplerate / 1000.0);
		if (sampstart >= sampstop)
		{
			fprintf(stderr, "ERROR: improper paste parameters.  Aborting.\n");
			exit(-1);
		}
	}
	if (sampstart > nsamples)
	{
		fprintf(stderr, "ERROR: improper paste parameter.  Skipping copy operation.\n");
		return;
	}
	if (bufsamples == NULL)
	{
		fprintf(stderr, "Nothing to paste! Aborting...\n");
		exit(-1);
	}

	newnsamples = nsamples;
	newsamples = (short *)calloc(newnsamples, sizeof(short));
	memcpy(newsamples, samples, nsamples * sizeof(short));

	towrite = nsamples - sampstart;
	if (sampstop != 0)
		towrite = MIN(sampstop - sampstart, towrite);
	else
		towrite = MIN(bufnsamples, towrite);

	{
		short *src = bufsamples;
		short *dest = newsamples + sampstart;
		int i, j;

		for (i=0, j=0; i < towrite; i++, j++)
			dest[i] = src[j % bufnsamples];
	}

	free(samples);
	*samples_p = newsamples;
	*nsamples_p = newnsamples;
}

static void reverse_samples(short *samples, int nsamples)
{
	short *s1, *s2, stmp;

	s1 = samples;
	s2 = &(samples[nsamples - 1]);
	while (s1 < s2)
	{
		stmp = *s1;
		*s1 = *s2;
		*s2 = stmp;
		s1++;
		s2--;
	}
}

static void rectify_samples(short *samples, int nsamples)
{
	int i;

	for (i=0; i < nsamples; i++)
		samples[i] = ABS(samples[i]);
}

static void sum_samples(short *samples, int nsamples, int lasttime)
{
	static int firsttime = 1, summed_nsamples = 0;
	static long *summed_samples = NULL;
	int i, clipped;

	if (firsttime == 1)
	{
		firsttime = 0;
		summed_samples = (long *)calloc(nsamples, sizeof(long));
		summed_nsamples = nsamples;
	}
	else
	{
		if (nsamples != summed_nsamples)
		{
			fprintf(stderr, "ERROR: When summing, all inputs must have same duration.\n");
			exit(-1);
		}
		for (i=0; i < nsamples; i++)
			summed_samples[i] += samples[i];
		if (lasttime == 1)
		{
			clipped = 0;
			for (i=0; i < nsamples; i++)
			{
				if ((summed_samples[i] > 32767) || (summed_samples[i] < -32768))
				{
					clipped++;
					summed_samples[i] = MIN(MAX(summed_samples[i], -32768), 32767);
				}
				samples[i] = summed_samples[i];
			}
			if (clipped > 0)
				fprintf(stderr, "WARNING: clipped %d samples.\n", clipped);
			free(summed_samples);
			summed_samples = 0;
			summed_nsamples = 0;
			firsttime = 1;
		}
	}
}

static void avg_samples(short *samples, int nsamples, int lasttime)
{
	static int firsttime = 1, summed_nsamples = 0, num_reps = 0;
	static long *summed_samples = NULL;
	int i, clipped;

	if (firsttime == 1)
	{
		firsttime = 0;
		summed_samples = (long *)calloc(nsamples, sizeof(long));
		summed_nsamples = nsamples;
		num_reps = 0;
	}
	else
	{
		num_reps++;
		if (nsamples != summed_nsamples)
		{
			fprintf(stderr, "ERROR: When summing, all inputs must have same duration.\n");
			exit(-1);
		}
		for (i=0; i < nsamples; i++)
			summed_samples[i] += samples[i];
		if (lasttime == 1)
		{
			clipped = 0;
			for (i=0; i < nsamples; i++)
			{
				summed_samples[i] /= num_reps;
				if ((summed_samples[i] > 32767) || (summed_samples[i] < -32768))
				{
					clipped++;
					summed_samples[i] = MIN(MAX(summed_samples[i], -32768), 32767);
				}
				samples[i] = summed_samples[i];
			}
			if (clipped > 0)
				fprintf(stderr, "WARNING: clipped %d samples.\n", clipped);
			free(summed_samples);
			summed_samples = 0;
			summed_nsamples = 0;
			firsttime = 1;
			num_reps = 0;
		}
	}
}

static void stripDC_samples(short *samples, int nsamples)
{
	int i;
	long long dcoff;

	dcoff = 0LL;
	for (i=0; i < nsamples; i++)
		dcoff += samples[i];
	dcoff /= nsamples;

	for (i=0; i < nsamples; i++)
		samples[i] -= (short)dcoff;
}

static void sexconv_samples(short *samples, int nsamples)
{
	short *ptr, *end;

#define CONVERTSHORT(in) ((((in) & 0x00ff) << 8 ) | (( (in) & 0xff00) >> 8))
	ptr = samples;
	end = ptr + nsamples;
	while (ptr < end)
	{
		*ptr = CONVERTSHORT(*ptr);
		ptr++;
	}
#undef CONVERTSHORT
}

static void scale_samples(short *samples, int nsamples, double scaledB)
{
	int i, dcoff, new_dcoff, minval, maxval, clipped, abs_maxval;
	double old_rms, new_rms, rmsdb, scalefactor, new_maxval;
	long long rms;

	if (scaledB > 96.0)
		scaledB = 96.0;

	for (dcoff = 0, minval = 33000, maxval = -33000, i=0; i < nsamples; i++)
	{
		dcoff += samples[i];
		if (samples[i] > maxval) maxval = samples[i];
		if (samples[i] < minval) minval = samples[i];
	}
	dcoff = (double)(dcoff) / (double)(nsamples);
	for (rms = 0LL, i=0; i < nsamples; i++)
		rms += SQR(samples[i] - dcoff);
	old_rms = sqrt(rms / (double)nsamples);
	rmsdb = 96.32961 + (20.0 * log10(old_rms / 32767.0));

	printf(" initial: DCoff=%d  min=%d  max=%d  RMS=%.2f dB\n", dcoff, minval, maxval, rmsdb);

	abs_maxval = MAX(ABS(minval - dcoff), ABS(maxval - dcoff));
	new_maxval = 32767.0 / pow(10.0, (96.32961 - scaledB) / 20.0);
	scalefactor = new_maxval / abs_maxval;

	clipped = 0;
	for (new_dcoff = 0, minval = 33000, maxval = -33000, i=0; i < nsamples; i++)
	{
		samples[i] = MAX(MIN(scalefactor * (double)(samples[i] - dcoff), 32767.0), -32768.0);
		if ((samples[i] == 32767) || (samples[i] == -32768)) clipped++;
		new_dcoff += samples[i];
		if (samples[i] > maxval) maxval = samples[i];
		if (samples[i] < minval) minval = samples[i];
	}
	new_dcoff = (double)(new_dcoff) / (double)(nsamples);
	for (rms = 0LL, i=0; i < nsamples; i++)
		rms += SQR(samples[i] - new_dcoff);
	new_rms = sqrt(rms / (double)nsamples);
	rmsdb = 96.32961 + (20.0 * log10(new_rms / 32767.0));

	if (clipped > 0) printf(" WARNING: %d %s clipped.\n", clipped, (clipped == 1) ? ("sample") : ("samples"));
	printf(" final  : DCoff=%d  min=%d  max=%d  RMS=%.2f dB\n", new_dcoff, minval, maxval, rmsdb);
}

static void rmsscale_samples(short *samples, int nsamples, double scaledB)
{
	int i, dcoff, new_dcoff, minval, maxval, clipped;
	double old_rms, new_rms, rmsdb, scalefactor;
	long long rms;

	if (scaledB > 96.0)
		scaledB = 96.0;

	for (dcoff = 0, minval = 33000, maxval = -33000, i=0; i < nsamples; i++)
	{
		dcoff += samples[i];
		if (samples[i] > maxval) maxval = samples[i];
		if (samples[i] < minval) minval = samples[i];
	}
	dcoff = (double)(dcoff) / (double)(nsamples);
	for (rms = 0LL, i=0; i < nsamples; i++)
		rms += SQR(samples[i] - dcoff);
	old_rms = sqrt(rms / (double)nsamples);
	rmsdb = 96.32961 + (20.0 * log10(old_rms / 32767.0));

	printf(" initial: DCoff=%d  min=%d  max=%d  RMS=%.2f dB\n", dcoff, minval, maxval, rmsdb);

	new_rms = 32767.0 / pow(10.0, (96.32961 - scaledB) / 20.0);
	scalefactor = new_rms / old_rms;

	clipped = 0;
	for (new_dcoff = 0, minval = 33000, maxval = -33000, i=0; i < nsamples; i++)
	{
		samples[i] = MAX(MIN(scalefactor * (double)(samples[i] - dcoff), 32767.0), -32768.0);
		if ((samples[i] == 32767) || (samples[i] == -32768)) clipped++;
		new_dcoff += samples[i];
		if (samples[i] > maxval) maxval = samples[i];
		if (samples[i] < minval) minval = samples[i];
	}
	new_dcoff = (double)(new_dcoff) / (double)(nsamples);
	for (rms = 0LL, i=0; i < nsamples; i++)
		rms += SQR(samples[i] - new_dcoff);
	new_rms = sqrt(rms / (double)nsamples);
	rmsdb = 96.32961 + (20.0 * log10(new_rms / 32767.0));

	if (clipped > 0) printf(" WARNING: %d %s clipped.\n", clipped, (clipped == 1) ? ("sample") : ("samples"));
	printf(" final  : DCoff=%d  min=%d  max=%d  RMS=%.2f dB\n", new_dcoff, minval, maxval, rmsdb);
}

static void linramp_samples(short *samples, int nsamples, int samplerate, double onramp, double offramp)
{
	int i, range_samples = 0;
	double scale;

	if (onramp > 0.0)
	{
		range_samples = onramp * (samplerate / 1000.0);
		if (range_samples < nsamples)
		{
			for (i=0; i < range_samples; i++)
			{
				scale = (i * 1.0) / (range_samples - 1.0);			/* Scale goes from 0.0 to 1.0 */
				samples[i] *= scale;
			}
		}
	}

	if (offramp > 0.0)
	{
		range_samples = offramp * (samplerate / 1000.0);
		if (range_samples < nsamples)
		{
			for (i=0; i < range_samples; i++)
			{
				scale = (i * 1.0) / (range_samples - 1.0);			/* Scale goes from 0.0 to 1.0 */
				samples[nsamples - 1 - i] *= scale;
			}
		}
	}

}

static void expramp_samples(short *samples, int nsamples, int samplerate, double onramp, double offramp)
{
	int i, range_samples = 0;
	double scale;

	if (onramp > 0.0)
	{
		range_samples = onramp * (samplerate / 1000.0);
		if (range_samples < nsamples)
		{
			for (i=0; i < range_samples; i++)
			{
				/*
				** Want to generate a value that goes from 0 to 1 exponentially
				** from 0 to range_samples
				*/
				scale = (i * 1.0) / (range_samples - 1.0);			/* Scale goes from 0.0 to 1.0 */
				scale = (exp(scale) - 1.0) / (exp(1.0) - 1.0);
				samples[i] *= scale;
			}
		}
	}

	if (offramp > 0.0)
	{
		range_samples = offramp * (samplerate / 1000.0);
		if (range_samples < nsamples)
		{
			for (i=0; i < range_samples; i++)
			{
				/*
				** Want to generate a value that goes from 0 to 1 exponentially
				** from 0 to range_samples
				*/
				scale = (i * 1.0) / (range_samples - 1.0);			/* Scale goes from 0.0 to 1.0 */
				scale = (exp(scale) - 1.0) / (exp(1.0) - 1.0);
				samples[nsamples - 1 - i] *= scale;
			}
		}
	}

}


static void agcscale_samples(short *samples, int nsamples, double scaledB, int winsamples)
{
	int i, j, dcoff;
	double old_rms, rmsdb, scalefactor;
	long long rms;
	short *window;

	if (scaledB > 96.0)
		scaledB = 96.0;

	window = (short *)calloc(winsamples, sizeof(short));
	memcpy(window, samples, winsamples);

	dcoff = 0;
	for (i=0; i < nsamples; i++)
		dcoff += samples[i];
	dcoff = (double)(dcoff) / (double)(nsamples);

	for (rms = 0LL, i=0; i < winsamples; i++)
		rms += SQR(window[i] - dcoff);
	old_rms = sqrt(rms / (double)winsamples);
	rmsdb = 96.32961 + (20.0 * log10(old_rms / 32767.0));
	scalefactor = (32767.0 / pow(10.0, (96.32961 - scaledB) / 20.0)) / old_rms;

	printf(" initial:  RMS=%.2f dB\n", rmsdb);

	for (i=0; i < winsamples / 2; i++)
		samples[i] = MAX(MIN(scalefactor * (double)(samples[i] - dcoff), 32767.0), -32768.0);
	for (; i < nsamples - winsamples / 2; i++)
	{
		j = (winsamples + (i - winsamples / 2)) % winsamples;
		rms -= SQR(window[j] - dcoff);
		window[j] = samples[i + winsamples / 2];
		rms += SQR(window[j] - dcoff);
		old_rms = sqrt(rms / (double)winsamples);
		scalefactor = (32767.0 / pow(10.0, (96.32961 - scaledB) / 20.0)) / old_rms;
		samples[i] = MAX(MIN(scalefactor * (double)(samples[i] - dcoff), 32767.0), -32768.0);
	}
	for (; i < nsamples; i++)
	{
		samples[i] = MAX(MIN(scalefactor * (double)(samples[i] - dcoff), 32767.0), -32768.0);
	}
	rmsdb = 96.32961 + (20.0 * log10(old_rms / 32767.0));
	printf(" final:  RMS=%.2f dB\n", rmsdb);
	free(window);
}

static void lpfilt_samples(short *samples, int nsamples, int freq)
{
	int i;
	double tg, rs, dl[3][3], err, y;

	/*
	** Initialize the filter coefficients
	*/
	lpcoef(freq / 20000.0, &tg);

	/*
	** Initialize the low and hi-pass delay groups
	*/
	dl[1][1] = 0.0;
	dl[1][2] = 0.0;
	dl[2][1] = 0.0;
	dl[2][2] = 0.0;

	/*
	** Now do the actual lo-pass filtering
	** IIR (recursive) filter.  Two cascaded 2nd order sections.
	*/
	for (i=0; i < nsamples; i++)
	{
		rs = samples[i];
		err = rs + c1[1][1] * dl[1][1] + c1[2][1] * dl[2][1];
		y = c1[5][1] * (err + (c1[3][1] * dl[1][1]) + (c1[4][1] * dl[2][1]));
		dl[2][1] = dl[1][1];
		dl[1][1] = err;
		err = y + c1[1][2] * dl[1][2] + c1[2][2] * dl[2][2];
		rs = c1[5][2] * (err + (c1[3][2] * dl[1][2]) + (c1[4][2] * dl[2][2]));
		dl[2][2] = dl[1][2];
		dl[1][2] = err;
		samples[i] = rs;
	}
}

static void hpfilt_samples(short *samples, int nsamples, int freq)
{
	int i,j,k;
	double tg, rs, dh[3][3], err, y;

	/*
	** Initialize the filter coefficients
	*/
	hpcoef(freq / 20000.0, &tg);

	for (j = 1; j < 6; j++) {
	  for (k = 1; k < 3; k++) {
	    fprintf(stdout, "%3.4f ", c2[j][k]);
	  }
	}
	fprintf(stdout, "\n");

	/*
	** Initialize the low and hi-pass delay groups
	*/
	dh[1][1] = 0.0;
	dh[1][2] = 0.0;
	dh[2][1] = 0.0;
	dh[2][2] = 0.0;

	/*
	** Now do the actual hi-pass filtering
	** IIR (recursive) filter.  Two cascaded 2nd order sections.
	*/
	for (i=0; i < nsamples; i++)
	{
		rs = samples[i];
		err = rs + c2[1][1] * dh[1][1] + c2[2][1] * dh[2][1];
		y = c2[5][1] * (err + (c2[3][1] * dh[1][1]) + (c2[4][1] * dh[2][1]));
		dh[2][1] = dh[1][1];
		dh[1][1] = err;
		err = y + c2[1][2] * dh[1][2] + c2[2][2] * dh[2][2];
		rs = c2[5][2] * (err + (c2[3][2] * dh[1][2]) + (c2[4][2] * dh[2][2]));
		dh[2][2] = dh[1][2];
		dh[1][2] = err;
		samples[i] = rs;
	}
}

/*
** Butterworth low-pass filter coefficients routine.
** Specialized to the 4th order
** fc = cuttoff frequency / sampling frequency
** c1 = coefficients array of dimension 10
** tg = group delay
** nsects = # of 2nd order sections required == 2
*/
static void lpcoef(double fc, double *tg_p)
{
	double tg, omega, osq, rep[5];
	int n, nhalf, i;

	omega = sin(PI * fc) / cos(PI * fc);
	osq = omega * omega;
	n = 4;
	nhalf = n / 2;
	tg = 0.0;

	for (i=1; i <= nhalf; i++)
	{
		rep[i] = omega * cos(PI * (i-0.5) / n);
		tg += rep[i] / osq;
		c1[1][i] = -2.0 * (osq - 1.0) / (1.0 + (2.0 * rep[i]) + osq);
		c1[2][i] = -1.0 * (1.0 - (2.0 * rep[i]) + osq) / (1.0 + 2.0 * rep[i] + osq);
		c1[3][i] = 2.0;
		c1[4][i] = 1.0;
		c1[5][i] = osq / (1.0 + (2.0 * rep[i]) + osq);
	}
}

/*
** Hi-pass filter elliptic coefficients routine.
** Specialized to the 4th order
** fc = cutoff frequency / sampling frequency
** c2 = coefficients array of dimension 10
*/
static void hpcoef(double fc, double *tg_p)
{
	double alp, tg, fcc, c[6][3];
	int istg;

	c[1][1] = 1.395377;
	c[2][1] = -0.7949162;
	c[3][1] = -0.4326008;
	c[4][1] = 1.0;
	c[5][1] = 0.1143855;
	c[1][2] = 1.323464;
	c[2][2] = -0.4867756;
	c[3][2] = 1.10924;
	c[4][2] = 1.0;
	c[5][2] = 0.1143855;

	tg = 0.0;
	fcc = 0.1;
	alp = -1 * cos(PI * (fcc + fc)) / cos(PI * (fcc - fc));
	for (istg = 1; istg <= 2; istg++)
	{
		hphelp(c[1][istg], c[2][istg], 1, alp, &(c2[1][istg]), &(c2[2][istg]));
		hphelp(c[3][istg], c[4][istg], -1, alp, &(c2[3][istg]), &(c2[4][istg]));
		c2[5][istg] = (1.0 + c2[1][istg] - c2[2][istg]) / (1.0 - c2[3][istg] + c2[4][istg]);
		tg = tg + 1.0 + ((-1 * c2[1][istg]) + (2.0 * c2[2][istg])) / (1.0 + c2[1][istg] - c2[2][istg]);
	}
	*tg_p = tg;
}

static void hphelp(double ccc1, double ccc2, int isgn, double alp, double *cc1_p, double *cc2_p)
{
	double c1, c2, cc1, cc2;

	cc1 = *cc1_p;
	cc2 = *cc2_p;
	c1 = isgn*ccc1;
	c2 = isgn*ccc2;
	cc1 = -1 * isgn * ((2.0 * alp) + c1 * (1.0 + (alp * alp)) - (2.0 * alp * c2));
	cc1 = cc1 / (1.0 + (alp * c1) - (alp * alp * c2));
	cc2 = -isgn * ((alp * alp) + (alp * c1) - c2) / (1.0 + (alp * c1) - (alp * alp * c2));
	*cc1_p = cc1;
	*cc2_p = cc2;
	return;
}

