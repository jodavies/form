/*
	First the fixed variables
*/
#define MAXPRENAMESIZE 128
/*
	The following variables are default sizes. They can be changed
	into values read from the setup file
*/
#ifdef WORDSIZE32
#define MAXPOWER 500000000
#define MAXVARIABLES 200000050
#define MAXDOLLARVARIABLES 1000000000L
#define WILDOFFSET 200000100
#define MAXINNAMETREE 2000000000
#define MAXDUMMIES 100000000
#define WORKBUFFER 20000000
#define MAXTER 40000
#define HALFMAX 0x10000
#else
#define MAXPOWER 10000
#define MAXVARIABLES 8050
#define MAXDOLLARVARIABLES 32000
#define WILDOFFSET 6100
#define MAXINNAMETREE 32768
#define MAXDUMMIES 1000
#define WORKBUFFER 1000000
#define MAXTER 10000
#define HALFMAX 0x100
#endif
#define MAXENAME 16

#define MAXPARLEVEL 100
#define MAXNUMBERSIZE 200

#define MAXREPEAT 100
#define NORMSIZE 1000

#define INITNODESIZE 10
#define INITNAMESIZE 100

#define NUMFIXED 128
#define MAXNEST 100
#define MAXMATCH 30
#define MAXIF 20
#define SIZEFACS 640L
#define NUMFACS 50
#define MAXLOOPS 30
#define MAXLABELS 20
#define COMMERCIALSIZE 24
/*
	The next quantities should still be eliminated from the program
	This should be together with changes in setfile!
*/
#define COMPRESSBUFFER 90000
#define FORTRANCONTINUATIONLINES 15
#define MAXLEVELS 2000
#define MAXLHS 400
#define MAXWILDC 100
#define NUMTABLEENTRIES 1000
#define COMPILERBUFFER 20000

#define SMALLBUFFER   10000000L
#define SMALLOVERFLOW 20000000L
#define TERMSSMALL      100000L
#define LARGEBUFFER   50000000L
#define MAXPATCHES 256
#define MAXFPATCHES 256
#define SORTIOSIZE 100000L

#define SSMALLBUFFER 500000L
#define SSMALLOVERFLOW 800000L
#define STERMSSMALL 10000L
#define SLARGEBUFFER 4000000L
#define SMAXPATCHES 64
#define SMAXFPATCHES 64
#define SSORTIOSIZE 32768L

#define SCRATCHSIZE 50000000L

#define MAXFLEVELS 30

#define COMPINC 2
 
#define MAXNUMSIZE 10

#define MAXBRACKETBUFFERSIZE 200000

#define ZIPBUFFERSIZE 32768L
#define SFHSIZE 40

#define SLAVEPATCHSIZE 10
#define SHMWINSIZE     65536L

#define TABLEEXTENSION 6

#define GZIPDEFAULT 0
#define DEFAULTTHREADS 0
#define DEFAULTTHREADBUCKETSIZE 500
#define DEFAULTTHREADLOADBALANCING 1
#define THREADSCRATCHSIZE 100000L
#define THREADSCRATCHOUTSIZE 2500000L

#ifdef WORDSIZE32
#define MAXTABLECOMBUF 100000000000L
#define MAXCOMBUFRHS 1000000000L
#else
#define MAXTABLECOMBUF 1000000L
#define MAXCOMBUFRHS 32500L
#endif

#define NUMSTORECACHES 4
#define SIZESTORECACHE 32768
