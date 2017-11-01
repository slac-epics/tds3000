/* file: tdsGSub.c
 * Collection of aSub subroutines for reprocessing wave form data.
 * Usage of inputs/outputs:
 *  inpa	channel preamble string
 *  inpb	channel waveform data
 *  inpc	channel enable
 *  inpd	channel position
 *  outa	channl1 scaled waveform data
 *  outb	channel reformatted preamble string
 * This code is used to process any of the four channels.
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <dbDefs.h>
#include <registryFunction.h>
#include "aSubRecord.h"
#include <epicsExport.h>
#include <epicsEndian.h>

#define PRLEN	256
#define WFLEN	10000

static char str[PRLEN*2];

static long inWFLen(){ return(WFLEN);}
static long inPRLen(){ return(PRLEN);}
static int nbyt,nbit,len,ptof,chn,endian;
static char enc[20],bfmt[20],bord[20],ids[80],ptfm[20],xunt[20],yunt[20];
static double xinc,xzr,ymult,yzr,yof;
static char chs[16];
static float hscl[]={   1,   2,   4,  10};	/* units of ns */
static int hsx0[]={   225, 200, 150,   0};
static int hsnp[]={    50, 100, 200, 500};

static long tdsInit( aSubRecord* p){
/*------------------------------------------------------------------------------
 * aSub initialization.  Here we set waveform to +10.0 so that it is not
 * visible initialy for all channels.
 *----------------------------------------------------------------------------*/
  int i,npt=p->nova;
  float* pw=(float*)p->vala;
  for( i=0; i<npt; i++,pw++) *pw=10.0;
  return(1);
}
static void getHSParams( double hs,int* x0,int* np){
/*------------------------------------------------------------------------------
 * Returns starting point in x0 and number of points in np where data should be
 * extracted from the input waveform.  Where hs is the horizontal scale in units
 * of a second for which x0 and np are obtained.
 *----------------------------------------------------------------------------*/
  int ix,h=((hs*1e9)+0.5);
  for( ix=0; ix<4&&(hscl[ix]!=h); ix++);
  if(ix>=4) ix=3;
  *x0=hsx0[ix];
  *np=hsnp[ix];
}
static int tdsWPre( char* p,float pos,int dbg){
/*------------------------------------------------------------------------------
 * parses the preamble string in p and decodes various fields of interest.
 *----------------------------------------------------------------------------*/
  int i=0,n,l=strlen(p),wd;
  if(dbg==1){ printf( "+-------------\n%s\n",p);
    printf( "strlen(p)=%d\n-------------\n",l);}
  n=sscanf( p,"%d;%d;%19[^;];%19[^;];%19[^;];%d;%79[^;];%19[^;];"
		"%lg;%d;%lg;%19[^;];%lg;%lg;%lg;%19[^;];%n",
		&nbyt,&nbit,enc,bfmt,bord,&len,ids,ptfm,&xinc,&ptof,&xzr,xunt,
		&ymult,&yzr,&yof,yunt,&i);
  if(n!=16){
    if(n!=5){
      printf( "tdsWPre: n=%d, failed to unpack preamble\n",n); return(-1);}
    return(0);
  }
  // i is the # bytes read from p already. wd reads # of bytes of data
  if((n=sscanf( &p[i],"#%1d",&wd))!=1){
    printf( "tdsWPre: n=%d, failed to get width\n",n);
    return(-1);
  }
  // i is bytes already read + 1 for the '#' + data length
  i+=2+wd;
  strncpy( chs,&ids[1],3); chs[3]=0; sscanf( chs,"Ch%d",&chn);
  if(bord[0] == 'L') endian = EPICS_ENDIAN_LITTLE; else endian = EPICS_ENDIAN_BIG;
  if(dbg==2){
    sprintf( str,"nbytes=%d,  nbits=%d\n",nbyt,nbit); l=strlen(str);
    sprintf( &str[l],"encode=%s, bfmt=%s\n",enc,bfmt); l=strlen(str);
    sprintf( &str[l],"bytOrd=%s,recLen=%d\n",bord,len); l=strlen(str);
    sprintf( &str[l],"%s\n",ids); l=strlen(str);
    sprintf( &str[l],"ptFmt=%s,  ptOfst=%d\n",ptfm,ptof); l=strlen(str);
    sprintf( &str[l],"xInc=%f, xZero=%f %s\n",xinc,xzr,xunt); l=strlen(str);
    sprintf( &str[l],"yMult=%f,  yZero=%f\n",ymult,yzr); l=strlen(str);
    sprintf( &str[l],"yOfst=%f %s, pos=%f\n",yof,yunt,pos); l=strlen(str);
    printf( "%s",str);
    printf( "strlength=%d\n",l);
  }
  else if( dbg==3){}
  return(i);
}
static long tdsWFScale( aSubRecord* p){
/*------------------------------------------------------------------------------
 * Scale the waveform data based on the scope setup.
 * Inputs:
 *  a		number of points in waveform
 *  b		waveform data
 *  c		channel enable
 *  d		trace position
 *  e		volts per division
 *  f		debug flag
 *  g		horizontal scale
 * Outputs:
 *  vala	normalized scaled waveform data
 *  valb	number of data points in vala.
 *----------------------------------------------------------------------------*/
  int nbt=(*(int*)p->a);	char* pr=(char*)p->b;  int on=(*(int*)p->c);
  float pos=(*(float*)p->d);	float vdiv=(*(float*)p->e);
  int dbg=(*(int*)p->f);	double hs=(*(double*)p->g);
  int nwin=p->nob,nwout=p->nova;

  float* pwf=(float*)p->vala;	long* pnp=(long*)p->valb;

  int i,j,x0,np,n,stat=1; float ftmp; char* pc; short* pw; unsigned char* puc; unsigned short* puw;

/*  printf( "tdsWFScale: %s: nbt=%d\n",p->name,nbt);*/
  if(nbt>0&&nbt<=nwin) pr[nbt]=0;
  if((i=tdsWPre( pr,pos,dbg))<0) return(0);
  if(vdiv==0.0) vdiv=1.0;
  getHSParams( hs,&x0,&np);
  pc=(&pr[i]);
  pw=(short*)&pr[i];
  puc=(unsigned char*)(&pr[i]);
  puw=(unsigned short*)&pr[i];
  n=nwin<nwout?nwin:nwout;
  n=len<n?len:n;
  n=n<0?0:n;
  for( i=j=0; i<n; i++,pc++,pw++,puc++,puw++){
	  if(on)
	  {
		  if(nbyt==1)
		  {
			  if(bfmt[1]=='I') // format string is RI (signed) or RP (unsigned)
				  ftmp=(float)(*pc);
			  else
				  ftmp=(float)(*puc);
		  }
		  else
		  {
			  if(endian == EPICS_BYTE_ORDER)
				  if(bfmt[1]=='I') // format string is RI (signed) or RP (unsigned)
					  ftmp=(float)(*pw);
				  else
					  ftmp=(float)(*puw);
			  else // swap byte order, being careful about signed/unsigned
				  if(bfmt[1]=='I')
					  ftmp=(float)(((short)*(pc++) << 8)| *(++puc)); // might do something cleaner than the ++ here
				  else
					  ftmp=(float)(((unsigned short)*puc << 8)| *(++puc));
		  }

		  if(i>=x0&&j<=np)
		  {
			  *pwf=((ftmp-yof)*ymult+yzr)/vdiv+pos;
			  pwf++; j++;
		  }
	  } else *pwf++=(10.0);
	  *pnp=np;
	  /* if(i<10) printf( "tdsWFScale: on=%d, *pwf=%f\n",on,*pwf);*/
  }
  return(stat);
}
epicsRegisterFunction( tdsInit);
epicsRegisterFunction( tdsWFScale);
epicsRegisterFunction( inWFLen);
epicsRegisterFunction( inPRLen);
