/* file: tdsGSub.c
 * Collection of genSub subroutines for reprocessing wave form data.
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
#include "genSubRecord.h"
#include <epicsExport.h>

#define PRLEN	256
#define WFLEN	10000

static char str[PRLEN*2];

static long inWFLen(){ return(WFLEN);}
static long inPRLen(){ return(PRLEN);}
static int nbyt,nbit,len,ptof,chn;
static char enc[20],bfmt[20],bord[20],ids[80],ptfm[20],xunt[20],yunt[20];
static double xinc,xzr,ymult,yzr,yof;
static char chs[16];

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
  if((n=sscanf( &p[i],"#%1d",&wd))!=1){
    printf( "tdsWPre: n=%d, failed to get width\n",n);
    return(-1);
  }
  i+=2+wd;
  strncpy( chs,&ids[1],3); chs[3]=0; sscanf( chs,"Ch%d",&chn);
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
static long tdsWFScale( genSubRecord* p){
/*------------------------------------------------------------------------------
 * Scale the waveform data based on the scope setup.
 *----------------------------------------------------------------------------*/
  int nbt=(*(int*)p->a); char* pr=(char*)p->b; float* pwf=(float*)p->vala;
  int on=(*(int*)p->c); float pos=(*(float*)p->d); float vdiv=(*(float*)p->e);
  int dbg=(*(int*)p->f); int nwin=p->nob,nwout=p->nova;
  int i,n,stat=1; float ftmp;
  char* pc; short* pw;
/*  printf( "tdsWFScale: nbt=%d\n",nbt); */
  if(nbt>0&&nbt<=nwin) pr[nbt]=0;
  if((i=tdsWPre( pr,pos,dbg))<0) return(0);
  if(vdiv==0.0) vdiv=1.0;
  pc=(&pr[i]);
  pw=(short*)&pr[i];
  n=nwin<nwout?nwin:nwout;
  n=len<n?len:n;
  n=n<0?0:n;
  for( i=0; i<n; i++,pc++,pw++,pwf++){
    if(nbyt==1) ftmp=(*pc); else ftmp=(*pw);
    if(on) *pwf=((ftmp-yof)*ymult+yzr)/vdiv+pos; else *pwf=0;
  }
  return(stat);
}
epicsRegisterFunction( tdsWFScale);
epicsRegisterFunction( inWFLen);
epicsRegisterFunction( inPRLen);
