/* file: GenSub.c
 * Collection of aSub subroutines...
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <dbDefs.h>
#include <registryFunction.h>
#include "aSubRecord.h"
#include <epicsExport.h>

#define c_nova_len	60
#define c_novb_len	120

static char str[c_novb_len];

static long outALen(){ return(c_nova_len);}
static long outBLen(){ return(c_novb_len);}

static long command( aSubRecord* pr){
/*------------------------------------------------------------------------------
 * Examine the command pass it to either write/read or write only.
 *----------------------------------------------------------------------------*/
  char* pa=(char*)pr->vala; char* pb=(char*)pr->valb; char* pc=(char*)pr->valc;
  char* p= (char*)pr->a; char* pch=strchr( p,'?');
  sprintf( str,"command: cmnd=%s\nand one more\n",p);
  strcpy( pc,str);
  if(!pch){
    strcpy( pb,p); pr->novb=strlen(p); pa[0]=0; pr->nova=0;
  }
  else{
    strcpy( pa,p); pr->nova=strlen(p); pb[0]=0; pr->novb=0;
  }
  return((int)pch);
}
static long errMsg( aSubRecord* pr){
/*------------------------------------------------------------------------------
 * Reformat the error message from the instrument by adding date/time.
 * NOTE: this routine is not presently used.  The reason is that it is designed
 * to work with edm Messsage Box widget.  I have a modified version which
 * works with this routine.  The standard distribution MBox does not add a new
 * line character after each line of text received and is also limited to 39
 * charaters per line.
 *----------------------------------------------------------------------------*/
  time_t tm=time(0); char tstr[28]; char* ptm; int i,sl;
  char* pa=(char*)pr->vala; char* p= (char*)pr->a;

  if(tm<0) strcpy( tstr,"Failed to get time");
  else{ ptm=ctime( &tm); strncpy( tstr,ptm,28); tstr[27]=0;}
  sl=strlen(tstr);
  if(tstr[sl-1]=='\n') tstr[sl-1]=0;
  strcpy( pa,tstr); strcat( pa," "); strcat( pa,p); /* strcat( pa,"\n"); */
  return(1);
}
static long dbgPrint( aSubRecord* pr){
/*------------------------------------------------------------------------------
 * diagnostic print...
 *----------------------------------------------------------------------------*/
  char* p= (char*)pr->a;

  printf( "dbgPrint: %s\n",p);
  return(1);
}
epicsRegisterFunction( command);
epicsRegisterFunction( dbgPrint);
epicsRegisterFunction( errMsg);
epicsRegisterFunction( outALen);
epicsRegisterFunction( outBLen);
