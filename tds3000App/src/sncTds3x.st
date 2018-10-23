/* file: sncTds3x.st
 *----------------------------------------------------------------------------*/
program sncTds3x( "P=TDS:,FN=/afs/slac/g/esb/ilc/restore/tds3000.txt")

%%#include <stdio.h>
%%#include <string.h>

#define SPATH 128
#define SETUP 4000

short dochan[4];
assign dochan to {"{P}SQ_CH1_UPD_WF.PROC","{P}SQ_CH2_UPD_WF.PROC",
                  "{P}SQ_CH3_UPD_WF.PROC","{P}SQ_CH4_UPD_WF.PROC"};
short chanrd[4];
assign chanrd to {"{P}R_CH1_ENABLE","{P}R_CH2_ENABLE",
                  "{P}R_CH3_ENABLE","{P}R_CH4_ENABLE"};
monitor chanrd;
evflag chRdFlg; sync chanrd chRdFlg;
short chanwt[4];
assign chanwt to {"{P}W_CH1_ENABLE","{P}W_CH2_ENABLE",
                  "{P}W_CH3_ENABLE","{P}W_CH4_ENABLE"};
short acqur; assign acqur to "{P}COUNT"; monitor acqur;
evflag acqEv; sync acqur acqEv;
short clrR; assign clrR to "{P}CLR_RATE"; monitor clrR;
evflag clrREv; sync clrR clrREv;
short ops; assign ops to "{P}1PER_SECOND.VAL";
short cts; assign cts to "{P}COUNT_GETS.VAL";
int busy[4]; assign busy to {	"{P}C_CH1_BUSY.VAL","{P}C_CH2_BUSY.VAL",
				"{P}C_CH3_BUSY.VAL","{P}C_CH4_BUSY.VAL"};
char spath[SPATH]; assign spath to "{P}SAV_PATH"; monitor spath;
char setup[SETUP]; assign setup to "{P}R_SETUP"; monitor setup;
evflag setpEv; sync setup setpEv;
short rstor; assign rstor to "{P}RESTORE"; monitor rstor;
evflag rstrEv; sync rstor rstrEv;
char wfRst[SETUP]; assign wfRst to "{P}W_SETUP";

int count=0;
int i;
int j;
int chRdOld[4];
int flg1[4];
int flg2[4];
char* fname;

%%static int rdFile( char* msg,char* fn);
%%static void wtFile( char* msg,char* fn);

ss savRest{
  state init{
    when(1){
      fname=macValueGet( "FN");
%%    strncpy( spath,fname,SPATH-1);
      spath[SPATH-1]=0;
      pvPut( spath);
      efClear( setpEv);
      efClear( rstrEv);
      printf( "sncTds3x:savRest: file path initialized\n");
    } state save
  }
  state save{
    when( efTestAndClear( setpEv)){
      wtFile( "sncTds3x",spath);
    } state restore
    when( delay(0.1)){} state restore
  }
  state restore{
    when( efTestAndClear(rstrEv)){
      if(!rdFile( "sncTds3x",spath)) pvPut( wfRst);
    } state save
    when( delay(0.5)){} state save
  }
}

ss acquire{
/* This sequence is used to select next waveform for processing,
 * it is done by looking for next enabled channel...
 *----------------------------------------------------------------------------*/
  state init{
    when(1){
      for( i=0; i<4; i++){
	flg1[i]=flg2[i]=0;
	busy[i]=0; pvPut(busy[0]);
      }
      printf( "sncTds3x:acquire:init\n");
    } state doIt
  }
  state doIt{
    when( efTestAndClear( acqEv)){
      for( i=0; i<4; i++){
/*	pvGet( chanrd[count]); */
/* printf( "doIt: i=%d,count=%d,enable=%d,busy=%d\n",
  i,count,chanrd[count],busy[count]); */
        if(chanrd[count]) {flg1[count]=1; break;}
	else if(flg2[count]) break;
        if((++count)>=4) count=0;
      }
      pvGet( busy[count]);
/* printf( "sncTds3x:acquire:busy[%d]=%d\n",count,busy[count]); */
      if(!busy[count]){
	if(flg1[count]||flg2[count]){
	  busy[count]=1; pvPut( busy[count]);
	  dochan[count]=1; pvPut( dochan[count]);}
	if(flg1[count]){ flg1[count]=0; flg2[count]=1;}
	else if(flg2[count]) flg2[count]=0;
/*
printf( "acquire: Chan%d,busy=%d,Enable=%d,flg1=%d,flg2=%d\n",
	count+1,busy[count],chanrd[count],flg1[count],flg2[count]);
*/
      }
      if((++count)>=4) count=0;
    } state doIt
    when( delay(0.05)){} state doIt
  }
}
ss chkChan{
/* Check the state of channel enable button by reading the PV and if state
 * changed, force the bo enable record to the new state.
 *----------------------------------------------------------------------------*/
  state chkChanInit{
    when(1){
      for( j=0; j<4; j++){
	pvGet( chanrd[j]);
	if(chanrd[j]){ chanwt[j]=1; pvPut( chanwt[j]);}
	chRdOld[j]=chanrd[j];
      }
    } state chkChanDo
  }
  state chkChanDo{
    when( efTestAndClear( chRdFlg)){
      for( j=0; j<4; j++){
	if(chanrd[j]!=chRdOld[j]){
	  chanwt[j]=chanrd[j];
	  pvPut( chanwt[j]);
	  chRdOld[j]=chanrd[j];
	}
      }
    } state chkChanDo
    when( delay(0.05)){} state chkChanDo
  }
}
ss doRate{
/* This sequence is used to clear both counters used in rate calculation.
 *----------------------------------------------------------------------------*/
  state doClr{
    when( efTestAndClear( clrREv)){
      ops=1; pvPut( ops);
      cts=1; pvPut( cts);
    } state doClr
    when( delay(0.5)){} state doClr
  }
}

%{
static int rdFile( char* msg,char* fn){
/*------------------------------------------------------------------------------
 * read configuration string data from fname file.
 *----------------------------------------------------------------------------*/
  FILE* fd; char str[32];

  strcpy( str,msg); strcat( str,":rdFile");
  fd=fopen( fn,"r");
  if(!fd) printf( "%s:rdFile: failed to open %s\n",msg,fn);
  else{
    printf( "%s:rdFile:opened %s\n",msg,fn);
    fgets( wfRst,SETUP-1,fd); wfRst[SETUP-1]=0;
printf( "sncTds3x:rdFile:strlen(wfRst)=%d\n",strlen(wfRst));
    if(ferror( fd)){ perror( str); return(-1);}
  }
  return(0);
}
static void wtFile( char* msg,char* fn){
/*------------------------------------------------------------------------------
 * saves a string of data in a file.
 *----------------------------------------------------------------------------*/
  FILE* fd; int stat;

  fd=fopen( fn,"w");
  if(!fd) printf( "%s:wtFile: failed to open %s\n",msg,fn);
  else{
    setup[SETUP-1]=0;
    stat=fputs( setup,fd);
    if(stat==EOF) printf( "%s:wtFile:fputs: failed to write\n",msg);
    else printf( "%s:wtFile: wrote to %s\n",msg,fn);
    fclose( fd);
  }
}
}%

