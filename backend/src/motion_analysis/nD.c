/* nD.c
   Dennis M. Freeman       July 1998    */

#include "nD.h"
#include <linux/limits.h>
#include <string.h>
#include <errno.h>

void perr(char* msg,char* name){fprintf(stderr,"%s \"%s",msg,name); perror("\"");
  exit(1);}

void err(char* msg,char* name){fprintf(stderr,"%s \"%s\"\n",msg,name); exit(1);}

void datetime(char* dstr,char* tstr){struct tm *tmt; time_t tt;
  tt = time(0);  tmt = localtime(&tt); tmt->tm_mon++;
  sprintf(dstr,"%d-%d%d-%d%d",tmt->tm_year+1900,tmt->tm_mon/10,tmt->tm_mon%10,
    tmt->tm_mday/10,tmt->tm_mday%10);
  sprintf(tstr,"%d%d:%d%d:%d%d",tmt->tm_hour/10,tmt->tm_hour%10,
    tmt->tm_min/10,tmt->tm_min%10,tmt->tm_sec/10,tmt->tm_sec%10);}

char* commandline(int argc,char **argv){int i,n; char* t;
  for(i=n=0; i<argc; i++) n += 1+strlen(argv[i]);
  if((t=malloc(n))==0) perr("Commandline malloc error",argv[0]);
  strcpy(t,argv[0]);
  for(i=1; i<argc; i++) {strcat(t," "); strcat(t,argv[i]);}
  return t;}

nD* open_nD(char* name){nD *id; int i;
  if((id=malloc(sizeof (nD)+strlen(name)+1))==0)
    perr("open_nD malloc error",name);
  id->name = (char*)(id)+sizeof (nD);  strcpy(id->name,name);
  if((id->fd=open(name,O_RDONLY,0)) < 0) perr("open_nD open error",name);
  id->mapsize = lseek(id->fd,0,2);
  if((long)(id->file=mmap(NULL,id->mapsize,PROT_READ,MAP_PRIVATE,id->fd,0))
    ==(-1)) perr("open_nD mmap error",name);
  if(sscanf(id->file,"nD %d %d",&(id->dataoffset),&(id->notesoffset))==2){
    if((id->dataoffset>id->notesoffset)||(id->notesoffset>id->mapsize))
      err("Bad nD file: file too short",name);
    id->data = id->file+id->dataoffset;
    id->notes = id->file+id->notesoffset;
    id->prefix = index(id->file,'\n'); (id->prefix)++;}
  else {
    if(id->file[0]!='P'||id->file[1]!='5'||id->file[2]!='\n')
      err("File is neither nD or pgm",name);
    for(i=3; i<id->mapsize; i++){
      if(id->file[i]!='\n') continue;
      if(id->file[i+1]!='2') continue;
      if(id->file[i+2]!='5') continue;
      if(id->file[i+3]!='5') continue;
      if(id->file[i+4]=='\n') break;}
    id->data = id->file+i+5; id->dataoffset = i+5;
    id->notes = id->file+id->mapsize; id->notesoffset = id->mapsize;
    id->prefix = id->file;}
  id->suffixhead = id->suffixtail = 0;
  id->merged = 0;
  id->width= id->height= id->times= id->planes= 0;
  return(id);}

nD* create_nD(char* name,char* prefix,int datasize,char* notes){
  int i,j,nprefix,nnotes;  nD *id;  char header[100];
  if((id=malloc(sizeof (nD)+strlen(name)+1))==0)
    perr("create_nD malloc error",name);
  id->name = (char*)(id)+sizeof (nD);  strcpy(id->name,name);
  nprefix = strlen(prefix);  nnotes = strlen(notes);
  id->dataoffset = nprefix;
  id->notesoffset = id->dataoffset + datasize;
  sprintf(header,"nD %d %d\n",id->dataoffset,id->notesoffset);
  while(id->dataoffset!=strlen(header)+nprefix){
    id->dataoffset = strlen(header) + nprefix;
    id->notesoffset = id->dataoffset + datasize;
    sprintf(header,"nD %d %d\n",id->dataoffset,id->notesoffset);}
  id->mapsize = id->notesoffset+strlen(notes);
  if((id->fd=open(name,O_CREAT|O_RDWR|O_TRUNC,0666)) < 0)
    perr("create_nD open error",name);
  lseek(id->fd,id->mapsize-1,0);
  if((write(id->fd," ", 1)) < 0) perr("create_nD write error",name);
  if((unsigned long)(id->file=
            mmap(NULL,id->mapsize,PROT_WRITE,MAP_SHARED,id->fd,0))==-1) 
    perr("create_nD mmap error",name);
  i = 0;  j = 0; while(header[j]!=0) id->file[i++] = header[j++];
  id->prefix = &(id->file[i]);
  j = 0; while(prefix[j]!=0) id->file[i++] = prefix[j++];
  id->data = id->file+id->dataoffset;
  id->notes = id->file+id->notesoffset;
  for(i=0; i<=nnotes; i++) id->notes[i] = notes[i];
  id->suffixhead = id->suffixtail = 0;
  id->merged = 0;
  id->width= id->height= id->times= id->planes= 0;
  return(id);}

void close_nD(nD* id){
  suffixlink *t,*link;

  munmap(id->file,id->mapsize);  close(id->fd);
  if(id->suffixhead){
    if((id->fd=open(id->name,O_RDWR)) < 0) perr("nD reopen error",id->name);
    lseek(id->fd,0,2);
    link = id->suffixhead;
    while(link){
      write(id->fd,link->line,link->nline);
      t = link->next;
      free(link);
      link = t;
    }
    close(id->fd);
  }
  if(id->merged) free(id->merged);
  free(id);
}

void putnotes(nD* id,char* line){int n; suffixlink *link;
  if(id->merged) {free(id->merged); id->merged = 0;}
  n = strlen(line);
  if((link=malloc(n+1+sizeof (suffixlink)))==0)
    perr("putnotes malloc error",id->name);
  link->next = 0;
  link->nline = n;
  link->line = (char*)(link)+sizeof (suffixlink);
  strcpy(link->line,line);
  if(id->suffixtail==0) id->suffixhead = link;
    else (id->suffixtail)->next = link;
  id->suffixtail = link;}

char* getnotes(nD* id){int n; suffixlink *link; char *notes;
  if(id->merged) return(id->merged);
  n = (id->data - id->prefix) + (id->mapsize - id->notesoffset);
  link = id->suffixhead; while(link){n += (link->nline); link = link->next;}
  if((id->merged=malloc(n+1))==0) perr("getnotes malloc error",id->name);
  notes = id->merged;
  memcpy(notes,id->prefix,id->data - id->prefix);
  notes += id->data - id->prefix;
  memcpy(notes,id->notes,id->mapsize - id->notesoffset);
  notes += id->mapsize - id->notesoffset;
  link = id->suffixhead;
  while(link){
    strcpy(notes,link->line);
    notes += link->nline;
    link = link->next;}
  notes[0] = 0;
  return id->merged;}

char* getkey(nD* id,char* key,char* ans,int nans){char* t; int i;
  if(nans<strlen(key)+4) err("getkey answer string too short",id->name);
  sprintf(ans,"\n=%s=",key);
  t = strstr(getnotes(id),ans);
  if(t==0) return 0;
  t += 3+strlen(key);
  for(i=0; (i<nans-1)&&(t[i]!='\n')&&(t[i]!='\0'); i++) ans[i] = t[i];
  ans[i] = 0;
  return ans;}

void putkey(nD* id,char* key,char* val){char* t;
  if((t=malloc(strlen(key)+strlen(val)+4))==0)
    perr("putkey malloc error",id->name);
  sprintf(t,"=%s=%s\n",key,val);
  putnotes(id,t);
  free(t);}

float getU8(nD* id, unsigned int i, unsigned int j){
  return((float)((id->d.rowb[j])[i]));}

float getU10(nD* id, unsigned int i, unsigned int j){
  unsigned int p; unsigned short xs;
  p = (i*5)/4;
  xs = (((unsigned short)((id->d.rowb[j])[p+1]))<<8) + (id->d.rowb[j])[p];
  if((i&3)==1) xs = xs>>2;
  if((i&3)==2) xs = xs>>4;
  if((i&3)==3) xs = xs>>6;
  return (float)(xs&0x3FF);}

float getU12(nD* id, unsigned int i, unsigned int j){
  unsigned int p;  unsigned short xs;
  p = (i*3)/2;  
  xs = (((unsigned short)((id->d.rowb[j])[p+1]))<<8) + (id->d.rowb[j])[p];
  if(i&1) return (float)((xs>>4)&0xFFF);
  return (float)(xs&0xFFF);}

float getI16(nD* id, unsigned int i, unsigned int j){
  return((float)((id->d.rows[j])[i]));}

float getF32(nD* id, unsigned int i, unsigned int j){
  return((float)((id->d.rowf[j])[i]));}

float rgetC64(nD* id, unsigned int i, unsigned int j){
  return((float)((id->d.rowf[j])[i+i]));}

complex cgetC64(nD* id, unsigned int i, unsigned int j){complex t;
  t.r = (id->d.rowf[j])[i+i];  t.i = (id->d.rowf[j])[i+i+1];  return t;}

complex cgetr(nD* id, unsigned int i, unsigned int j){complex t;
  t.r = id->get.n2D(id,i,j); t.i = 0; return t;}

float rgetRGB24(nD* id, unsigned int i, unsigned int j){
  return((float)(((id->d.rowb[j])[i+i+i]+(id->d.rowb[j])[i+i+i+1]
    +(id->d.rowb[j])[i+i+i+2])/3.));}

rgb rgbgetRGB24(nD* id, unsigned int i, unsigned int j){rgb t;
  t.r = (float)((id->d.rowb[j])[i+i+i]);
  t.g = (float)((id->d.rowb[j])[i+i+i+1]);
  t.b = (float)((id->d.rowb[j])[i+i+i+2]);
  return t;}

rgb rgbgetr(nD* id, unsigned int i, unsigned int j){rgb t;
  t.r = t.g = t.b = id->get.n2D(id,i,j);  return t;}

float rgetRGB48(nD* id, unsigned int i, unsigned int j){
  return((float)(((id->d.rows[j])[i+i+i]+(id->d.rows[j])[i+i+i+1]
    +(id->d.rows[j])[i+i+i+2])/3.));}

rgb rgbgetRGB48(nD* id, unsigned int i, unsigned int j){rgb t;
  t.r = (float)((id->d.rows[j])[i+i+i]);
  t.g = (float)((id->d.rows[j])[i+i+i+1]);
  t.b = (float)((id->d.rows[j])[i+i+i+2]);
  return t;}

void badput2(nD* id, float value, unsigned int i, unsigned int j){
  err("Attempted to write to read-only file",id->name);}

nD* open2D(char* name){nD* id; int i,rowsize=0;
  id = open_nD(name);
  if(sscanf(id->prefix,"# 2D %5s %d %d",id->pixel,&(id->width),&(id->height))==3){
    id->cget.n2D = cgetr;  id->rgbget.n2D = rgbgetr;
    if(strcmp(id->pixel,"U8")==0)
      {id->get.n2D = getU8; rowsize = id->width;}
    else if(strcmp(id->pixel,"U10")==0)
      {id->get.n2D = getU10; rowsize = ceil(id->width*10./8.);}
    else if(strcmp(id->pixel,"U12")==0)
      {id->get.n2D = getU12; rowsize = ceil(id->width*12./8.);}
    else if(strcmp(id->pixel,"I16")==0)
      {id->get.n2D = getI16; rowsize = id->width+id->width;}
    else if(strcmp(id->pixel,"F32")==0)
      {id->get.n2D = getF32; rowsize = 4*id->width;}
    else if(strcmp(id->pixel,"C64")==0)
      {id->get.n2D = rgetC64; id->cget.n2D = cgetC64; rowsize = 8*id->width;}
    else if(strcmp(id->pixel,"RGB24")==0)
      {id->get.n2D = rgetRGB24; id->rgbget.n2D = rgbgetRGB24; rowsize=3*id->width;}
    else if(strcmp(id->pixel,"RGB48")==0)
      {id->get.n2D = rgetRGB48; id->rgbget.n2D = rgbgetRGB48; rowsize=2*3*id->width;}
    else err("Bad 2D pixel type",id->name);}
  else {
    if(id->prefix[0]!='P'||id->prefix[1]!='5'||id->prefix[2]!='\n')
      err("Bad 2D file",name);
    for(i=id->dataoffset-6; i>1; i--) if(id->prefix[i]=='\n') break;
    if(i==1) {
      err("No size specifier in pgm file",name);
    }
    if(sscanf(id->prefix+i+1,"%d %d",&(id->width),&(id->height))!=2)
      err("Bad size specifier in pgm file",name);
    strcpy(id->pixel,"U8"); id->get.n2D = getU8; rowsize = id->width;}
  if((rowsize*(id->height))!=(id->notesoffset)-(id->dataoffset))
    err("2D file has wrong length",id->name);
  id->put.n2D = badput2;
  if((id->d.rowb=malloc((id->height)*sizeof (char*)))==0)
    perr("open2D malloc error",name);
  id->d.rowb[(id->height)-1] = id->data;
  for(i=(id->height)-2; i>=0; i--) id->d.rowb[i] = id->d.rowb[i+1]+rowsize;
  return(id);}

void putU8(nD* id, float value, unsigned int i, unsigned int j){
  if(value>255) id->d.rowb[j][i]= 255;
  else (id->d.rowb[j])[i] = value;}

void putU10(nD* id, float value, unsigned int i, unsigned int j){
  unsigned int p; unsigned short t, xs;
  t = ((unsigned short)(value)) & 0x03FF;  p = (i*5)/4;
  xs = (((unsigned short)((id->d.rowb[j])[p+1]))<<8) + (id->d.rowb[j])[p];
  if((i&3)==0) xs &= 0xFC00;
  else if((i&3)==1) {t = (t<<2); xs &= 0xF003;}
  else if((i&3)==2) {t = (t<<4); xs &= 0xC00F;}
  else if((i&3)==3) {t = (t<<6); xs &= 0x003F;}
  xs |= t; (id->d.rowb[j])[p] = xs&0xFF; (id->d.rowb[j])[p+1] = xs>>8;}

void putU12(nD* id, float value, unsigned int i, unsigned int j){
  unsigned int p; unsigned short t, xs;
  t = ((unsigned short)(value)) & 0x0FFF;  p = (i*3)/2;
  xs = (((unsigned short)((id->d.rowb[j])[p+1]))<<8) + (id->d.rowb[j])[p];
  if((i&1)==0) {xs &= 0xF000;} else {t = t<<4; xs &= 0x000F;}
  xs |= t; (id->d.rowb[j])[p] = xs&0xFF; (id->d.rowb[j])[p+1] = xs>>8;}

void putI16(nD* id, float value, unsigned int i, unsigned int j){
  (id->d.rows[j])[i] = value;}

void putF32(nD* id, float value, unsigned int i, unsigned int j){
  (id->d.rowf[j])[i] = value;}

void rputC64(nD* id, float value, unsigned int i, unsigned int j){
  (id->d.rowf[j])[i+i] = value;  (id->d.rowf[j])[i+i+1] = 0;}

void cputC64(nD* id, complex value, unsigned int i, unsigned int j){
  (id->d.rowf[j])[i+i] = value.r;  (id->d.rowf[j])[i+i+1] = value.i;}

void cputr(nD* id, complex value, unsigned int i, unsigned int j){
  id->put.n2D(id,value.r,i,j);}

void rputRGB24(nD* id, float value, unsigned int i, unsigned int j){
  (id->d.rowb[j])[i+i+i] = value;
  (id->d.rowb[j])[i+i+i+1] = value;
  (id->d.rowb[j])[i+i+i+2] = value;}

void rgbputRGB24(nD* id, rgb value, unsigned int i, unsigned int j){
  (id->d.rowb[j])[i+i+i] = value.r;
  (id->d.rowb[j])[i+i+i+1] = value.g;
  (id->d.rowb[j])[i+i+i+2] = value.b;}

void rgbputr(nD* id, rgb value, unsigned int i, unsigned int j){
  id->put.n2D(id,(value.r+value.g+value.b)/3.,i,j);}

void rputRGB48(nD* id, float value, unsigned int i, unsigned int j){
  (id->d.rows[j])[i+i+i] = value;
  (id->d.rows[j])[i+i+i+1] = value;
  (id->d.rows[j])[i+i+i+2] = value;}

void rgbputRGB48(nD* id, rgb value, unsigned int i, unsigned int j){
  (id->d.rows[j])[i+i+i] = value.r;
  (id->d.rows[j])[i+i+i+1] = value.g;
  (id->d.rows[j])[i+i+i+2] = value.b;}

void rgbputr48(nD* id, rgb value, unsigned int i, unsigned int j){
  id->put.n2D(id,(value.r+value.g+value.b)/(256*3.),i,j);}

nD* createPGM(char* name,unsigned int width,unsigned int height){
  nD* id; char *prefix; int rowsize=1; char header[100]; int i;
  if((prefix=malloc(10000))==0) perr("createPGM malloc error",name);
  rowsize = width;
  if((id=malloc(sizeof(nD)+strlen(name)+1))==0)
    perr("createPGM malloc error",name);
  id->name = (char*)(id)+sizeof(nD);  strcpy(id->name,name);
  sprintf(header,"P5\n%d %d\n255\n",width,height);
  id->dataoffset = strlen(header);
  id->notesoffset = id->dataoffset + width*height;
  id->mapsize = id->notesoffset;
  if((id->fd=open(name,O_CREAT|O_RDWR|O_TRUNC,0666)) < 0)
    perr("createPGM open error",name);
  lseek(id->fd,id->mapsize-1,0);
  if((write(id->fd," ", 1)) < 0) perr("createPGM write error",name);
  if((long)(id->file=
    mmap(NULL,id->mapsize,PROT_WRITE,MAP_SHARED,id->fd,0))<0)
    perr("createPGM mmap error",name);
  for(i=0; header[i]!=0; i++) id->file[i] = header[i];
  id->prefix = &(id->file[i]);
  id->data = id->file+id->dataoffset;
  id->notes = id->file+id->notesoffset;
  id->suffixhead = id->suffixtail = 0;
  id->merged = 0;
  id->cget.n2D = cgetr; id->cput.n2D = cputr;
  id->rgbget.n2D = rgbgetr; id->rgbput.n2D = rgbputr;
  id->put.n2D = putU8; id->get.n2D = getU8;
  id->width = width;
  id->height = height;
  strcpy(id->pixel,"U8");
  free(prefix);
  if((id->d.rowb=malloc((id->height)*sizeof (char*)))==0)
    perr("open2D malloc error",name);
  id->d.rowb[(id->height)-1] = id->data;
  for(i=(id->height)-2; i>=0; i--) id->d.rowb[i] = id->d.rowb[i+1]+rowsize;
  return(id);}

nD* create2D(char* name,char* pixel,unsigned int width,unsigned int height,
  char* commandline,char* notes){
  nD* id; char *prefix; int rowsize=1; char date[100],time[100]; int i;
  if((prefix=malloc(10000))==0) perr("create2D malloc error",name);
  if(strcmp(pixel,"U8")==0) rowsize = width;
  else if(strcmp(pixel,"U10")==0) rowsize = ceil(width*10./8.);
  else if(strcmp(pixel,"U12")==0) rowsize = ceil(width*12./8.);
  else if(strcmp(pixel,"I16")==0) rowsize = width + width;
  else if(strcmp(pixel,"F32")==0) rowsize = 4*width;
  else if(strcmp(pixel,"C64")==0) rowsize = 8*width;
  else if(strcmp(pixel,"RGB24")==0) rowsize = 3*width;
  else if(strcmp(pixel,"RGB48")==0) rowsize = 3*2*width;
  else err("Bad 2D pixel type",name);
  datetime(date,time);
  sprintf(prefix,"# 2D %s %d %d\n=date=%s\n=time=%s\n=command=%s\n",
    pixel,width,height,date,time,commandline);
  id = create_nD(name,prefix,height*rowsize,notes);
  id->cget.n2D = cgetr; id->cput.n2D = cputr;
  id->rgbget.n2D = rgbgetr; id->rgbput.n2D = rgbputr;
  if(strcmp(pixel,"U8")==0)
    {id->put.n2D = putU8; id->get.n2D = getU8;}
  else if(strcmp(pixel,"U10")==0)
    {id->put.n2D = putU10; id->get.n2D = getU10;
    /* you get a segmentation error if you read a location that wasn't written.*/
    for(i=0; i<height*rowsize; i++) id->data[i] = 0;}
  else if(strcmp(pixel,"U12")==0)
    {id->put.n2D = putU12; id->get.n2D = getU12;
    for(i=0; i<height*rowsize; i++) id->data[i] = 0;}
  else if(strcmp(pixel,"I16")==0)
    {id->put.n2D = putI16; id->get.n2D = getI16;}
  else if(strcmp(pixel,"F32")==0)
    {id->put.n2D = putF32; id->get.n2D = getF32;}
  else if(strcmp(pixel,"C64")==0)
    {id->put.n2D = rputC64; id->cput.n2D = cputC64;
    id->get.n2D = rgetC64; id->cget.n2D = cgetC64;}
  else if(strcmp(pixel,"RGB24")==0)
    {id->put.n2D = rputRGB24; id->rgbput.n2D = rgbputRGB24;
    id->get.n2D = rgetRGB24; id->rgbget.n2D = rgbgetRGB24;}
  else if(strcmp(pixel,"RGB48")==0)
    {id->put.n2D = rputRGB48; id->rgbput.n2D = rgbputRGB48;
    id->get.n2D = rgetRGB48; id->rgbget.n2D = rgbgetRGB48;}
  id->width = width;
  id->height = height;
  strcpy(id->pixel,pixel);
  free(prefix);
  if((id->d.rowb=malloc((id->height)*sizeof (char*)))==0)
    perr("open2D malloc error",name);
  id->d.rowb[(id->height)-1] = id->data;
  for(i=(id->height)-2; i>=0; i--) id->d.rowb[i] = id->d.rowb[i+1]+rowsize;
  return(id);}

void close2D(nD* id){
  if(id) {
    if(id->d.rowb) {
      free(id->d.rowb);
    } 
    close_nD(id);
  }
}

float get3Dx(nD* id, unsigned int i, unsigned int j, unsigned int k){
  return(((id->d.z[k])->get.n2D)(id->d.z[k],i,j));}
complex cget3Dx(nD* id, unsigned int i, unsigned int j, unsigned int k){
  return(((id->d.z[k])->cget.n2D)(id->d.z[k],i,j));}
void badput3(nD* id, float value, unsigned int i, unsigned int j, unsigned int k){
  err("Attempted to write to read-only file",id->name);}

nD* open3D(char* name){nD *id; unsigned char *name2D,*t; int i,j;
  id = open_nD (name);
  if(sscanf(id->prefix,"# 3D %5s %d %d %d",id->pixel,&(id->width),
    &(id->height),&(id->planes))!=4) err("Bad 3D file",id->name);
  id->get.n3D = get3Dx;
  id->cget.n3D = cget3Dx;
  id->put.n3D = badput3;
  if((id->d.z=malloc((id->planes)*sizeof (nD*)))==0)
    perr("open3D malloc error",name);
  if((name2D=malloc(1+(id->notesoffset)-(id->dataoffset)))==0)
    perr("open3D malloc error",name);
  t = id->data;
  for(i=0; i<id->planes; i++){
    for(j=0; t<id->notes&&t[0]!='\n'; j++) {name2D[j] = t[0]; t++;}
    if(j==0) err("open3D error: Empty 2D filename in 3D image",name);
    name2D[j] = 0;  t++;
    id->d.z[i] = open2D(name2D);
  }
  for(; t<id->notes; t++)
    if((t[0]!='\n')&&(t[0]!=' '))
      err("open3D error: Too many 2D files in 3D image",name);
  free(name2D);
  return(id);}

nD* create3Dfrom2D(char* name,char* pixel,unsigned int width,unsigned int height,
  unsigned int planes,char* files,char* commandline,char* notes){
  nD* id; int k; char *prefix; char date[100],time[100];
  if((prefix=malloc(10000))==0) perr("create3Dfrom2D malloc error",name);
  datetime(date,time);
  sprintf(prefix,"# 3D %s %d %d %d\n=date=%s\n=time=%s\n=command=%s\n",
    pixel,width,height,planes,date,time,commandline);
  id = create_nD(name,prefix,strlen(files),notes);
  free(prefix);
  strncpy(id->data,files,strlen(files));
  close_nD(id);
  id = open3D(name);
  for(k=0; k<planes; k++){
    if(id->width!=(id->d.z[k])->width)
      perr("create3Dfrom2D error: Incompatible image widths in",name);
    if(id->height!=(id->d.z[k])->height)
      perr("create3Dfrom2D error: Incompatible image heights in",name);}
  return(id);}

nD* create4Dfrom3D(char *name,char* pixel,unsigned int width,unsigned int height,
                   unsigned int planes,unsigned int times,char *files,
                   char *commandline,char* notes){
  nD *id, *id3d; int t,k; char *prefix; char date[100],time[100];

  if((prefix=malloc(10000))==0) perr("create4Dfrom3D malloc error",name);
  datetime(date,time);
  sprintf(prefix,"# 4D %s %d %d %d %d\n=date=%s\n=time=%s\n=command=%s\n",
          pixel,width,height,planes,times,date,time,commandline);
  id = create_nD(name,prefix,strlen(files),notes);
  free(prefix);
  strncpy(id->data,files,strlen(files));
  close_nD(id);
  id = open4Dasnames(name);
  for(t=0; t<times;t++) {
    id3d = open3D(id->nDnames3D[t]);
    if(id->planes != (id3d->planes))
      perr("create4Dfrom3D error: Incompatible number of planes in",name);
    for(k=0;k<planes;k++) {
      if(id->width!=(id3d->d.z[k]->width))
        perr("create4Dfrom3D error: Incompatible image widths in",name);
      if(id->height!=(id3d->d.z[k]->height))
        perr("create4Dfrom3D error: Incompatible image widths in",name);
    }
    close3D(id3d);
  }
  return(id);
}

void close3D(nD* id){int i;
  for(i=0; i<id->planes; i++) close2D(id->d.z[i]);
  free(id->d.z); close_nD(id);}

float get4Dx(nD* id,unsigned int i,unsigned int j,unsigned int k,unsigned int l){
  return(((id->d.t[l])->get.n3D)(id->d.t[l],i,j,k));}
complex cget4Dx(nD* id,unsigned int i,unsigned int j,unsigned int k,unsigned int l){
  return(((id->d.t[l])->cget.n3D)(id->d.t[l],i,j,k));}
void badput4(nD* id,float value,unsigned int i,unsigned int j,
  unsigned int k,unsigned int l){
  err("Attempted to write to read-only file",id->name);}

nD* open4D(char* name){nD *id; char *name3D,*t; int i,j;
  id = open_nD (name);
  if(sscanf(id->prefix,"# 4D %5s %d %d %d %d",id->pixel,&(id->width),
    &(id->height),&(id->planes),&(id->times))!=5) err("Bad 4D file",id->name);
  id->get.n4D = get4Dx;
  id->cget.n4D = cget4Dx;
  id->put.n4D = badput4;
  if((id->d.t=malloc((id->times)*sizeof (nD*)))==0)
    perr("open4D malloc error",name);
  if((name3D=malloc(1+(id->notesoffset)-(id->dataoffset)))==0)
    perr("open4D malloc error",name);
  t = id->data;
  for(i=0; i<id->times; i++){
    for(j=0; t[0]!='\n'; j++) {name3D[j] = t[0]; t++;}
    name3D[j] = 0;  t++;
    id->d.t[i] = open3D(name3D);
  }
  free(name3D);
  return(id);}

nD* open3Dasnames(char *name)
{
  nD *id;
  char name2D[PATH_MAX];
  int i,j;
  unsigned char *t;

  id = open_nD (name);
  if(sscanf(id->prefix,"# 3D %5s %d %d %d",id->pixel,&(id->width),
    &(id->height),&(id->planes))!=4) err("Bad 3D file",id->name);
  id->get.n3D = get3Dx;
  id->cget.n3D = cget3Dx;
  id->put.n3D = badput3;
  if(!(id->nDnames=(char ***)malloc(sizeof(char **))))
    perr("open3Dasnames malloc error",name);
  if(!(id->nDnames[0]=(char **)malloc(id->planes*sizeof(char *))))
    perr("open3Dasnames malloc error",name);
  t = id->data;
  for(i=0; i<id->planes; i++){
    for(j=0; t<id->notes&&t[0]!='\n'; j++) {
      name2D[j] = t[0]; t++;
    }
    name2D[j] = 0; t++;
    if(!(id->nDnames[0][i]=(char *)malloc((strlen(name2D)+1)*sizeof(char))))
      perr("open3Dasnames malloc error",name);
    strcpy(id->nDnames[0][i],name2D);
  }
  return(id);
}

void close3Dasnames(nD* id){
  int i;
  for(i=0; i<id->planes; i++) {
    free(id->nDnames[0][i]);
  }
  free(id->nDnames[0]);
  free(id->nDnames);
  close_nD(id);
}


nD* open4Dasnames(char *name) 
{
  nD *id;
  nD *nD3D;
  char name3D[PATH_MAX];
  int i,j,k;
  char *t;
  id = open_nD (name);
  if(sscanf(id->prefix,"# 4D %5s %d %d %d %d",id->pixel,&(id->width),
    &(id->height),&(id->planes),&(id->times))!=5) err("Bad 4D file",id->name);
  id->get.n4D = get4Dx;
  id->cget.n4D = cget4Dx;
  id->put.n4D = badput4;
  if(!(id->nDnames3D=(char **)malloc((id->times)*sizeof(char *))))
    perr("open4Dasnames malloc error for nDnames3D",name);
  if(!(id->nDnames=(char ***)malloc((id->times)*sizeof(char **)))) 
    perr("open4Dasnames malloc error for nDnames",name);
  t = id->data;
  for(i=0; i<id->times; i++){
    for(j=0; t[0]!='\n'; j++) {name3D[j] = t[0]; t++;}
    name3D[j] = 0;  t++;
    if(!(id->nDnames3D[i]=(char *)malloc((j+1)*sizeof(char))))
      perr("open4Dasnames malloc error for nDnames3D[i]",name);
    strcpy(id->nDnames3D[i],name3D);
    nD3D = open3Dasnames(id->nDnames3D[i]);
    if(!(id->nDnames[i]=(char **)malloc((nD3D->planes)*sizeof(char *)))) 
      perr("open4Dasnames malloc error",name);
    for(k=0;k<nD3D->planes;k++) {
      if(!(id->nDnames[i][k]=
           (char *)malloc((strlen(nD3D->nDnames[0][k])+1)*sizeof(char)))) 
        perr("open4Dasnames malloc error",name);
      strcpy(id->nDnames[i][k],nD3D->nDnames[0][k]);
    }
    close3Dasnames(nD3D);
  }      
  return(id);
}  

void close4Dasnames(nD* id){
  int i,j;
  for(j=0;j<id->times;j++) {
    for(i=0; i<id->planes; i++) {
      free(id->nDnames[j][i]);
    }
    free(id->nDnames[j]);
    free(id->nDnames3D[j]);
  }
  free(id->nDnames);
  free(id->nDnames3D);
  close_nD(id);
}

float get4Dbyname(nD* id,unsigned int i,unsigned int j,unsigned int k,unsigned int l){
  nD *nD2D;
  float retval;

  nD2D = open2D(id->nDnames[l][k]);
  retval = get2D(nD2D,i,j);
  close2D(nD2D);
  return(retval);
}

float get3Dbyname(nD* id,unsigned int i,unsigned int j,unsigned int k){
  nD *nD2D;
  float retval;

  nD2D = open2D(id->nDnames[0][k]);
  retval = get2D(nD2D,i,j);
  close2D(nD2D);
  return(retval);
}

void close4D(nD* id){int i;
  for(i=0; i<id->times; i++) close3D(id->d.t[i]);
  free(id->d.t); close_nD(id);}
