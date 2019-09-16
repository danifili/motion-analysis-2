/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
/* MIT Shared memory extension */
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/time.h>
#include <signal.h>
#include <X11/extensions/XShm.h>
#include <X11/cursorfont.h>
#include <X11/extensions/Xrandr.h> //added 2010-11-07 (Scott) for Screen resolution detection
#include <linux/limits.h>
#include "vnd.h"
#include "nD.h"
#include "debug.h"
#include "ajerror.h"
#include "X11stuff.h"
#include "errhandler.h"
#include "read_qgm.h"

#define PI 3.14159265
#define CHANGETEXT 1
#define CREATEPOPUP 2
#define DESTROYPOPUP 3

#define MAXFRAMERATE 50

/* Segment size = length of dash and gap; line size = length of dash */
#define SEGMENT_SIZE 6
#define LINE_SIZE 2

#define VERT_TRIM 80  //Number of Pixels to allow for margins
#define WIDTH_TRIM 10


char *valid_keys="\n\
q     quits\n\
n,N   next pic (n stops at last pic, N goes on to first pic)\n\
p,P   previous pic (p stops at first pic, P goes on to last pic)\n\
a     autoscale display to each pic's min and max\n\
r     scale display to the global min and max\n\
g     set gamma value for display\n\
s     manually set scale for display\n\
t     set frame rate for showing multiple pictures\n\
f,F   show pictures continuously, moving forward\n\
b,B   show pictures continuously, moving backward\n\
z,x   set bottom and top of 3D roi, respectively\n\
e     change roifile\n\
keypad +        Next roi\n\
keypad -        Previous roi\n\
";

static XShmSegmentInfo shminfo;
int use_shm = 0;
extern X11Stuff xstf;
extern char *progname;
extern int debug;
extern double pixval_floor;
extern double pixval_slope;
extern double origmaxmax;
extern double origminmin;
extern int fullscreen;

/* Globals are necessary for alarm_handler */
extern int increment_pic;
X11Stuff *x_info_ptr;
int handlerwidth, handlerheight;
vndPic *vndpicptr;
int *picno;
int maxpicno;
int increment;
unsigned char *lookuptable;
int *mw, *mh;
int picheight,picwidth;
float widthratio, heightratio;


/* GetScreenSize function added by Scott 2010-11-07
 * based on code found at http://www.blitzbasic.com/Community/posts.php?topic=86911
 * which uses the libraries found in the original vnd, plus the librandr-dev library
 */

int GetScreenSize(int *width, int *height) {

    int num_sizes;
    Rotation original_rotation;

    Display *dpy = XOpenDisplay(NULL);
    Window root = RootWindow(dpy, 0);
    XRRScreenSize *xrrs = XRRSizes(dpy, 0, &num_sizes);

    XRRScreenConfiguration *conf = XRRGetScreenInfo(dpy, root);
    short original_rate          = XRRConfigCurrentRate(conf);
    SizeID original_size_id       = XRRConfigCurrentConfiguration(conf, &original_rotation);

    *width=xrrs[original_size_id].width - WIDTH_TRIM;
    *height=xrrs[original_size_id].height -VERT_TRIM;

    XCloseDisplay(dpy);
}


void killshmimages(void) {
    XShmDetach(xstf.theDisplay,&shminfo);
    XDestroyImage(xstf.theImage);
    shmdt(shminfo.shmaddr);
    shmctl(shminfo.shmid,IPC_RMID,0);
}

void handler(int errnum) {
    if(use_shm)
        killshmimages();
    if (errnum > 31)
        fprintf(stderr,"error %d received\n",errnum);
    else
        fprintf(stderr,"signal trapped: %s\n",errlist[errnum]);
    exit(errnum);
}

int load_pics(nD **pics, vndPic *vndpics, char **argv, int optind, int npics, int verticalmax, int widthmax, int savedots, char *dotfilename)
{
    int i,j,k;
    float min, max;
    char *tempstr;
    int iscplx=0;
    

    double *dx, *dy, *dz;
    char *startpos;
    double *xdispl, *ydispl, *zdispl;
    double xmag, xphase, ymag, yphase, zmag, zphase;    
    FILE *file;   
    int region,x1,y1,z1,x2,y2,z2;
    int roiregion[7];

    complex cval;

    pdebug(5,"load_pics()\n");
    for(k=0; k<npics; k++) {
        pdebug(10,"pic %d: %s\n",k,argv[optind+k]);
        /* Yep, that's a single equals-sign in the next line */
        if((vndpics[k].isqgm= is_qgm(argv[optind+k]))) {
            pdebug(5,"QGM\n");
            get_qgm_info(argv[optind+k],&(vndpics[k].width),&(vndpics[k].height));
            if (vndpics[k].width > widthmax && widthmax > 0) {
                vndpics[k].width = widthmax;
                ajerror("%s: Image width is larger than Screen/Commanded width and has been downsampled to fit. Use the -f option to display entire image!\n",progname);
            }
            if (vndpics[k].height > verticalmax && verticalmax > 0)
                vndpics[k].height = verticalmax;
            if(!(vndpics[k].pixvals=
                        (float **)malloc(vndpics[k].height*sizeof(float *)))) {
                ajerror("%s: couldn't malloc pointers to pic %d rows!\n",progname,k);
            }
            for(j=0; j<vndpics[k].height; j++) {
                if(!(vndpics[k].pixvals[j]=
                            (float *)malloc(vndpics[k].width*sizeof(float)))) {
                    ajerror("%s: couldn't malloc data for pic %d!\n",progname,k);
                }
            }
            if(read_qgm(vndpics[k].pixvals,argv[optind+k],
                        vndpics[k].width,vndpics[k].height)) {
                fprintf(stderr,"error reading QGM file!\n");
                return(-1);
            }
            min=max=vndpics[k].pixvals[0][0];
            for(j=0; j<vndpics[k].height; j++) {
                for(i=0; i<vndpics[k].width; i++) {
                    if(vndpics[k].pixvals[j][i] > max) max= vndpics[k].pixvals[j][i];
                    if(vndpics[k].pixvals[j][i] < min) min= vndpics[k].pixvals[j][i];
                }
            }
            if(!(vndpics[k].notes = (char *)malloc(10*sizeof(char)))) {
                ajerror("%s: couldn't malloc data for notes of pic %d\n",progname,k);
            }
            strcpy(vndpics[k].notes,"No Notes");
            pdebug(10,"      width=%d, height=%d, min=%f, max=%f\n",vndpics[k].width,
                   vndpics[k].height,vndpics[k].min,vndpics[k].max);

        } else { /* !is_qgm */
            pdebug(5,"NOT QGM\n");
            pics[k]=open2D(argv[optind+k]);
            if (pics[k]->width > widthmax && widthmax > 0) {

                fprintf(stderr,"%s: Image %d: width (%d pixels) > Screen/Commanded width (%d pixels) and has been downsampled to fit. Use -f option to display entire image!\n",progname,k, pics[k]->width,widthmax);
                //pics[k]->width = widthmax;
            }
            if (pics[k]->height > verticalmax && verticalmax > 0) {

                fprintf(stderr,"%s: Image %d: height (%d pixels) > Screen/Commanded height (%d pixels) and has been downsampled to fit. Use -f option to display entire image!\n",progname,k, pics[k]->height,verticalmax);
                //pics[k]->height = verticalmax;
            }
            if (verticalmax == 0)
                verticalmax = pics[k]->height;
            if (widthmax == 0)
                widthmax = pics[k]->width;
                 picheight = pics[k]->height;
            picwidth = pics[k]->width;

	    if (picheight < verticalmax)
                verticalmax = picheight;
            if (picwidth < widthmax)
                widthmax = picwidth;

            vndpics[k].width=widthmax;
            vndpics[k].height=verticalmax;
      
            if (picwidth == widthmax)
                 widthratio = 1;
            else
                 widthratio = (float)picwidth/(float)widthmax;
            if (picheight == verticalmax)
                 heightratio = 1;
            else
                 heightratio = (float)picheight/(float)verticalmax;
            

            
            printf("Pic width: %d, Max width: %d, Pic height: %d, Max height: %d\n",picwidth,widthmax,picheight,verticalmax);
            if(!(vndpics[k].pixvals=
                        (float **)malloc(vndpics[k].height*sizeof(float *)))) {
                ajerror("%s: couldn't malloc pointers to pic %d rows!\n",progname,k);
            }
            for(j=0; j<vndpics[k].height; j++) {
                if(!(vndpics[k].pixvals[j]=
                            (float *)malloc(vndpics[k].width*sizeof(float)))) {
                    ajerror("%s: couldn't malloc data for pic %d!\n",progname,k);
                }
            }
            if(!(strcmp(pics[k]->pixel,"C64"))) {
                pdebug(5,"C64\n");
                iscplx=1;
            }
            if(iscplx) {
                cval = cget2D(pics[k],0,0);
                min=max=sqrt(cval.r*cval.r+cval.i*cval.i);
            } else {
                min=max=get2D(pics[k],0,0);
            }/*
      for(j=0;j<pics[k]->height;j=j++) {
        for(i=0;i<pics[k]->width;i=i++) {
          if(iscplx) {
            cval = cget2D(pics[k],i,j);
            vndpics[k].pixvals[j][i] = sqrt(cval.r*cval.r+cval.i*cval.i);
          } else {
            vndpics[k].pixvals[j][i] = get2D(pics[k],i,j);
          }
          if(vndpics[k].pixvals[j][i]<min) {
            min=vndpics[k].pixvals[j][i];
            pdebug(99,"min %f at (%d,%d)\n",min,i,j);
          }
          if(vndpics[k].pixvals[j][i]>max) {
            max=vndpics[k].pixvals[j][i];
            pdebug(99,"max %f at (%d,%d)\n",max,i,j);
          }
        }
      }*/


            for(j=0; j<verticalmax; j++) {
                for(i=0; i<widthmax; i++) {
                    if(iscplx) {
                        cval = cget2D(pics[k],i,j);
                        vndpics[k].pixvals[j][i] = sqrt(cval.r*cval.r+cval.i*cval.i);
                    } else {


                        vndpics[k].pixvals[j][i] = get2D(pics[k],i* widthratio,j* heightratio);
                    }
                    if(vndpics[k].pixvals[j][i]<min) {
                        min=vndpics[k].pixvals[j][i];
                        pdebug(99,"min %f at (%d,%d)\n",min,i,j);
                    }
                    if(vndpics[k].pixvals[j][i]>max) {
                        max=vndpics[k].pixvals[j][i];
                        pdebug(99,"max %f at (%d,%d)\n",max,i,j);
                    }
                }
            }


            tempstr = getnotes(pics[k]);
            if(!(vndpics[k].notes = (char *)malloc((strlen(tempstr)+1)*sizeof(char)))) {
                ajerror("%s: couldn't malloc data for notes of pic %d\n",progname,k);
            }
            strcpy(vndpics[k].notes,tempstr);
            pdebug(10,"      width=%d, height=%d, min=%f, max=%f\n",vndpics[k].width,
                   vndpics[k].height,vndpics[k].min,vndpics[k].max);
            close2D(pics[k]);
        } /* if(is_qgm) ... else ... */
        vndpics[k].min = min;
        vndpics[k].max = max;
        vndpics[k].name = argv[optind+k];
    }
   region = 0;
 if(!(xdispl=(double *)malloc(npics*sizeof(double)))) {
            ajerror("couldn't malloc xdispl!\n");
        }
        if(!(ydispl=(double *)malloc(npics*sizeof(double)))) {
            ajerror("couldn't malloc ydispl!\n");
        }
        if(!(zdispl=(double *)malloc(npics*sizeof(double)))) {
            ajerror("couldn't malloc zdispl!\n");
        }




if (savedots == 1) {
    printf("\nSaving dot magnitude and phase to %s ...\n",dotfilename);
    if(!(dx=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dx\n");
    }
    if(!(dy=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dy\n");
    }
    if(!(dz=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dz\n");
    }
    /* Save dots information */

 

   
        if((file=fopen(dotfilename,"w"))==NULL) {
            ajerror("%s: couldn't open %s for writing",progname,dotfilename);
        }
      
        
    startpos = vndpics[0].notes;
    while((startpos=strchr(startpos,'r'))!=NULL) {
        if(!strncmp(startpos,"region ",7)) {
            if(sscanf(startpos,
                      "region %d (%d,%d,%d-%d,%d,%d): dx=%lg, dy=%lg, dz=%lg\n",
                      &region,&x1,&y1,&z1,&x2,&y2,&z2,dx,dy,dz)!=10) {
                pdebug(10,"2D data\n");
                
                sscanf(startpos,
                       "region %d (%d,%d-%d,%d): dx=%lg, dy=%lg\n",
                       &region,&x1,&y1,&x2,&y2,dx,dy);
                z1=z2=*dz=0;
                   
		   getmotions(vndpics, npics, x1+1, y1+1, &xmag, &xphase,  &ymag,  &yphase,  &zmag,&zphase,  roiregion);		
		fprintf(file, "Region: %d xmag: %f xphase: %f ymag: %f, yphase: %f\n",region, xmag,xphase, ymag,yphase);
            }           
        }
        startpos++;
    }
fclose(file);
    printf("%d Regions saved.\n",region);
}

    pdebug(5,"load_pics() done\n");
    return(0);
}

int compute_stats(vndPic *vndpics, int *maxwidth, int *maxheight,
                  double *minmin, double *maxmax, double svalue, double tvalue, int npics)
{
    int i;
    pdebug(5,"compute_stats()\n");
    if (tvalue >= 0)
	 *maxmax = tvalue;
    else
	  *maxmax = vndpics[0].max;
    if (svalue >= 0)
          *minmin = svalue;
    else
          *minmin = vndpics[0].min;

    *maxheight = vndpics[0].height;
    *maxwidth = vndpics[0].width;
    
    for(i=1; i<npics; i++) {
        if(vndpics[i].width > *maxwidth) *maxwidth = vndpics[i].width;
        if(vndpics[i].height > *maxheight) *maxheight = vndpics[i].height;
        if(vndpics[i].min < *minmin) {
		   if (svalue >= 0)
          		*minmin = svalue;
   		   else
          		*minmin = vndpics[i].min;
	}	
        if(vndpics[i].max > *maxmax) {
		 if (tvalue >= 0)
	 		*maxmax = tvalue;
    		else
	  		*maxmax = vndpics[i].max;
	}
    }
    pdebug(10,"maxwidth = %d, maxheight = %d, minmin = %f, maxmax = %f\n",
           *maxwidth, *maxheight, *minmin, *maxmax);
    pdebug(5,"compute_stats() done\n");
    return(0);
}

inline int pixelmap(double pixel)
{
    int i;
    i= (int)((pixel+pixval_floor)*pixval_slope);
    if(i<0) i= 0;
    if(i>LUT_SIZE-1) i= LUT_SIZE-1;
    return i;
}

int fill_lut(unsigned char *lut, double min, double max, double gamma,
             int ngrays)
{
    int i;
    double range;

    pdebug(5,"fill_lut()\n");
    range = pixelmap(max)-pixelmap(min);
    pdebug(5,"%f %f %f %d %d %d\n",min,max,range,pixelmap(min),pixelmap(max),
           pixelmap((max+min)/2));
    pdebug(6,"filling LUT below min\n");
    for(i=0; i<pixelmap(min); i++) {
        lut[i] = 0;
    }
    pdebug(6,"filling LUT\n");
    for(i=pixelmap(min); i<pixelmap(max); i++) {
        lut[i] = ngrays*pow((i-pixelmap(min))/range,1/gamma);
    }
    pdebug(6,"filling LUT above max\n");
    for(i=pixelmap(max); i<LUT_SIZE; i++) {
        lut[i] = ngrays-1;
    }
    pdebug(10,"min=%f, max=%f, ngrays=%d, gamma=%f\n",min,max,ngrays,gamma);
    pdebug(15,"value %f maps to %d\n",min+(max-min)/2,lut[(int)(pixelmap(min)/2+range/2)]);
    pdebug(5,"fill_lut() done\n");
    return(0);
}

int init_x_window(X11Stuff *xstf, int maxwidth, int maxheight,int useshm)
{
    int i;
    Cursor cur;
    XClassHint class_hints;
    char dpyvar[100],hostvar[100],tt[100];






    //maxwidth = 2452;
    //maxheight = 2056;

    pdebug(5,"init_x_window()\n");



    if((xstf->theDisplay=XOpenDisplay(NULL)) == NULL)
        ajerror("%s: can't open X server %s",progname,XDisplayName(NULL));
    pdebug(15,"Assigning screen\n");
    xstf->theScreen = DefaultScreen(xstf->theDisplay);

    pdebug(15,"Matching visual info\n");
    if(XMatchVisualInfo(xstf->theDisplay,xstf->theScreen,
                        8,PseudoColor,&(xstf->visual_info))) {
        pdebug(15,"Using pseudocolor\n");
        xstf->pseudocolor = 1;
    } else {
        if(!XMatchVisualInfo(xstf->theDisplay,xstf->theScreen,
                             24,TrueColor,&(xstf->visual_info))) {
            if(!XMatchVisualInfo(xstf->theDisplay,xstf->theScreen,
                                 32,TrueColor,&(xstf->visual_info))) {
                ajerror("This program only works in 8-bit PseudoColor and 24/32 bit TrueColor modes.\n");
            } else {
                pdebug(20,"32-bit\n");
            }
        } else {
            pdebug(20,"24-bit\n");
        }
        pdebug(15,"Using TrueColor\n");
        xstf->pseudocolor = 0;
    }
    pdebug(15,"Assigning visual\n");
    pdebug(15,"%x\n",xstf);
    pdebug(15,"%x\n",xstf->theVisual);
    pdebug(15,"%x\n",&(xstf->visual_info));
    pdebug(15,"%x\n",xstf->visual_info.visual);
    xstf->theVisual = xstf->visual_info.visual;
    pdebug(10,"Querying for shared memory extension\n");
    /* A MUCH better way to determine whether to use MIT SHM extensions */
    /* Unfortunately it doesn't work - it recognizes the SHM extension  */
    /* on remote machines, but can't access it and thereby fails.       */
    /* So I have to do this stupid fugly hack */
    if((getenv("DISPLAY") == NULL) || (getenv("HOSTNAME") == NULL)) {
        use_shm= 1;
    } else {
        strcpy(dpyvar,getenv("DISPLAY"));
        strcpy(hostvar,getenv("HOSTNAME"));
        if(!strcmp(dpyvar,":0.0")) use_shm=1;
        if(!strcmp(dpyvar,":0")) use_shm=1;
        if(!strcmp(dpyvar,"::")) use_shm=1;
        sprintf(tt,"%s:0.0",hostvar);
        if(!strcmp(dpyvar,tt)) use_shm=1;
        sprintf(tt,"%s:0",hostvar);
        if(!strcmp(dpyvar,tt)) use_shm=1;
        sprintf(tt,"%s::",hostvar);
        if(!strcmp(dpyvar,tt)) use_shm=1;
    }
    if (useshm == 0) {
    	use_shm=0; /* Set use_shm to 0 for remote X window xnd functionality */
    }
    /* Note: that really is supposed to be a single = on the next line  */
    if((use_shm)&&(use_shm= XShmQueryExtension(xstf->theDisplay))) {
        if(debug) fprintf(stderr,"using shared memory extensions\n");
        armhandler();

        pdebug(10,"Creating image\n");
        xstf->theImage =
            XShmCreateImage(xstf->theDisplay,xstf->theVisual,
                            xstf->visual_info.depth, ZPixmap, NULL, &(shminfo),
                            maxwidth, maxheight);
        if((shminfo.shmid= shmget(IPC_PRIVATE,((xstf->theImage)->bytes_per_line) *
                                  ((xstf->theImage)->height),IPC_CREAT|0777))== -1) {
            fprintf(stderr,"error in getting shared memory segment");
            killshmimages();
            exit(-1);
        }
        pdebug(10,"Allocating shared memory segment\n");
        if((shminfo.shmaddr = xstf->thePic = (xstf->theImage)->data =
                (char *)shmat(shminfo.shmid,0,0)) == 0) {
            fprintf(stderr,"couldn't allocate shared memory segment\n");
            killshmimages();
            exit(-1);
        }
        shminfo.readOnly = False;
        if(XShmAttach(xstf->theDisplay,&shminfo) == 0) {
            fprintf(stderr,"error attaching shared memory segment\n");
            killshmimages();
            exit(-1);
        }
    } else { /* Not using shared memory extensions */
        pdebug(10,"Not using shared memory\n");
        xstf->theImage =
            XCreateImage(xstf->theDisplay, xstf->theVisual,
                         xstf->visual_info.depth,
                         ZPixmap,0, xstf->thePic, maxwidth, maxheight, 8, 0);
        if(!(xstf->thePic=(unsigned char *)
                          malloc(maxheight*xstf->theImage->bytes_per_line)))
            ajerror("%s: Couldn't allocate memory for thePic",progname);
        xstf->theImage->data=xstf->thePic;
    }
    xstf->theAttributes.backing_store = Always;
    if (fullscreen == 1) {
        xstf->theAttributes.override_redirect = True;
    }
//xstf->theAttributes.event_mask = ResizeRedirectMask|KeyPressMask;
    pdebug(5,"Creating X Window\n");


    printf("Screen Width: %d Screen Height: %d\n", DisplayWidth(xstf->theDisplay, xstf->theScreen), DisplayHeight(xstf->theDisplay, xstf->theScreen));
   
    xstf->theWindow =
        XCreateWindow(xstf->theDisplay,RootWindow(xstf->theDisplay,xstf->theScreen),
                      0,0,maxwidth,maxheight,10,
                      DefaultDepth(xstf->theDisplay,xstf->theScreen),InputOutput,
                      xstf->theVisual,CWOverrideRedirect|CWBackingStore,&(xstf->theAttributes));
    xstf->theSizeHints.flags = PPosition | PSize;
    xstf->theSizeHints.x = xstf->theSizeHints.y = 0;
    xstf->theSizeHints.width = maxwidth;
    xstf->theSizeHints.height = maxheight;
    XSetStandardProperties(xstf->theDisplay,xstf->theWindow,
                           xstf->WindowTitle,xstf->WindowTitle,
                           None,0,0,&(xstf->theSizeHints));
    XSelectInput(xstf->theDisplay,xstf->theWindow,THEMASK);
    class_hints.res_class = "Untitled";
    class_hints.res_name = "vnd";
    XSetClassHint(xstf->theDisplay,xstf->theWindow,&class_hints);
    cur = XCreateFontCursor(xstf->theDisplay,XC_draft_small);
    XDefineCursor(xstf->theDisplay,xstf->theWindow,cur);
    XMapWindow(xstf->theDisplay,xstf->theWindow);
    XMoveWindow(xstf->theDisplay,xstf->theWindow,5,24);
    i=100;
    if (fullscreen == 1) {
        XGrabKeyboard(xstf->theDisplay, xstf->theWindow, 0, GrabModeAsync, GrabModeAsync, CurrentTime);
    }
    while((XGrabPointer(xstf->theDisplay,xstf->theWindow,False,0,GrabModeAsync,
                        GrabModeAsync,xstf->theWindow,None,CurrentTime)
            != GrabSuccess) && i!=0) {
        usleep(100);
        i--;
    }
    if(i==0)  {
        fprintf(stderr,"Couldn't grab pointer\n");
    }
    XUngrabPointer(xstf->theDisplay,CurrentTime);
    pdebug(5,"init_x_window() done\n");
    return(0);
}

int xcolors(X11Stuff *xstf, int color)
{
    int i,filled[256];
    static int order[240]= {
        0,239,119,59,179,29,149,89,209,14,104,194,44,134,224,74,164,6,81,156,231,
        21,96,171,36,111,186,51,126,201,66,141,216,3,70,138,205,33,100,168,235,17,
        85,153,220,48,115,183,9,78,145,213,63,130,198,25,93,160,228,55,123,190,40,
        108,175,2,61,121,181,57,117,177,237,20,80,140,200,11,72,132,192,38,98,158,
        218,8,68,128,188,23,83,143,203,35,95,155,215,5,65,125,185,16,76,136,196,50,
        110,170,230,27,87,147,207,31,91,151,211,46,106,166,226,53,113,173,233,42,
        102,162,222,34,69,103,137,172,206,32,67,101,135,169,204,238,15,49,84,118,
        152,187,221,7,41,75,109,144,178,212,24,58,92,127,161,195,229,12,45,79,114,
        148,182,217,19,54,88,122,157,191,225,28,62,97,131,165,199,234,1,37,71,105,
        139,174,208,13,47,82,116,150,184,219,22,56,90,124,159,193,227,30,64,99,133,
        167,202,236,4,39,73,107,142,176,210,18,52,86,120,154,189,223,26,60,94,129,
        163,197,232,10,43,77,112,146,180,214
    };
    static int reds[]= {0,0,0,0,2,2,2,2,1,1,1,1,3,3,3,3};
    static int greens[]= {0,0,2,2,0,0,1,2,1,1,3,3,1,1,3,3};
    static int blues[]= {0,2,0,2,0,2,0,2,1,3,1,3,1,3,1,3};

    pdebug(5,"xcolors()\n");
    for(i=0; i<256; i++) filled[i]=0;
    for(i=0; i<240; i++) {
        xstf->cdefs[i].flags = DoRed | DoGreen | DoBlue;
        xstf->cdefs[i].blue=xstf->cdefs[i].red=xstf->cdefs[i].green=65535.*i/239;
    }
    for(i=240; i<256; i++) {
        xstf->cdefs[i].flags = DoRed | DoGreen | DoBlue;
        xstf->cdefs[i].red=21845*reds[i-240];
        xstf->cdefs[i].green=21845*greens[i-240];
        xstf->cdefs[i].blue=21845*blues[i-240];
    }
    for(i=240; i<256; i++) {
        if(!XAllocColor(xstf->theDisplay,DefaultColormap(xstf->theDisplay,
                        xstf->theScreen),
                        &((xstf->cdefs)[i])))
            ajwarn("%s: Couldn't allocate 1 of 16 colors in color map",progname);
        else filled[i]=1;
    }
    for(i=0; i<240; i++)
        if(XAllocColor(xstf->theDisplay,DefaultColormap(xstf->theDisplay,
                       xstf->theScreen),
                       &((xstf->cdefs)[order[i]])))
            filled[order[i]]=1;
    for(i=xstf->ngrays=0; i<240; i++)
        if(filled[i]&&(xstf->cdefs[i].pixel!=xstf->cdefs[xstf->ngrays].pixel))
            xstf->cdefs[(xstf->ngrays)++]=xstf->cdefs[i];
    for(i=0; i<xstf->ngrays; i++) xstf->pixels[i]=xstf->cdefs[i].pixel;
    for(; i<240; i++) {
        xstf->pixels[i]=xstf->cdefs[xstf->ngrays-1].pixel;
    }
    for(i=240; i<256; i++)  {
        if(filled[i]) xstf->pixels[i]=xstf->cdefs[i].pixel;
        else xstf->pixels[i] = xstf->cdefs[xstf->ngrays-1].pixel;
    }
    if(color) {
        xstf->pixels[0]=xstf->pixels[GREEN];
        xstf->pixels[++(xstf->ngrays)]=xstf->pixels[RED];
        xstf->pixels[xstf->ngrays+1]=xstf->pixels[RED];
    }
    if(xstf->visual_info.class==PseudoColor) {
        if(xstf->ngrays<100) {
            ajwarn("only %d grays available\n",xstf->ngrays);
        }
    } else {
        xstf->ngrays = 256; /* All grays available in TrueColor */
    }
    pdebug(5,"xcolors() done\n");
    return(0);
}

int map_pic_to_image(X11Stuff *xstf,vndPic *vndpics,int picnum,
                     unsigned char *lut, int maxwidth, int maxheight)
{
    int i,j;

    pdebug(5,"map_pic_to_image()\n");
    pdebug(10,"%x\n",xstf->thePic);
    if(xstf->visual_info.class==PseudoColor) {
        if(!(vndpics[picnum].image_mapped)) {
            for(j=0; j<vndpics[picnum].height; j++) {
                for(i=0; i<vndpics[picnum].width; i++) {
                    //          vndpics[picnum].image[maxwidth*j+i] =
                    vndpics[picnum].image[(xstf->theImage)->bytes_per_line*j+i] =
                        xstf->pixels[(unsigned int)lut[pixelmap(vndpics[picnum].pixvals[vndpics[picnum].height-1-j][i])]];
                    //            lut[(int)vndpics[picnum].pixvals[j][i]];
                }
                memset((void *)(vndpics[picnum].image+(xstf->theImage)->bytes_per_line*j+vndpics[picnum].width),
                       xstf->pixels[lut[0]],(xstf->theImage)->bytes_per_line-vndpics[picnum].width);
            }
            memset((void *)(vndpics[picnum].image+(xstf->theImage)->bytes_per_line*vndpics[picnum].height),
                   xstf->pixels[lut[0]],(xstf->theImage)->bytes_per_line*(maxheight-vndpics[picnum].height));
            vndpics[picnum].image_mapped = 1;
        }
    } else { /* TrueColor */
        if(!(vndpics[picnum].image_mapped)) {
            for(j=0; j<vndpics[picnum].height; j++) {
                for(i=0; i<vndpics[picnum].width; i++) {
                    memset((void *)(vndpics[picnum].image+(xstf->theImage)->bytes_per_line*j+(xstf->theImage)->bits_per_pixel/8*i),
                           lut[pixelmap(vndpics[picnum].pixvals[vndpics[picnum].height-1-j][i])],3);
                }
                memset((void *)(vndpics[picnum].image+(xstf->theImage)->bytes_per_line*j+
                                (xstf->theImage)->bits_per_pixel/8*vndpics[picnum].width),
                       0,(xstf->theImage)->bytes_per_line-(xstf->theImage)->bits_per_pixel/8*vndpics[picnum].width);
            }
            for(j=vndpics[picnum].height; j<maxheight; j++) {
                memset((void *)(vndpics[picnum].image+(xstf->theImage)->bytes_per_line*j),0,(xstf->theImage)->bytes_per_line);
            }
            vndpics[picnum].image_mapped = 1;
        }
    }
    memcpy(xstf->thePic,vndpics[picnum].image,
           maxheight*((xstf->theImage)->bytes_per_line));
    pdebug(5,"map_pic_to_image() done\n");
    return(0);
}

int draw_image_to_screen(X11Stuff *xstf,int maxwidth,int maxheight,char *name)
{
    pdebug(5,"draw_image_to_screen()\n");
    if(use_shm) {
        pdebug(15,"XShmPutImage\n");
        XShmPutImage(xstf->theDisplay,xstf->theWindow,
                     DefaultGC(xstf->theDisplay,xstf->theScreen),
                     xstf->theImage,0,0,0,0,maxwidth,maxheight,False);
    } else {
        pdebug(15,"XPutImage\n");
        XPutImage(xstf->theDisplay,xstf->theWindow,
                  DefaultGC(xstf->theDisplay,xstf->theScreen),
                  xstf->theImage,0,0,0,0,maxwidth,maxheight);
    }
    if(name!=NULL) {
        pdebug(15,"XStoreName\n");
        XStoreName(xstf->theDisplay,xstf->theWindow,name);
    }
    pdebug(15,"XFlush\n");
    XFlush(xstf->theDisplay);
    pdebug(5,"draw_image_to_screen() done\n");
    return(0);
}

/* This code was lifted directly from qvs, so don't blame me.  -A.J. */
int infowindow(X11Stuff *xstf, int command, char text[], int xpos, int ypos)
{
    /* pop up window.*/
    static int open=0;
    static Window popup;
    static GC gc;
    static int width, height;
    static XFontStruct *font_info=NULL;
    char junk[2000];
    Window junk_win;
    int xroot, yroot;
    int i;
    char *foo,*startofline;
    int maxline=0;
    int nlines=0;

    pdebug(5,"infowindow()\n");
    if(command==DESTROYPOPUP) {
        pdebug(10,"DESTROYPOPUP\n");
        if(open) {
            XDestroyWindow(xstf->theDisplay,popup);
            open=0;
            return(0);
        } else {
            fprintf(stderr,"infowindow: can't destroy window, popup isn't open\n");
            return(-1);
        }
    } else if (command==CHANGETEXT) {
        pdebug(10,"CHANGETEXT\n");
        if(open) {
            sprintf(junk,"%s \n",text);
            startofline=junk;
            nlines=0;
            while((foo=strchr(startofline,'\n'))!=NULL) {
                XDrawImageString(xstf->theDisplay,popup,gc,font_info->max_bounds.width,
                                 (int)((1.2)* font_info->max_bounds.ascent+nlines*
                                       (font_info->max_bounds.ascent+font_info->max_bounds.descent)),
                                 startofline,foo-startofline);
                startofline=foo+1;
                nlines++;
            }
            //      XDrawImageString(xstf->theDisplay,popup,gc,font_info->max_bounds.width,
            //                       (int)(1.2* font_info->max_bounds.ascent),
            //                       text,strlen(text));
            for(i=0; i<8; i++) {
                XStoreBuffer(xstf->theDisplay,junk,strlen(junk),i);
            }
            XStoreBytes(xstf->theDisplay, junk, strlen(junk));
            XTranslateCoordinates(xstf->theDisplay,xstf->theWindow,
                                  RootWindow(xstf->theDisplay,xstf->theScreen),
                                  xpos,ypos,&xroot,&yroot,&junk_win);
            if(xroot+width>DisplayWidth(xstf->theDisplay,xstf->theScreen))
                xroot=DisplayWidth(xstf->theDisplay,xstf->theScreen)-width;
            if(yroot<0) yroot=0;
            XMoveWindow(xstf->theDisplay,popup,xroot,yroot);
            return(0);
        }
        else {
            fprintf(stderr,"infowindow: can't change text, popup isn't open\n");
            return(-1);
        }
    } else if (command==CREATEPOPUP) {
        pdebug(10,"CREATEPOPUP\n");
        if(!open) {
            XSetWindowAttributes att;
            char *fontname1 ="10x20", *fontname2 ="9x15", *fontname3 ="fixed";
            char *userfont;
            XGCValues *gcval=NULL;

            XTranslateCoordinates(xstf->theDisplay,xstf->theWindow,
                                  RootWindow(xstf->theDisplay,xstf->theScreen),
                                  xpos,ypos,&xroot,&yroot,&junk_win);
            open=1;
            att.override_redirect=1;
            att.background_pixel=BlackPixel(xstf->theDisplay,xstf->theScreen);
            att.event_mask = ResizeRedirectMask|KeyPressMask;
            if((userfont=XGetDefault(xstf->theDisplay,"qvs","font"))!=NULL)
                font_info = XLoadQueryFont(xstf->theDisplay,userfont);
            if(font_info==NULL)
                if ((font_info=XLoadQueryFont(xstf->theDisplay,fontname1)) == NULL)
                    if ((font_info=XLoadQueryFont(xstf->theDisplay,fontname2)) == NULL)
                        if ((font_info=XLoadQueryFont(xstf->theDisplay,fontname3)) == NULL) {
                            ajerror("Couldn't load any of the following fonts %s %s %s\n",
                                    fontname1, fontname2, fontname3);
                            exit(-1);
                        }
            sprintf(junk,"%s \n",text);
            startofline=junk;
            while((foo=strchr(startofline,'\n'))!=NULL) {
                if(maxline<foo-startofline) maxline=foo-startofline;
                startofline=foo+1;
                nlines++;
            }
            width=(3+maxline)* font_info->max_bounds.width;
            height=(nlines+0.5)*(font_info->max_bounds.ascent+font_info->max_bounds.descent);
            if(xroot+width>DisplayWidth(xstf->theDisplay,xstf->theScreen))
                xroot=DisplayWidth(xstf->theDisplay,xstf->theScreen)-width;
            if(yroot<0) yroot=0;
            popup=XCreateWindow(xstf->theDisplay,
                                RootWindow(xstf->theDisplay,xstf->theScreen),
                                xroot,yroot,width,height,0,
                                DefaultDepth(xstf->theDisplay,xstf->theScreen),
                                InputOutput,xstf->theVisual,
                                CWBackPixel|CWOverrideRedirect|CWEventMask|CWBorderPixel,&att);
            gc = XCreateGC(xstf->theDisplay,popup,0,gcval);
            XSetFont(xstf->theDisplay,gc,font_info->fid);
            XSetForeground(xstf->theDisplay,gc,
                           WhitePixel(xstf->theDisplay,xstf->theScreen));
            XSetBackground(xstf->theDisplay,gc,
                           BlackPixel(xstf->theDisplay,xstf->theScreen));
            XMapWindow(xstf->theDisplay,popup);
            startofline=junk;
            nlines=0;
            while((foo=strchr(startofline,'\n'))!=NULL) {
                XDrawImageString(xstf->theDisplay,popup,gc,font_info->max_bounds.width,
                                 (int)((1.2)* font_info->max_bounds.ascent+nlines*
                                       (font_info->max_bounds.ascent+font_info->max_bounds.descent)),
                                 startofline,foo-startofline);
                startofline=foo+1;
                nlines++;
            }
            //      sprintf(junk,"%s \n",text);
            for(i=0; i<8; i++) {
                XStoreBuffer(xstf->theDisplay,junk,strlen(junk),i);
            }
            XStoreBytes(xstf->theDisplay,junk,strlen(junk));
            return(0);
        } else {
            fprintf(stderr,"infowindow: can't create window, popup already open\n");
            return(-1);
        }
    } else {
        fprintf(stderr,"infowindow: don't know command %d\n",command);
        return(-1);
    }
    pdebug(5,"infowindow() done\n");
    return(0);
}

int makecursor(X11Stuff *xstf, int x, int y, int maxheight)
{
    pdebug(5,"makecursor()\n");
    y=maxheight-y;
    XDrawLine(xstf->theDisplay,xstf->theWindow,
              DefaultGC(xstf->theDisplay,xstf->theScreen),x+2,y,x+6,y);
    XDrawLine(xstf->theDisplay,xstf->theWindow,
              DefaultGC(xstf->theDisplay,xstf->theScreen),x-2,y,x-6,y);
    XDrawLine(xstf->theDisplay,xstf->theWindow,
              DefaultGC(xstf->theDisplay,xstf->theScreen),x,y+2,x,y+6);
    XDrawLine(xstf->theDisplay,xstf->theWindow,
              DefaultGC(xstf->theDisplay,xstf->theScreen),x,y-2,x,y-6);
    XFlush(xstf->theDisplay);
    pdebug(5,"makecursor() done\n");
    return(0);
}

int makesmallcursor(X11Stuff *xstf, int x, int y, int maxheight)
{
    pdebug(5,"makesmallcursor()\n");
    XSetForeground(xstf->theDisplay,
                   DefaultGC(xstf->theDisplay,xstf->theScreen),
                   xstf->pixels[LT_BLUE]);
    y=maxheight-y;
    XDrawLine(xstf->theDisplay,xstf->theWindow,
              DefaultGC(xstf->theDisplay,xstf->theScreen),x-6,y,x+6,y);
    XDrawLine(xstf->theDisplay,xstf->theWindow,
              DefaultGC(xstf->theDisplay,xstf->theScreen),x,y-6,x,y+6);
    XSetForeground(xstf->theDisplay,
                   DefaultGC(xstf->theDisplay,xstf->theScreen),
                   xstf->pixels[YELLOW]);
    XFlush(xstf->theDisplay);
    pdebug(5,"makesmallcursor() done\n");
    return(0);
}

int draw_roi(X11Stuff *xstf, Roi roi, int picnum, int maxwidth, int maxheight)
{
    static XSegment *segments=NULL;
    int nsegments=0;
    int i;

    pdebug(5,"draw_roi()\n");
    if(segments==NULL) {
        if(!(segments=malloc((MAX(maxwidth,maxheight)/SEGMENT_SIZE+1)*
                             sizeof(XSegment)))) {
            ajerror("%s: error mallocing line segments",progname);
        }
    }
    if((picnum<roi.zmin)||(picnum>roi.zmax)) { /* If not within z span */
        nsegments=(roi.xmax  - roi.xmin)/ widthratio /SEGMENT_SIZE+1;
        /* Draw bottom horiz line */
        for(i=0; i<nsegments; i++) {
            segments[i].x1 = roi.xmin / widthratio  + SEGMENT_SIZE*i;
            segments[i].x2 = roi.xmin / widthratio + SEGMENT_SIZE*i + LINE_SIZE-1;
            segments[i].y1 = (picheight - roi.ymin) / heightratio;
            segments[i].y2 = (picheight - roi.ymin) / heightratio;
        }
        XDrawSegments(xstf->theDisplay,xstf->theWindow,
                      DefaultGC(xstf->theDisplay,xstf->theScreen),
                      segments,nsegments);
        /* Draw top horiz line */
        for(i=0; i<nsegments; i++) {
            segments[i].x1 = roi.xmin / widthratio + SEGMENT_SIZE*i;
            segments[i].x2 = roi.xmin / widthratio + SEGMENT_SIZE*i + LINE_SIZE-1;
            segments[i].y1 = (picheight - roi.ymax) / heightratio;
            segments[i].y2 = (picheight - roi.ymax) / heightratio;
        }
        XDrawSegments(xstf->theDisplay,xstf->theWindow,
                      DefaultGC(xstf->theDisplay,xstf->theScreen),
                      segments,nsegments);
        nsegments=(roi.ymax - roi.ymin)/SEGMENT_SIZE+1;
        /* Draw left vert line */
        for(i=0; i<nsegments; i++) {
            segments[i].x1 = roi.xmin / widthratio;
            segments[i].x2 = roi.xmin / widthratio;
            segments[i].y1 = (picheight  - (roi.ymin  + SEGMENT_SIZE*i)) / heightratio ;
            segments[i].y2 = (picheight - (roi.ymin  + SEGMENT_SIZE*i + LINE_SIZE-1)) / heightratio ;
        }
        XDrawSegments(xstf->theDisplay,xstf->theWindow,
                      DefaultGC(xstf->theDisplay,xstf->theScreen),
                      segments,nsegments);
        /* Draw right vert line */
        for(i=0; i<nsegments; i++) {
            segments[i].x1 = roi.xmax / widthratio;
            segments[i].x2 = roi.xmax / widthratio;
            segments[i].y1 = (picheight  - (roi.ymin + SEGMENT_SIZE*i)) / heightratio ;
            segments[i].y2 = (picheight  - (roi.ymin  + SEGMENT_SIZE*i + LINE_SIZE-1)) / heightratio ;
        }
        XDrawSegments(xstf->theDisplay,xstf->theWindow,
                      DefaultGC(xstf->theDisplay,xstf->theScreen),
                      segments,nsegments);
    } else { /* Within z span, so can draw solid lines */
        XDrawLine(xstf->theDisplay,xstf->theWindow,
                  DefaultGC(xstf->theDisplay,xstf->theScreen),
                  roi.xmin / widthratio,(picheight-roi.ymin) / heightratio,roi.xmax / widthratio,(picheight-roi.ymin) / heightratio);
        XDrawLine(xstf->theDisplay,xstf->theWindow,
                  DefaultGC(xstf->theDisplay,xstf->theScreen),
                  roi.xmin / widthratio,(picheight-roi.ymax) / heightratio,roi.xmax / widthratio,(picheight-roi.ymax) / heightratio);
        XDrawLine(xstf->theDisplay,xstf->theWindow,
                  DefaultGC(xstf->theDisplay,xstf->theScreen),
                  roi.xmin / widthratio,(picheight-roi.ymin) / heightratio,roi.xmin / widthratio,(picheight-roi.ymax) / heightratio);
        XDrawLine(xstf->theDisplay,xstf->theWindow,
                  DefaultGC(xstf->theDisplay,xstf->theScreen),
                  roi.xmax / widthratio,(picheight-roi.ymin) / heightratio,roi.xmax / widthratio,(picheight-roi.ymax) / heightratio);
    }
    XFlush(xstf->theDisplay);
    pdebug(5,"draw_roi() done\n");
    return(0);
}

int draw_rois(X11Stuff *xstf, Roi *rois, int nrois, int currentroi, int picnum,
              int maxwidth, int maxheight)
{
    int i;

    pdebug(5,"draw_rois()\n");
    XSetForeground(xstf->theDisplay,
                   DefaultGC(xstf->theDisplay,xstf->theScreen),
                   xstf->pixels[LT_BLUE]);
    for(i=0; i<nrois; i++) {
        if(i==currentroi) {
            pdebug(5,"currentroi=%d\n",i);
            XSetForeground(xstf->theDisplay,
                           DefaultGC(xstf->theDisplay,xstf->theScreen),
                           xstf->pixels[YELLOW]);
        }
        pdebug(5,"rois[i].inuse = %d\n",rois[i].inuse);
        if(rois[i].inuse) {
            draw_roi(xstf,rois[i],picnum,maxwidth,maxheight);
        }
        if(i==currentroi) {
            XSetForeground(xstf->theDisplay,
                           DefaultGC(xstf->theDisplay,xstf->theScreen),
                           xstf->pixels[LT_BLUE]);
        }
    }
    pdebug(5,"draw_rois() done\n");
    return(0);
}

int redisplay_pic(X11Stuff *xstf, int autoscale, unsigned char *lut,
                  double minmin, double maxmax, double gamma, int picnum,
                  int maxwidth, int maxheight,vndPic *vndpics,
                  Roi *rois, int nrois, int currentroi)
{
    int retval=0;
    static int ntimes=0;

    ntimes++;
    pdebug(5,"redisplay_pic()\n");
    if(!(vndpics[picnum].image_mapped)) {
        if(autoscale) {
            pdebug(10,"calling fill_lut()\n");
            retval |= fill_lut(lut,vndpics[picnum].min,vndpics[picnum].max,gamma,xstf->ngrays);
        }
    }
    pdebug(10,"calling map_pic_to_image()\n");
    retval |= map_pic_to_image(xstf,vndpics,picnum,lut,maxwidth,maxheight);
    pdebug(10,"calling draw_image_to_screen()\n");
    retval |= draw_image_to_screen(xstf,maxwidth,maxheight,vndpics[picnum].name);
    retval |= draw_rois(xstf,rois,nrois,currentroi,picnum,maxwidth,maxheight);
    pdebug(5,"redisplay_pic() done %d\n",ntimes);
    if (fullscreen == 1) {
        XGrabKeyboard(xstf->theDisplay, xstf->theWindow, 0, GrabModeAsync, GrabModeAsync, CurrentTime);
    }
    return(retval);
}

/* Another one lifted straight from Q's code */
void getvalue(X11Stuff *xstf, char text[], double *ans)
{
    Window popup;
    XSetWindowAttributes att;
    int width = 620;
    int height = 100, i;
    char echo[200],ansstring[200];
    XFontStruct *font_info;
    char *fontname1 ="10x20", *fontname2 ="9x15", *fontname3 ="fixed";
    char *userfont;
    GC gc;
    XGCValues *gcval=NULL;
    static int firsttime =1;
    static XColor cdefs;
    static unsigned long pixel;
    static int reds[]= {0,0,0,0,2,2,2,2,1,1,1,1,3,3,3,3};
    static int greens[]= {0,0,2,2,0,0,1,2,1,1,3,3,1,1,3,3};
    static int blues[]= {0,2,0,2,0,2,0,2,1,3,1,3,1,3,1,3};

    pdebug(5,"getvalue()\n");
    if(firsttime) {
        firsttime=0;
        cdefs.red=21845*reds[244-240];
        cdefs.green=21845*greens[244-240];
        cdefs.blue=21845*blues[244-240];

        if(!XAllocColor(xstf->theDisplay,
                        DefaultColormap(xstf->theDisplay,xstf->theScreen),&cdefs))
            pixel=WhitePixel(xstf->theDisplay,xstf->theScreen);
        else pixel = cdefs.pixel;
    }

    ansstring[0]=0;
    att.override_redirect=1;
    att.background_pixel=BlackPixel(xstf->theDisplay,xstf->theScreen);
    att.border_pixel=pixel;
    att.event_mask = ResizeRedirectMask|KeyPressMask;
    popup=
        XCreateWindow(xstf->theDisplay,RootWindow(xstf->theDisplay,xstf->theScreen),
                      (DisplayWidth(xstf->theDisplay,xstf->theScreen)-width)/2,
                      (DisplayHeight(xstf->theDisplay,xstf->theScreen)-height)/2,
                      width,height,10,DefaultDepth(xstf->theDisplay,xstf->theScreen),
                      InputOutput,xstf->theVisual,
                      CWBackPixel|CWOverrideRedirect|CWEventMask|CWBorderPixel,
                      &att);
    font_info=NULL;
    if((userfont=XGetDefault(xstf->theDisplay,"qvs","font"))!=NULL)
        font_info = XLoadQueryFont(xstf->theDisplay,userfont);
    if(font_info==NULL)
        if ((font_info = XLoadQueryFont(xstf->theDisplay,fontname1)) == NULL)
            if ((font_info = XLoadQueryFont(xstf->theDisplay,fontname2)) == NULL)
                if ((font_info = XLoadQueryFont(xstf->theDisplay,fontname3)) == NULL) {
                    fprintf(stderr,
                            "Couldn't load any of the following fonts %s %s %s\n",
                            fontname1, fontname2, fontname3);
                    exit(-2);
                }
    gc = XCreateGC(xstf->theDisplay,popup,0,gcval);
    XSetFont(xstf->theDisplay,gc,font_info->fid);
    XSetForeground(xstf->theDisplay,gc,
                   WhitePixel(xstf->theDisplay,xstf->theScreen));
    XSetBackground(xstf->theDisplay,gc,
                   BlackPixel(xstf->theDisplay,xstf->theScreen));
    XMapWindow(xstf->theDisplay,popup);
    if(XGrabKeyboard(xstf->theDisplay,popup,False,GrabModeAsync,GrabModeAsync,
                     CurrentTime) != GrabSuccess) {
        fprintf(stderr,"Couldn't grab keyboard");
        exit(-88);
    }

    {   i=100;
        while((XGrabPointer(xstf->theDisplay,popup,False,0,
                            GrabModeAsync,
                            GrabModeAsync,popup,None,CurrentTime)
                != GrabSuccess)
                && i!=0)    i--;
        if(i==0)  {
            fprintf(stderr,"Couldn't grab pointer\n");
        }

        XUngrabPointer(xstf->theDisplay,CurrentTime);
    }
    sprintf(echo,"%s (%g):",text,*ans);
    XDrawImageString(xstf->theDisplay,popup,gc,10,50,echo,strlen(echo));
    while(1) {
        XEvent e;
        char key[10];
        KeySym keysym;
        int length;
        for(i=0; i<sizeof(key); i++) key[i]='\0';
        XNextEvent(xstf->theDisplay,&e);
        switch(e.type) {
        case KeyPress:
            if(XLookupString((XKeyEvent *)&e,key,sizeof(key),&keysym,NULL)) {
                if ((keysym==XK_Return)||(keysym==XK_KP_Enter)||(keysym==XK_Linefeed)||
                        (keysym==XK_space)) {
                    if((length=strlen(ansstring))!=0) *ans = atof(ansstring);
                    XDestroyWindow(xstf->theDisplay,popup);
                    return;
                }
                else if (keysym==XK_Escape) {
                    XDestroyWindow(xstf->theDisplay,popup);
                    return;
                }
                else if (((keysym <= XK_9) && (keysym >= XK_0)) ||
                         (keysym == XK_period)||(keysym==XK_minus))
                    strcat(ansstring,key);
                else if ((keysym==XK_BackSpace)||(keysym==XK_Delete)) {
                    if((length=strlen(ansstring))>0) {
                        ansstring[length-1]='\0';
                        XClearWindow(xstf->theDisplay,popup);
                    }
                    else XBell(xstf->theDisplay,10);
                }
                else XBell(xstf->theDisplay,10);
                sprintf(echo,"%s (%g):%s",text,*ans,ansstring);
                XDrawImageString(xstf->theDisplay,popup,gc,10,50,echo,strlen(echo));
            }
        }
    }
    pdebug(5,"getvalue() done\n");
}

void getstring(X11Stuff *xstf, char text[], char *ans)
{
    Window popup;
    XSetWindowAttributes att;
    int width = 620;
    int height = 100, i;
    char echo[200],ansstring[200];
    XFontStruct *font_info;
    char *fontname1 ="10x20", *fontname2 ="9x15", *fontname3 ="fixed";
    char *userfont;
    GC gc;
    XGCValues *gcval=NULL;
    static int firsttime =1;
    static XColor cdefs;
    static unsigned long pixel;
    static int reds[]= {0,0,0,0,2,2,2,2,1,1,1,1,3,3,3,3};
    static int greens[]= {0,0,2,2,0,0,1,2,1,1,3,3,1,1,3,3};
    static int blues[]= {0,2,0,2,0,2,0,2,1,3,1,3,1,3,1,3};

    pdebug(5,"getstring()\n");
    if(firsttime) {
        firsttime=0;
        cdefs.red=21845*reds[244-240];
        cdefs.green=21845*greens[244-240];
        cdefs.blue=21845*blues[244-240];

        if(!XAllocColor(xstf->theDisplay,
                        DefaultColormap(xstf->theDisplay,xstf->theScreen),&cdefs))
            pixel=WhitePixel(xstf->theDisplay,xstf->theScreen);
        else pixel = cdefs.pixel;
    }

    ansstring[0]=0;
    att.override_redirect=1;
    att.background_pixel=BlackPixel(xstf->theDisplay,xstf->theScreen);
    att.border_pixel=pixel;
    att.event_mask = ResizeRedirectMask|KeyPressMask;
    popup=
        XCreateWindow(xstf->theDisplay,RootWindow(xstf->theDisplay,xstf->theScreen),
                      (DisplayWidth(xstf->theDisplay,xstf->theScreen)-width)/2,
                      (DisplayHeight(xstf->theDisplay,xstf->theScreen)-height)/2,
                      width,height,10,DefaultDepth(xstf->theDisplay,xstf->theScreen),
                      InputOutput,xstf->theVisual,
                      CWBackPixel|CWOverrideRedirect|CWEventMask|CWBorderPixel,
                      &att);
    font_info=NULL;
    if((userfont=XGetDefault(xstf->theDisplay,"qvs","font"))!=NULL)
        font_info = XLoadQueryFont(xstf->theDisplay,userfont);
    if(font_info==NULL)
        if ((font_info = XLoadQueryFont(xstf->theDisplay,fontname1)) == NULL)
            if ((font_info = XLoadQueryFont(xstf->theDisplay,fontname2)) == NULL)
                if ((font_info = XLoadQueryFont(xstf->theDisplay,fontname3)) == NULL) {
                    fprintf(stderr,
                            "Couldn't load any of the following fonts %s %s %s\n",
                            fontname1, fontname2, fontname3);
                    exit(-2);
                }
    gc = XCreateGC(xstf->theDisplay,popup,0,gcval);
    XSetFont(xstf->theDisplay,gc,font_info->fid);
    XSetForeground(xstf->theDisplay,gc,
                   WhitePixel(xstf->theDisplay,xstf->theScreen));
    XSetBackground(xstf->theDisplay,gc,
                   BlackPixel(xstf->theDisplay,xstf->theScreen));
    XMapWindow(xstf->theDisplay,popup);
    if(XGrabKeyboard(xstf->theDisplay,popup,False,GrabModeAsync,GrabModeAsync,
                     CurrentTime) != GrabSuccess) {
        fprintf(stderr,"Couldn't grab keyboard");
        exit(-88);
    }

    {   i=100;
        while((XGrabPointer(xstf->theDisplay,popup,False,0,
                            GrabModeAsync,
                            GrabModeAsync,popup,None,CurrentTime)
                != GrabSuccess)
                && i!=0)    i--;
        if(i==0)  {
            fprintf(stderr,"Couldn't grab pointer\n");
        }

        XUngrabPointer(xstf->theDisplay,CurrentTime);
    }
    sprintf(echo,"%s (%s):",text,ans);
    XDrawImageString(xstf->theDisplay,popup,gc,10,50,echo,strlen(echo));
    while(1) {
        XEvent e;
        char key[10];
        KeySym keysym;
        int length;
        for(i=0; i<sizeof(key); i++) key[i]='\0';
        XNextEvent(xstf->theDisplay,&e);
        switch(e.type) {
        case KeyPress:
            if(XLookupString((XKeyEvent *)&e,key,sizeof(key),&keysym,NULL)) {
                if ((keysym==XK_Return)||(keysym==XK_KP_Enter)||(keysym==XK_Linefeed)||
                        (keysym==XK_space)) {
                    if((length=strlen(ansstring))!=0) strcpy(ans,ansstring);
                    XDestroyWindow(xstf->theDisplay,popup);
                    return;
                }
                else if (keysym==XK_Escape) {
                    XDestroyWindow(xstf->theDisplay,popup);
                    return;
                }
                else if ((keysym==XK_BackSpace)||(keysym==XK_Delete)) {
                    if((length=strlen(ansstring))>0) {
                        ansstring[length-1]='\0';
                        XClearWindow(xstf->theDisplay,popup);
                    }
                    else XBell(xstf->theDisplay,10);
                }
                else strcat(ansstring,key);
                sprintf(echo,"%s (%s):%s",text,ans,ansstring);
                XDrawImageString(xstf->theDisplay,popup,gc,10,50,echo,strlen(echo));
            }
        }
    }
    pdebug(5,"getstring() done\n");
}

void alarm_handler(int errnum)
{
    static XEvent exp_event;
    static int firsttime=1;
    static int ntimes=0;

    pdebug(5,"alarm_handler()\n");
    ntimes++;
    if(firsttime) {
        exp_event.type=Expose;
        exp_event.xexpose.type=Expose;
        exp_event.xexpose.send_event=0;
        exp_event.xexpose.display=x_info_ptr->theDisplay;
        exp_event.xexpose.window=x_info_ptr->theWindow;
        exp_event.xexpose.x=0;
        exp_event.xexpose.y=0;
        exp_event.xexpose.width=handlerwidth;
        exp_event.xexpose.height=handlerheight;
        exp_event.xexpose.count=0;
        firsttime=0;
    }
    //  increment_pic++;
    (*picno)+=increment;
    //  if((*picno)==maxpicno) (*picno)=0;
    //  (*picno)+=increment;
    while((*picno)>=maxpicno) (*picno)-=maxpicno;
    while((*picno)<0) (*picno)+=maxpicno;
    //  if(increment_pic==1) {
    map_pic_to_image(x_info_ptr,vndpicptr,*picno,lookuptable,*mw,*mh);
    draw_image_to_screen(x_info_ptr,*mw,*mh,vndpicptr[*picno].name);
    //    XSendEvent(x_info_ptr->theDisplay,x_info_ptr->theWindow,False,0,&exp_event);
    XFlush(x_info_ptr->theDisplay);
    //  }
    pdebug(5,"alarm_handler() done %d\n",ntimes);
    return;
}

int find_roi(int x, int y, Roi *rois, int nrois, int *currentroi)
{
    int i;

    pdebug(5,"find_roi()\n");
    for(i=0; i<nrois; i++) {
        if((rois[i].xmin<x)&&(rois[i].xmax>x)&&(rois[i].ymin<y)&&(rois[i].ymax>y)) {
            (*currentroi)=i;
        }
    }
    pdebug(5,"find_roi() done\n");
    return(*currentroi);
}


int save_rois(Roi *rois, char *roifile, int nrois)
{
    int i;
    FILE *file;
    char roiline[200];

    pdebug(5,"save_rois()\n");
    if(nrois) {
        if((file=fopen(roifile,"w"))==NULL) {
            ajerror("%s: couldn't open %s for writing",progname,file);
        }
        for(i=0; i<nrois; i++) {
            if(rois[i].inuse) {
                sprintf(roiline,"%d,%d,%d-%d,%d,%d\n",rois[i].xmin,rois[i].ymin,
                        rois[i].zmin,rois[i].xmax,rois[i].ymax,rois[i].zmax);
                fprintf(file,"%s",roiline);
            }
        }
        fclose(file);
    }
    pdebug(5,"save_rois() done\n");
    return(0);
}

int read_rois(Roi *rois, char *roifile, int *nrois)
{
    int i;
    char roitext[200];
    FILE *file;

    pdebug(5,"read_rois()\n");
    (*nrois)=0;
    if((file=fopen(roifile,"r")) != NULL) {
        for(i=0; i<MAXNUMROIS; i++) {
            if(fgets(roitext,199,file)==NULL) break;
            sscanf(roitext,"%d,%d,%d-%d,%d,%d",&(rois[i].xmin),&(rois[i].ymin),
                   &(rois[i].zmin),&(rois[i].xmax),&(rois[i].ymax),&(rois[i].zmax));
            rois[i].inuse=1;
            (*nrois)++;
        }
    }
    pdebug(5,"read_rois() done\n");
    return(*nrois);
}

int getdispls(vndPic *vndpics, int npics, int x, int y, double *xdispl,
              double *ydispl, double *zdispl, int *roiregion)
{
    double *dx, *dy, *dz;
    int i;
    char *startpos;
    int found;
    int region,x1,y1,z1,x2,y2,z2;
    int lastpos;
    int nx,ny;
    int ii,jj;
    int only2d=0;
    double xave=0, yave=0, zave=0;

 
 //int region,x1,y1,z1,x2,y2,z2;
//int roiregion[7];

    pdebug(5,"getdispls()\n");

    if(!(dx=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dx\n");
    }
    if(!(dy=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dy\n");
    }
    if(!(dz=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dz\n");
    }
    /* Figure out nx and ny, so we can jump to roughly the right place */
    startpos = vndpics[0].notes;
    while((startpos=strchr(startpos,'r'))!=NULL) {
        if(!strncmp(startpos,"region ",7)) {
            if(sscanf(startpos,
                      "region %d (%d,%d,%d-%d,%d,%d): dx=%lg, dy=%lg, dz=%lg\n",
                      &region,&x1,&y1,&z1,&x2,&y2,&z2,dx,dy,dz)!=10) {
                pdebug(10,"2D data\n");
                only2d=1;
                sscanf(startpos,
                       "region %d (%d,%d-%d,%d): dx=%lg, dy=%lg\n",
                       &region,&x1,&y1,&x2,&y2,dx,dy);
                z1=z2=*dz=0;
            }
            break;
        }
        startpos++;
    }
    nx = (int)(vndpics[0].width*widthratio)/(x2-x1)-3;
    ny = (int)(vndpics[0].height*heightratio)/(y2-y1)-3;
    ii = x/(x2-x1)-2;
    jj = y/(y2-y1)-2;
    if(ii<0) ii=0;
    if(jj<0) jj=0;
    if(ii>nx-1) ii=nx-1;
    if(jj>ny-1) jj=ny-1;
    lastpos = ((ii*(ny-1))+jj)*strlen(vndpics[0].notes)/(nx*ny);
    for(i=0; i<npics; i++) {
        pdebug(10,"Searching pic %d\n",i);
        found=0;
        startpos = vndpics[i].notes + lastpos - 1000;
        if(startpos < vndpics[i].notes) startpos = vndpics[i].notes;
        while((startpos=strchr(startpos,'r'))!=NULL) {
            if(!strncmp(startpos,"region ",7)) {
                if(only2d) {
                    sscanf(startpos,
                           "region %d (%d,%d-%d,%d): dx=%lg, dy=%lg\n",
                           &region,&x1,&y1,&x2,&y2,dx+i,dy+i);
                    z1=z2=*(dz+i)=0;
                } else {
                    sscanf(startpos,
                           "region %d (%d,%d,%d-%d,%d,%d): dx=%lg, dy=%lg, dz=%lg\n",region
                          &region,&x1,&y1,&z1,&x2,&y2,&z2,dx+i,dy+i,dz+i);
                }
                if((x1<=x)&&(x2>=x)&&(y1<=y)&&(y2>=y)) {
                    pdebug(10,"Found one!\n");
                    found=1;
                    break;
                }
                if((startpos=strchr(startpos,'\n'))==NULL) return 0;
            }
            startpos++;
        }
        xave+=dx[i];
        yave+=dy[i];
        zave+=dz[i];
        if(found==0) return 0;
    }
    xave/=npics;
    yave/=npics;
    zave/=npics;
    for(i=0; i<npics; i++) {
        //    xdispl[i]=dx[i]-xave;
        //    ydispl[i]=dy[i]-yave;
        //    zdispl[i]=dz[i]-zave;
        xdispl[i]=dx[i];
        ydispl[i]=dy[i];
        zdispl[i]=dz[i];
    }

    roiregion[0]=region;
    roiregion[1]=x1;
    roiregion[2]=y1;
    roiregion[3]=z1;
    roiregion[4]=x2;
    roiregion[5]=y2;
    roiregion[6]=z2;

    pdebug(5,"getdispls() done\n");


    return 1;
}

int getmotions(vndPic *vndpics, int npics, int x, int y, double *xmag,
               double *xphase, double *ymag, double *yphase, double *zmag,
               double *zphase, int *roiregion)
{
    double *dx, *dy, *dz;
    int i;
    char *startpos;
    int found;
    int region,x1,y1,z1,x2,y2,z2;
    int lastpos;
    int nx,ny;
    int ii,jj;
    int only2d=0;
    double xave=0, yave=0, zave=0;
    double xreal=0, ximag=0, yreal=0, yimag=0, zreal=0, zimag=0;
    double cosf, sinf;

    pdebug(5,"getmotions()\n");

    if(!(dx=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dx\n");
    }
    if(!(dy=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dy\n");
    }
    if(!(dz=(double *)malloc(npics*sizeof(double)))) {
        ajerror("couldn't malloc dz\n");
    }
    /* Figure out nx and ny, so we can jump to roughly the right place */
    startpos = vndpics[0].notes;
    
    while((startpos=strchr(startpos,'r'))!=NULL) {
        if(!strncmp(startpos,"region ",7)) {
            if(sscanf(startpos,
                      "region %d (%d,%d,%d-%d,%d,%d): dx=%lg, dy=%lg, dz=%lg\n",
                      &region,&x1,&y1,&z1,&x2,&y2,&z2,dx,dy,dz)!=10) {
                pdebug(10,"2D data\n");
                only2d=1;
                sscanf(startpos,"region %d (%d,%d-%d,%d): dx=%lg, dy=%lg\n",
                       &region,&x1,&y1,&x2,&y2,dx,dy);
                z1=z2=(*dz)=0;
            }
            break;
        }
        startpos++;
    }
    
    nx = (int)(vndpics[0].width*widthratio)/(x2-x1)-3;
    ny = (int)(vndpics[0].height*heightratio)/(y2-y1)-3;
    pdebug(5,"nx=%d, ny=%d\n",nx,ny);
    ii = x/(x2-x1)-2;
    jj = y/(y2-y1)-2;
    if(ii<0) ii=0;
    if(jj<0) jj=0;
    if(ii>nx-1) ii=nx-1;
    if(jj>ny-1) jj=ny-1;
    lastpos = ((ii*(ny-1))+jj)*strlen(vndpics[0].notes)/(nx*ny);
    for(i=0; i<npics; i++) {
        pdebug(10,"Searching pic %d\n",i);
        found=0;
        startpos = vndpics[i].notes + lastpos - 1000;
        if(startpos < vndpics[i].notes) startpos = vndpics[i].notes;
        while((startpos=strchr(startpos,'r'))!=NULL) {
            if(!strncmp(startpos,"region ",7)) {
                if(only2d) {
                    sscanf(startpos,"region %d (%d,%d-%d,%d): dx=%lg, dy=%lg\n",
                           &region,&x1,&y1,&x2,&y2,dx+i,dy+i);
                    z1=z2=(*(dz+i))=0;
                } else {
                    sscanf(startpos,
                           "region %d (%d,%d,%d-%d,%d,%d): dx=%lg, dy=%lg, dz=%lg\n",
                           &region,&x1,&y1,&z1,&x2,&y2,&z2,dx+i,dy+i,dz+i);
                }
                if((x1<=x)&&(x2>=x)&&(y1<=y)&&(y2>=y)) {
                    pdebug(10,"Found one!\n");
                    found=1;
                    break;
                }
                if((startpos=strchr(startpos,'\n'))==NULL) return 0;
            }
            startpos++;
        }
        xave+=dx[i];
        yave+=dy[i];
        zave+=dz[i];
        if(found==0) return 0;
    }
    xave/=npics;
    yave/=npics;
    zave/=npics;
    for(i=0; i<npics; i++) {
        dx[i]-=xave;
        dy[i]-=yave;
        dz[i]-=zave;
        cosf = cos(2*PI*i/npics);
        sinf = sin(2*PI*i/npics);
        xreal += dx[i]*cosf;
        ximag -= dx[i]*sinf;
        yreal += dy[i]*cosf;
        yimag -= dy[i]*sinf;
        zreal += dz[i]*cosf;
        zimag -= dz[i]*sinf;
    }
    xreal /= (npics/2);
    yreal /= (npics/2);
    zreal /= (npics/2);
    ximag /= (npics/2);
    yimag /= (npics/2);
    zimag /= (npics/2);
    *xmag = sqrt(xreal*xreal+ximag*ximag);
    *xphase = 180/PI*atan2(ximag,xreal);
    *ymag = sqrt(yreal*yreal+yimag*yimag);
    *yphase = 180/PI*atan2(yimag,yreal);
    *zmag = sqrt(zreal*zreal+zimag*zimag);
    *zphase = 180/PI*atan2(zimag,zreal);
    roiregion[0]=region;
    roiregion[1]=x1;
    roiregion[2]=y1;
    roiregion[3]=z1;
    roiregion[4]=x2;
    roiregion[5]=y2;
    roiregion[6]=z2;

    pdebug(5,"getmotions() done\n");
    return 1;
}

int handle_xevents(X11Stuff *xstf,int *picnum,unsigned char *lut,int *maxwidth,
                   int *maxheight,double *minmin,double *maxmax,double *gamma,
                   vndPic *vndpics,int npics, int *autoscale,
                   Roi *rois, int *nrois, int *currentroi, char *roifile)
{
    char text[1500];
    char junktext[20];
    int xx,yy,dxx,dyy;
    static int first_x=0, first_y=0, dfirst_x=0, dfirst_y=0;
    static int roi_x=0, roi_y=0;
    int gotkey;
    int i;
    char key[10];
    KeySym keysym;
    static struct itimerval tim,junk;
    static sigset_t myset, nothing;
    static struct sigaction myalarmhandler = {alarm_handler,0,0,NULL};
    static int firsttime=1;
    static int interval_us=70000;
    static int showing_movie=0;
    static double frame_rate;
    static int buttondown=0;
    static int shiftdown=0;
    static int ctrldown=0;
    double xmag, xphase, ymag, yphase, zmag, zphase;
    char tempfilename[PATH_MAX];
    int tempint;
    int roiregion[7];
    double *xdispl, *ydispl, *zdispl;
    double scaletemp;
    Cursor cur;

    pdebug(5,"handle_xevents()\n");

    if(firsttime) {
        frame_rate = 1000000/interval_us;
        x_info_ptr = xstf;
        handlerwidth = *maxwidth;
        handlerheight = *maxheight;
        vndpicptr=vndpics;
        picno=picnum;
        maxpicno=npics;
        lookuptable=lut;
        mw = maxwidth;
        mh = maxheight;
        if(!(xdispl=(double *)malloc(npics*sizeof(double)))) {
            ajerror("couldn't malloc xdispl!\n");
        }
        if(!(ydispl=(double *)malloc(npics*sizeof(double)))) {
            ajerror("couldn't malloc ydispl!\n");
        }
        if(!(zdispl=(double *)malloc(npics*sizeof(double)))) {
            ajerror("couldn't malloc zdispl!\n");
        }
        pdebug(5,"firsttime\n");
        if(sigemptyset(&nothing)) {
            ajerror("Error emptying sigemptyset");
        }
        if(sigemptyset(&myset)) {
            ajerror("Error emptying myset");
        }
        if(sigaddset(&myset,SIGALRM)) {
            ajerror("Error adding SIGALRM to myset");
        }
        sigaction(SIGALRM,&myalarmhandler,NULL);
        //    sigprocmask(SIG_BLOCK,&myset,NULL);
        junk.it_interval.tv_sec=junk.it_interval.tv_usec=0;
        junk.it_value.tv_sec=junk.it_value.tv_usec=0;
        tim.it_interval.tv_sec=tim.it_value.tv_sec=(interval_us/1000000);
        tim.it_interval.tv_usec=tim.it_value.tv_usec=(interval_us%1000000);
        firsttime=0;
    }
    XNextEvent(xstf->theDisplay,&(xstf->event));

   
    switch(xstf->event.type) {
    case Expose:
        pdebug(5,"Expose event\n");
        redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                      *picnum,*maxwidth,*maxheight,vndpics,rois,*nrois,*currentroi);
        //    draw_image_to_screen(xstf,*maxwidth,*maxheight,vndpics[*picnum].name);
        break;
    case MotionNotify:
        while(XCheckWindowEvent(xstf->theDisplay,xstf->theWindow,
                                Button1MotionMask|Button2MotionMask|Button3MotionMask,
                                &(xstf->event))==True);

        if(buttondown==1) {
            dxx=xstf->event.xmotion.x * widthratio;
            dyy= (*maxheight-1-xstf->event.xmotion.y) * heightratio;
            xx=xstf->event.xmotion.x ;
            yy= (*maxheight-1-xstf->event.xmotion.y);
            sprintf(text,"%3d,%3d  %6.3f  dx=%3d  dy=%3d  d %.1f    ",dxx,dyy,
                    vndpics[*picnum].pixvals[yy][xx],dxx-dfirst_x,dyy-dfirst_y,
                    sqrt((double)(dfirst_x-dxx)*(dfirst_x-dxx)+
                         (double)(dfirst_y-dyy)*(dfirst_y-dyy)));
            infowindow(xstf,CHANGETEXT,text,xstf->event.xbutton.x,
                       xstf->event.xbutton.y-50);
        } else if(buttondown==2) {
            if(ctrldown) {
                if(getdispls(vndpics,npics,(int)((xstf->event.xbutton.x)*widthratio),(int)((*maxheight-1-xstf->event.xbutton.y)*heightratio),xdispl,ydispl,zdispl,roiregion)) {
                    if((roiregion[3]==0)&&(roiregion[6]==0)) {
                        sprintf(text,"region %d: (%d,%d-%d,%d)\nx: ",
                                roiregion[0],roiregion[1],roiregion[2],
                                roiregion[4],roiregion[5]);
                    } else {
                        sprintf(text,"region %d: (%d,%d,%d-%d,%d,%d)\nx: ",
                                roiregion[0],roiregion[1],roiregion[2],
                                roiregion[3],roiregion[4],roiregion[5],roiregion[6]);
                    }
                    for(i=0; i<npics; i++) {
                        sprintf(junktext,"%5.3f ",xdispl[i]);
                        strcat(text,junktext);
                    }
                    sprintf(junktext,"   \ny: ");
                    strcat(text,junktext);
                    for(i=0; i<npics; i++) {
                        sprintf(junktext,"%5.3f ",ydispl[i]);
                        strcat(text,junktext);
                    }
                    if((roiregion[3])||(roiregion[6])) {
                        sprintf(junktext,"   \nz: ");
                        strcat(text,junktext);
                        for(i=0; i<npics; i++) {
                            sprintf(junktext,"%5.3f ",zdispl[i]);
                            strcat(text,junktext);
                        }
                        sprintf(junktext,"   ");
                        strcat(text,junktext);
                    }
                    //          sprintf(junktext,"\n");
                    //          strcat(text,junktext);
                    infowindow(xstf,CHANGETEXT,text,xstf->event.xbutton.x,
                               xstf->event.xbutton.y-50);
                }
            } else { 		
                if(getmotions(vndpics,npics,(int)((xstf->event.xbutton.x)*widthratio),(int)((*maxheight-1-xstf->event.xbutton.y)*heightratio),&xmag,&xphase,&ymag,&yphase,&zmag,&zphase,roiregion)) {
                    if((roiregion[3]==0)&&(roiregion[6]==0)) {
                        sprintf(text,"region %d: (%d,%d-%d,%d)\nx: %6.3f pixels at %8.3f deg\ny: %6.3f pixels at %8.3f deg\n",roiregion[0],roiregion[1],roiregion[2],
                                roiregion[4],roiregion[5],xmag,xphase,ymag,yphase);
                    } else {
                        sprintf(text,"region %d: (%d,%d,%d-%d,%d,%d)\nx: %6.3f pixels at %8.3f deg\ny: %6.3f pixels at %8.3f deg\nz: %6.3f pixels at %8.3f deg\n",roiregion[0],roiregion[1],roiregion[2],
                                roiregion[3],roiregion[4],roiregion[5],roiregion[6],
                                xmag,xphase,ymag,yphase,zmag,zphase);
                    }
                    infowindow(xstf,CHANGETEXT,text,xstf->event.xbutton.x,
                               xstf->event.xbutton.y-50);
                }
            }
        } else if (buttondown==3) {
            pdebug(5,"%d %d %d %d\n",roi_x,xstf->event.xmotion.x,
                   roi_y,xstf->event.xmotion.y);
            rois[*currentroi].xmin=MIN(roi_x,xstf->event.xmotion.x * widthratio) ;
            rois[*currentroi].xmax=MAX(roi_x,xstf->event.xmotion.x * widthratio) ;
            rois[*currentroi].ymax=(*maxheight -1-MIN(roi_y / heightratio,xstf->event.xmotion.y)) * heightratio ;
            rois[*currentroi].ymin=(*maxheight -1-MAX(roi_y / heightratio,xstf->event.xmotion.y)) * heightratio;
            pdebug(5,"rois[%d].inuse = %d, nrois = %d\n",
                   *currentroi,rois[*currentroi].inuse,*nrois);
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
        }
        break;
    case ButtonPress:
        pdebug(5,"ButtonPress event\n");
        if(xstf->event.xbutton.button==1) {
            buttondown=1;
            if(XGrabPointer(xstf->theDisplay,xstf->theWindow,False,
                            Button1MotionMask|ButtonReleaseMask|ButtonPressMask,
                            GrabModeAsync,GrabModeAsync,xstf->theWindow,None,
                            CurrentTime)!=GrabSuccess)
                ajerror("couldn't grab pointer!\n");
            cur = XCreateFontCursor(xstf->theDisplay,XC_crosshair);
            XDefineCursor(xstf->theDisplay,xstf->theWindow,cur);
            first_x=xstf->event.xbutton.x;
            first_y=*maxheight-1-xstf->event.xbutton.y;
            dfirst_x=first_x * widthratio;
            dfirst_y=first_y * heightratio;

            XSetForeground(xstf->theDisplay,
                           DefaultGC(xstf->theDisplay,xstf->theScreen),
                           xstf->pixels[YELLOW]);
            makecursor(xstf,first_x,first_y,*maxheight);
            sprintf(text,"%3d,%3d  %6.3f                                    ",
                    first_x,first_y,
                    vndpics[*picnum].pixvals[first_y][first_x]);
            infowindow(xstf,CREATEPOPUP,text,xstf->event.xbutton.x,
                       xstf->event.xbutton.y-50);
        } else if(xstf->event.xbutton.button==2) {
            if(ctrldown) {
                if(getdispls(vndpics,npics,(int)((xstf->event.xbutton.x)*widthratio),(int)((*maxheight-1-xstf->event.xbutton.y)*heightratio),xdispl,ydispl,zdispl,roiregion)) {
                    buttondown=2;
                    if(XGrabPointer(xstf->theDisplay,xstf->theWindow,False,
                                    Button2MotionMask|ButtonReleaseMask|ButtonPressMask,
                                    GrabModeAsync,GrabModeAsync,xstf->theWindow,None,
                                    CurrentTime)!=GrabSuccess)
                        ajerror("couldn't grab pointer!\n");
                    makesmallcursor(xstf,xstf->event.xbutton.x,
                                    *maxheight-1-xstf->event.xbutton.y,*maxheight);
                    if((roiregion[3]==0)&&(roiregion[6]==0)) {
                        sprintf(text,"region %d: (%d,%d-%d,%d)\nx: ",
                                roiregion[0],roiregion[1],roiregion[2],
                                roiregion[4],roiregion[5]);
                    } else {
                        sprintf(text,"region %d: (%d,%d,%d-%d,%d,%d)\nx: ",
                                roiregion[0],roiregion[1],roiregion[2],
                                roiregion[3],roiregion[4],roiregion[5],roiregion[6]);
                    }
                    for(i=0; i<npics; i++) {
                        sprintf(junktext,"%5.3f ",xdispl[i]);
                        strcat(text,junktext);
                    }
                    sprintf(junktext,"\ny: ");
                    strcat(text,junktext);
                    for(i=0; i<npics; i++) {
                        sprintf(junktext,"%5.3f ",ydispl[i]);
                        strcat(text,junktext);
                    }
                    if((roiregion[3])||(roiregion[6])) {
                        sprintf(junktext,"\nz: ");
                        strcat(text,junktext);
                        for(i=0; i<npics; i++) {
                            sprintf(junktext,"%5.3f ",zdispl[i]);
                            strcat(text,junktext);
                        }
                    }
                    //          sprintf(junktext,"\n");
                    //          strcat(text,junktext);
                    infowindow(xstf,CREATEPOPUP,text,xstf->event.xbutton.x,
                               xstf->event.xbutton.y-50);
                }
            } else {
		pdebug(5,"X: %d Y: %d\n",(int)((xstf->event.xbutton.x)*widthratio),(int)((*maxheight-1-xstf->event.xbutton.y)*heightratio));
                if(getmotions(vndpics,npics,(int)((xstf->event.xbutton.x)*widthratio),(int)((*maxheight-1-xstf->event.xbutton.y)*heightratio),&xmag,&xphase,&ymag,&yphase,&zmag,&zphase,roiregion)) {
                    buttondown=2;
                    if(XGrabPointer(xstf->theDisplay,xstf->theWindow,False,
                                    Button2MotionMask|ButtonReleaseMask|ButtonPressMask,
                                    GrabModeAsync,GrabModeAsync,xstf->theWindow,None,
                                    CurrentTime)!=GrabSuccess)
                        ajerror("couldn't grab pointer!\n");
                    makesmallcursor(xstf,xstf->event.xbutton.x,
                                    *maxheight-1-xstf->event.xbutton.y,*maxheight);
                    if((roiregion[3]==0)&&(roiregion[6]==0)) {
                        sprintf(text,"region %d: (%d,%d-%d,%d)\nx: %6.3f pixels at %8.3f deg\ny: %6.3f pixels at %8.3f deg\n",roiregion[0],roiregion[1],roiregion[2],
                                roiregion[4],roiregion[5],xmag,xphase,ymag,yphase);
                    } else {
                        sprintf(text,"region %d: (%d,%d,%d-%d,%d,%d)\nx: %6.3f pixels at %8.3f deg\ny: %6.3f pixels at %8.3f deg\nz: %6.3f pixels at %8.3f deg\n",roiregion[0],roiregion[1],roiregion[2],
                                roiregion[3],roiregion[4],roiregion[5],roiregion[6],
                                xmag,xphase,ymag,yphase,zmag,zphase);
                    }
                    infowindow(xstf,CREATEPOPUP,text,xstf->event.xbutton.x,
                               xstf->event.xbutton.y-50);
                }
            }
        } else if(xstf->event.xbutton.button==3) {
            buttondown=3;
            pdebug(5,"buttondown=3\n");
            if(XGrabPointer(xstf->theDisplay,xstf->theWindow,False,
                            Button3MotionMask|ButtonReleaseMask|ButtonPressMask,
                            GrabModeAsync,GrabModeAsync,xstf->theWindow,None,
                            CurrentTime)!=GrabSuccess)
                ajerror("couldn't grab pointer!\n");
            roi_x = xstf->event.xbutton.x * widthratio;
            roi_y = (xstf->event.xbutton.y) * heightratio;
            rois[*currentroi].inuse=1;
            if(*nrois==0) (*nrois)++;
        } else {
            buttondown=0;
        }
        break;
    case ButtonRelease:
        pdebug(5,"ButtonRelease event\n");
        if(xstf->event.xbutton.button==1) {
            buttondown=0;
            cur = XCreateFontCursor(xstf->theDisplay,XC_draft_small);
            XDefineCursor(xstf->theDisplay,xstf->theWindow,cur);
            XUngrabPointer(xstf->theDisplay,CurrentTime);
            infowindow(xstf,DESTROYPOPUP,text,0,0);
        } else if(xstf->event.xbutton.button==2) {
            if(buttondown==2) {
                buttondown=0;
                XUngrabPointer(xstf->theDisplay,CurrentTime);
                infowindow(xstf,DESTROYPOPUP,text,0,0);
            }
        } else if(xstf->event.xbutton.button==3) {
            buttondown=0;
            XUngrabPointer(xstf->theDisplay,CurrentTime);
        }
        break;
    case KeyPress:
        pdebug(5,"Keypress event\n");
        gotkey=XLookupString((XKeyEvent *)&(xstf->event),key,sizeof(key),
                             &keysym,NULL);
        if((keysym==XK_Shift_L)||(keysym==XK_Shift_R)) {
            pdebug(5,"shift down\n");
            shiftdown=1;
        } else if((keysym==XK_Control_L)||(keysym==XK_Control_R)) {
            pdebug(5,"ctrl down\n");
            ctrldown=1;
        } else if(((keysym==XK_Left)||(keysym==XK_KP_4))&&(rois[*currentroi].inuse)) {
            if(shiftdown) {
                pdebug(5,"left, shift down\n");
                rois[*currentroi].xmin++;
                if(rois[*currentroi].xmin>*maxwidth) rois[*currentroi].xmin=*maxwidth;
                if(rois[*currentroi].xmin > rois[*currentroi].xmax) {
                    tempint = rois[*currentroi].xmin;
                    rois[*currentroi].xmin = rois[*currentroi].xmax;
                    rois[*currentroi].xmax = tempint;
                }
            } else if(ctrldown) {
                pdebug(5,"left, ctrl down\n");
                rois[*currentroi].xmin--;
                rois[*currentroi].xmax--;
                if(rois[*currentroi].xmin < 0) rois[*currentroi].xmin=0;
                if(rois[*currentroi].xmax < 0) rois[*currentroi].xmax=0;
            } else {
                pdebug(5,"left, nothing down\n");
                rois[*currentroi].xmin--;
                if(rois[*currentroi].xmin < 0) rois[*currentroi].xmin=0;
            }
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
        } else if(((keysym==XK_Right)||(keysym==XK_KP_6))&&(rois[*currentroi].inuse)) {
            if(shiftdown) {
                rois[*currentroi].xmax--;
                if(rois[*currentroi].xmax<0) rois[*currentroi].xmax=0;
                if(rois[*currentroi].xmin > rois[*currentroi].xmax) {
                    tempint = rois[*currentroi].xmin;
                    rois[*currentroi].xmin = rois[*currentroi].xmax;
                    rois[*currentroi].xmax = tempint;
                }
            } else if(ctrldown) {
                rois[*currentroi].xmin++;
                rois[*currentroi].xmax++;
                if(rois[*currentroi].xmin > *maxwidth) rois[*currentroi].xmin=*maxwidth;
                if(rois[*currentroi].xmax > *maxwidth) rois[*currentroi].xmax=*maxwidth;
            } else {
                rois[*currentroi].xmax++;
                if(rois[*currentroi].xmax > *maxwidth) rois[*currentroi].xmax=*maxwidth;
            }
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
        } else if(((keysym==XK_Up)||(keysym==XK_KP_8))&&(rois[*currentroi].inuse)) {
            if(shiftdown) {
                rois[*currentroi].ymax--;
                if(rois[*currentroi].ymax<0) rois[*currentroi].ymax=0;
                if(rois[*currentroi].ymin > rois[*currentroi].ymax) {
                    tempint = rois[*currentroi].ymin;
                    rois[*currentroi].ymin = rois[*currentroi].ymax;
                    rois[*currentroi].ymax = tempint;
                }
            } else if(ctrldown) {
                rois[*currentroi].ymin++;
                rois[*currentroi].ymax++;
                if(rois[*currentroi].ymin>*maxheight) rois[*currentroi].ymin=*maxheight;
                if(rois[*currentroi].ymax>*maxheight) rois[*currentroi].ymax=*maxheight;
            } else {
                rois[*currentroi].ymax++;
                if(rois[*currentroi].ymax>*maxheight) rois[*currentroi].ymax=*maxheight;
            }
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
        } else if(((keysym==XK_Down)||(keysym==XK_KP_2))&&(rois[*currentroi].inuse)) {
            if(shiftdown) {
                rois[*currentroi].ymin++;
                if(rois[*currentroi].ymin>*maxheight) rois[*currentroi].ymin=*maxheight;
                if(rois[*currentroi].ymin > rois[*currentroi].ymax) {
                    tempint = rois[*currentroi].ymin;
                    rois[*currentroi].ymin = rois[*currentroi].ymax;
                    rois[*currentroi].ymax = tempint;
                }
            } else if(ctrldown) {
                rois[*currentroi].ymin--;
                rois[*currentroi].ymax--;
                if(rois[*currentroi].ymin < 0) rois[*currentroi].ymin=0;
                if(rois[*currentroi].ymax < 0) rois[*currentroi].ymax=0;
            } else {
                rois[*currentroi].ymin--;
                if(rois[*currentroi].ymin < 0) rois[*currentroi].ymin=0;
            }
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
        } else if((keysym==XK_KP_Add)&&(rois[*currentroi].inuse)) {
            if(*currentroi==(*nrois)-1) (*nrois)++;
            (*currentroi)++;
            pdebug(5,"current roi = %d\n",*currentroi);
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
            save_rois(rois,roifile,*nrois);
        } else if((keysym==XK_KP_Subtract)&&(*currentroi > 0)) {
            (*currentroi)--;
            pdebug(5,"current roi = %d\n",*currentroi);
            redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                          *picnum,*maxwidth,*maxheight,vndpics,rois,
                          *nrois,*currentroi);
            save_rois(rois,roifile,*nrois);
        }
        if(gotkey) {
            switch(key[0]) {
            case 'q': /* Quit */
                pdebug(5,"key q\n");
                save_rois(rois,roifile,*nrois);
                if(use_shm) killshmimages();
                exit(0);
                break;
            case 'N': /* Next pic, wrapping around to first pic */
                pdebug(5,"key N\n");
                if(npics!=1) {
                    if(*picnum>=npics-1) {
                        (*picnum)=-1;
                    }
                }
                /* Note: no break here */
            case 'n': /* Next pic, stopping on last pic */
                pdebug(5,"key n\n");
                if(showing_movie) {
                    tim.it_interval.tv_sec=tim.it_value.tv_sec=0;
                    tim.it_interval.tv_usec=tim.it_value.tv_usec=0;
                    setitimer(ITIMER_REAL,&tim,&junk);
                    showing_movie=0;
                };
                if(npics!=1) {
                    pdebug(10,"*picnum=%d\n",*picnum);
                    if((*picnum)<npics-1) {
                        (*picnum)++;
                        pdebug(10,"*picnum=%d\n",*picnum);
                        redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                                      *picnum,*maxwidth,*maxheight,vndpics,
                                      rois,*nrois,*currentroi);
                    }
                }
                break;
            case 'P': /* Previous pic, wrapping around to last pic */
                pdebug(5,"key P\n");
                if(npics!=1) {
                    if(*picnum<=0) {
                        (*picnum)=npics; /* We'll subtract one in the next section */
                    }
                }
                /* Note: no break here */
            case 'p': /* Previous pic, stopping on first pic */
                pdebug(5,"key p\n");
                if(showing_movie) {
                    tim.it_interval.tv_sec=tim.it_value.tv_sec=0;
                    tim.it_interval.tv_usec=tim.it_value.tv_usec=0;
                    setitimer(ITIMER_REAL,&tim,&junk);
                    showing_movie=0;
                };
                if(npics!=1) {
                    if(*picnum>0) {
                        (*picnum)--;
                        redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                                      *picnum,*maxwidth,*maxheight,vndpics,
                                      rois,*nrois,*currentroi);
                    }
                }
                break;
            case 'a':
                pdebug(5,"key a\n");
                (*autoscale)=1;
                for(i=0; i<npics; i++) {
                    vndpics[i].image_mapped=0;
                }
                //        compute_stats(vndpics,maxwidth,maxheight,minmin,maxmax,npics);
                *maxmax= origmaxmax;
                *minmin= origminmin;
                pixval_floor= 0.0 - *minmin;
                pixval_slope= (LUT_SIZE-2)/(*maxmax - *minmin);
                for(i=0; i<npics; i++) {
                    vndpics[i].image_mapped=0;
                }
                fill_lut(lut,*minmin,*maxmax,*gamma,xstf->ngrays);
                redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                              *picnum,*maxwidth,*maxheight,vndpics,
                              rois,*nrois,*currentroi);
                break;
            case 'r':
                pdebug(5,"key r\n");
                (*autoscale)=0;
                for(i=0; i<npics; i++) {
                    vndpics[i].image_mapped=0;
                }
                compute_stats(vndpics,maxwidth,maxheight,minmin,maxmax,-1.0, -1.0, npics);
                pixval_floor= 0.0 - *minmin;
                pixval_slope= (LUT_SIZE-2)/(*maxmax - *minmin);
                for(i=0; i<npics; i++) {
                    vndpics[i].image_mapped=0;
                }
                fill_lut(lut,*minmin,*maxmax,*gamma,xstf->ngrays);
                redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                              *picnum,*maxwidth,*maxheight,vndpics,
                              rois,*nrois,*currentroi);
                break;
            case 'g':
                pdebug(5,"key g\n");
                getvalue(xstf,"Enter gamma value",gamma);
                for(i=0; i<npics; i++) {
                    vndpics[i].image_mapped=0;
                }
                fill_lut(lut,*minmin,*maxmax,*gamma,xstf->ngrays);
                redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                              *picnum,*maxwidth,*maxheight,vndpics,
                              rois,*nrois,*currentroi);
                break;
            case 's':
                pdebug(5,"key s\n");
                (*autoscale)=0;
                fprintf(stderr,"present min value and range: %g, %g\n",
                        *minmin,*maxmax-*minmin);
                getvalue(xstf,"Enter min value",minmin);
                getvalue(xstf,"Enter max value",maxmax);
                if(*maxmax < *minmin) {
                    scaletemp= *maxmax;
                    *maxmax= *minmin;
                    *minmin= scaletemp;
                }
                pixval_floor= 0.0 - *minmin;
                pixval_slope= (LUT_SIZE-2)/(*maxmax - *minmin);
                for(i=0; i<npics; i++) {
                    vndpics[i].image_mapped=0;
                }
                fill_lut(lut,*minmin,*maxmax,*gamma,xstf->ngrays);
                redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                              *picnum,*maxwidth,*maxheight,vndpics,
                              rois,*nrois,*currentroi);
                break;
            case 't':
                if(showing_movie) {
                    tim.it_interval.tv_sec=tim.it_value.tv_sec=0;
                    tim.it_interval.tv_usec=tim.it_value.tv_usec=0;
                    setitimer(ITIMER_REAL,&tim,&junk);
                    showing_movie=0;
                };
                getvalue(xstf,"Enter frame rate (pics/sec)",&frame_rate);
                if(frame_rate>MAXFRAMERATE) frame_rate=MAXFRAMERATE;
                interval_us = (int)(1000000/frame_rate);
                tim.it_interval.tv_sec=tim.it_value.tv_sec=(interval_us/1000000);
                tim.it_interval.tv_usec=tim.it_value.tv_usec=(interval_us%1000000);
                break;
            case 'f':
            case 'F':
                showing_movie=1;
                increment=1;
                if((tim.it_interval.tv_sec==0)&&(tim.it_interval.tv_usec==0)) {
                    interval_us = (int)(1000000/frame_rate);
                    tim.it_interval.tv_sec=tim.it_value.tv_sec=(interval_us/1000000);
                    tim.it_interval.tv_usec=tim.it_value.tv_usec=(interval_us%1000000);
                }
                pdebug(5,"setting itimer to %d.%d\n",tim.it_interval.tv_sec,
                       tim.it_interval.tv_usec);
                setitimer(ITIMER_REAL,&tim,&junk);
                break;
            case 'b':
            case 'B':
                showing_movie=1;
                increment=-1;
                if((tim.it_interval.tv_sec==0)&&(tim.it_interval.tv_usec==0)) {
                    interval_us = (int)(1000000/frame_rate);
                    tim.it_interval.tv_sec=tim.it_value.tv_sec=(interval_us/1000000);
                    tim.it_interval.tv_usec=tim.it_value.tv_usec=(interval_us%1000000);
                }
                pdebug(5,"setting itimer to %d.%d\n",tim.it_interval.tv_sec,
                       tim.it_interval.tv_usec);
                setitimer(ITIMER_REAL,&tim,&junk);
                break;
            case 'e':
                save_rois(rois,roifile,*nrois);
                strcpy(tempfilename,roifile);
                getstring(xstf,"Roifile",tempfilename);
                *currentroi = 0;
                for(i=0; i<MAXNUMROIS; i++) rois[i].inuse=0;
                read_rois(rois,tempfilename,nrois);
                strcpy(roifile,tempfilename);
                redisplay_pic(xstf,*autoscale,lut,*minmin,*maxmax,*gamma,
                              *picnum,*maxwidth,*maxheight,vndpics,
                              rois,*nrois,*currentroi);
                break;
            case 'z':
                rois[*currentroi].zmin = *picnum;
                if(rois[*currentroi].zmin > rois[*currentroi].zmax) {
                    tempint = rois[*currentroi].zmin;
                    rois[*currentroi].zmin = rois[*currentroi].zmax;
                    rois[*currentroi].zmax = tempint;
                }
                break;
            case 'x':
                rois[*currentroi].zmax = *picnum;
                if(rois[*currentroi].zmin > rois[*currentroi].zmax) {
                    tempint = rois[*currentroi].zmin;
                    rois[*currentroi].zmin = rois[*currentroi].zmax;
                    rois[*currentroi].zmax = tempint;
                }
                break;
            default:
                if((keysym==XK_KP_Add)||(keysym==XK_KP_Subtract)||(keysym==XK_Return)||
                        (keysym==XK_KP_Enter)||(keysym==XK_Linefeed)||(keysym==XK_space));
                else fprintf(stderr,valid_keys);
                break;
            }
        }
        break;
    case KeyRelease:
        pdebug(5,"Keyrelease event\n");
        gotkey=XLookupString((XKeyEvent *)&(xstf->event),key,sizeof(key),
                             &keysym,NULL);
        if((keysym==XK_Shift_L)||(keysym==XK_Shift_R)) {
            shiftdown=0;
            break;
        } else if((keysym==XK_Control_L)||(keysym==XK_Control_R)) {
            ctrldown=0;
            break;
        }
    }
    pdebug(5,"handle_xevents() done\n");
    return(0);
}


