: $Id: matrix.mod,v 1.64 2010/05/10 15:14:46 billl Exp $

:* COMMENT
COMMENT
NB: only minimal error checking done
NB: no dynamic allocation - eg m1.transpose(m1) will give a wrong result 
NB: matrix and vec sizes must be correct before using: use .resize()
NB: representation: first col values of vec are the first row

================            USAGE               ================
objref mat
mat = new Vector(rows*cols)
mat.mprintf(M,N)     // print out as M rows and N columns
mat2.transpose(mat,col)  // transpose of matrix
mat.revrows(col)  // reverse row order
y.mmult(mat,x)       // y = mat*x
y.spmult(pre,post,mat,x) // y = mat*x using "sparse matrix"
w.spget(pre,post,row,col) // ie pre,post,post,pre!!
spidget(pre,post,prid,poid ...) // match vals in 4 vecs to 4 vals
wt.mkspcp/chkspcp(pre,post) // copy the indices into integer arrays
mat.outprod(x,y)     // mat = outer product of vectors x and y
// All of the following can be replaced with vec.copy commands
mat.mget(i,j,cols)   // i=row#; j=col#
mat.mset(i,j,cols,val)
y.sector(mat,i,j,width,cols)
y.mrow(mat,i,cols)   // replace with y.resize(cols) y.copy(mat,0,j*cols,(j+1)*cols-1,1,1)
y.mcol(mat,j,cols)   // replace with y.resize(rows) y.copy(mat,0,j,rows*cols-1,1,cols)
y.msetrow(mat,i,cols)
y.msetcol(mat,j,cols)
================================================================
ENDCOMMENT

NEURON {
  SUFFIX nothing
  GLOBAL ROWS,COLS
}
 
PARAMETER {
  MATRIX_INSTALLED=0
}

ASSIGNED {
  ROWS
  COLS
}

VERBATIM
#include "misc.h"
ENDVERBATIM

:* ind.symmclean([COLS,flag]) flag to get rid of values from diagonal
VERBATIM
static double symmclean(void* vv) {
  int ix,dg,cols,i,j,n,flag,nx;
  double *x, *d;
  nx = vector_instance_px(vv, &x);
  cols =(ifarg(1)?((int)*getarg(1)):(int)COLS); 
  flag =(ifarg(2)?(int)*getarg(2):0);
  d=dcrset(nx);
  for (i=0;i<nx;i++) {
    if (x[i]!=floor(x[i]) || x[i]<0. || x[i]>cols*cols) {
      printf("symmclean ERR: OOB %g (%d x %d)\n",x[i],cols,cols); return -1; }
    ix=(int)x[i];
    dg=(n=ix/(cols+1))*(cols+1);
    if (ix<dg || ix>=dg+(cols-n) || (flag && ix==dg)) x[i]=-1; 
  }
  for (i=0,j=0;i<nx;i++) if (x[i]!=-1) d[j++]=x[i]; // keep these
  for (i=0;i<j;i++) x[i]=d[i]; // copy back
  vector_resize((IvocVect*)vv,j);
  return (double)j;
}
ENDVERBATIM

:* mat.outprod(x,y) // mat = outer product of vectors x and y
VERBATIM
static double outprod(void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z;
  /* this will be the outer product */
  nx = vector_instance_px(vv, &x);
	
  /* these are the two vectors that make it up */
  COLS=ny = vector_arg_px(1, &y); // will be number of columns
  ROWS=nz = vector_arg_px(2, &z); // will be number of rows
  if (nx != ny*nz) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<ny;i++) {
    for (j=0;j<nz;j++) {
      x[i*nz+j] = y[i]*z[j];
    }
  }
  return nx;
}
ENDVERBATIM
 
:* mmult
VERBATIM
static double mmult(void* vv) {
  int i, j, nx, ny, nz;
  double *x, *y, *z;
  /* x will be the product of matrix y and vec z */
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  nz = vector_arg_px(2, &z);
  if (ny != nx*nz) {
    hoc_execerror("Vector size mismatch", 0);
  }
  for (i=0;i<nx;i++) {
    x[i] = 0.;
    for (j=0;j<nz;j++) {
      x[i] += y[i*nz+j]*z[j];
    }
  }
  return nx;
}
ENDVERBATIM
 
:* ST[PO].spltp(pr,po,wt,ST[PRE])
VERBATIM
static double spltp(void* vv) {
  int ii, jj, nstpr, nstpo, nw, npr, npo, flag, cnt;
  double *stpr, *stpo, *w, *pr, *po;
  extern double hoc_call_func(Symbol*, int narg);

  char func[4] = "ltp";
  Symbol* s = hoc_lookup(func);
  if (! s) { hoc_execerror("Can't find ltp() func", 0); }
  nstpo = vector_instance_px(vv, &stpo);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nw = vector_arg_px(3, &w);
  nstpr = vector_arg_px(4, &stpr);
  for (ii=0,jj=0,cnt=0;ii<nstpo;ii++) {
    if (stpo[ii]==1.0) { /* connections to these will be changed */ 
      for (;po[jj]<ii;jj++) ; /* move forward till find a po */
      for (;po[jj]==ii;jj++) { /* move through these po's */
	if (stpr[(int)pr[jj]]==1.) { /*  did the presyn spike? */
	  cnt++; hoc_pushx(1.0);
	} else { 
	  cnt--; hoc_pushx(-1.0);
	}
        hoc_pushx(w[jj]);
        w[jj]=hoc_call_func(s, 2);
      }
    }
  }
  return cnt;
}
ENDVERBATIM
 
VERBATIM
/* Maintain a parallel vector of ints to avoid the slowness of repeated casts in spmult */
static int *pr_int;
static int *po_int;
static int cpfl=0;
ENDVERBATIM

:* wt.mkspcp(pr,po)
VERBATIM
static double mkspcp(void* vv) {
  int j, nw, npr, npo;
  double *w, *pr, *po;
  if (! ifarg(1)) { 
    cpfl=0; 
    if (po_int!=NULL) free(po_int); 
    if (pr_int!=NULL) free(pr_int);
    po_int=(int *)NULL; pr_int=(int *)NULL; 
    return 0;
  }
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  pr_int=(int *)ecalloc(nw, sizeof(int));
  po_int=(int *)ecalloc(nw, sizeof(int));
  for (j=0;j<nw;j++) {
    po_int[j]=(int)po[j];
    pr_int[j]=(int)pr[j];
  }
  cpfl=nw;
  return cpfl;
}
ENDVERBATIM

:* wt.chkspcp(pr,po)
VERBATIM
static double chkspcp(void* vv) {
  int j, nw, npr, npo, flag;
  double *w, *pr, *po;
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  flag=1;
  if (po_int==NULL || pr_int==NULL) { cpfl=0; return 0; }
  if (cpfl!=nw) { flag=0;
  } else for (j=0;j<nw;j++) {
    if (po_int[j]!=(int)po[j] || pr_int[j]!=(int)pr[j]) {flag=0; continue;}
  }
  if (flag==0) {
    cpfl=0; free(po_int); free(pr_int); 
    po_int=(int *)NULL; pr_int=(int *)NULL; 
  }
  return flag;
}
ENDVERBATIM

:* y.spmult(pr,po,wt,x[,flag])
: y=W*x, y will be the product of matrix w with pre/post indices and vec x
: optional flag (5th arg present) - do not clear dest vector initially
VERBATIM
static double spmult(void* vv) {
  int i, j, nx, ny, nw, npr, npo, flag;
  double *x, *y, *w, *pr, *po, xx;
  ny = vector_instance_px(vv, &y);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nw = vector_arg_px(3, &w);
  nx = vector_arg_px(4, &x);
  if (ifarg(5)) {flag=1;} else {flag=0;}
  if (nw!=npr || nw!=npo) {
    hoc_execerror("Sparse mat must have 3 identical size vecs for pre/post/wt", 0); 
  }
  if (flag==0) for (i=0;i<ny;i++) y[i] = 0.; // clear dest vec
  if (cpfl==0) {
    for (j=0;j<nw;j++) y[(int)po[j]] += (x[(int)pr[j]]*w[j]);
  } else if (cpfl!=nw) { hoc_execerror("cpfl!=nw in spmult", 0); } else {
    for (j=0;j<nw;j++) if (x[pr_int[j]]!=0) { y[po_int[j]] += ((x[pr_int[j]])*w[j]); }
  }
  return nx;
}
ENDVERBATIM
 

:* wt.spget(pr,po,row,col) returns weight value
VERBATIM
static double spget(void* vv) {
  int j, nw, npr, npo;
  double *w, *pr, *po, row, col;
  nw = vector_instance_px(vv, &w);
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  row = *getarg(3);
  col = *getarg(4);
  for (j=0;j<nw;j++) if (row==po[j]&&col==pr[j]) break;
  if (j==nw) return 0.; else return w[j];
}
ENDVERBATIM

:* spidget(prv,pov,pridv,poidv,pr,po,prid,poid) returns an index
FUNCTION spidget() { : gaussian distribution around 0
VERBATIM
{
  int j, npr, npo, nprid, npoid;
  double *pr, *po, *prid, *poid, pri, poi, pridi, poidi;
  npr = vector_arg_px(1, &pr);
  npo = vector_arg_px(2, &po);
  nprid = vector_arg_px(3, &prid);
  npoid = vector_arg_px(4, &poid);
  pri= *getarg(5);
  poi= *getarg(6);
  pridi= *getarg(7);
  poidi= *getarg(8);
  for (j=0;j<npr;j++) { 
    if (poi==po[j]&&pri==pr[j]&&pridi==prid[j]&&poidi==poid[j]) break;
  }
  if (j==npr) _lspidget=-1.0; else _lspidget=(double)j;
}
ENDVERBATIM
}

:* transpose
VERBATIM
static double transpose (void* vv) {
  int i, j, nx, ny, rows, cols, flag;
  double *x, *y;
  /* x will be the transpose of matrix y */
  nx = vector_instance_px(vv, &x);
  if (ifarg(2)) {
    ny = vector_arg_px(1, &y);
    COLS=cols=(ifarg(2))?(int)*getarg(2):(int)COLS;
    if (ny!=nx) hoc_execerror("Vector size mismatch", 0);
    flag=0;
  } else { 
    COLS=cols=(ifarg(1))?(int)*getarg(1):(int)COLS;
    y = (double *)ecalloc(nx, sizeof(double));
    flag=1;
  }
  ROWS=rows = nx/cols;
  if (flag) {
    for (i=0;i<rows;i++) for (j=0;j<cols;j++) y[j*rows+i] = x[i*cols+j];
    for (i=0;i<nx;i++) x[i] = y[i];
    free(y);
  } else {
    for (i=0;i<rows;i++) for (j=0;j<cols;j++) x[j*rows+i] = y[i*cols+j];
  }
  ROWS=cols; COLS=rows;
  return nx;
}
ENDVERBATIM
 
:* revrows() first row becomes last etc.
VERBATIM
static double revrows(void* vv) {
  int i, j, k, nx, rows, cols;
  double *x, tmp;
  nx = vector_instance_px(vv, &x);
  COLS=cols=(ifarg(1))?(int)*getarg(1):(int)COLS;
  ROWS=rows = nx/cols;
  for (i=0;i<rows/2;i++) {
    k=rows-i-1;
    for (j=0;j<cols;j++) {
      tmp=x[i*cols+j];
      x[i*cols+j]=x[k*cols+j];
      x[k*cols+j]=tmp;
    }
  }
  return nx;
}
ENDVERBATIM
 
:* mprintf() // rows cols stored in ROWS COLS
: mprintf(rows[,cols,format,prows,pcols]) // only print prows and pcols
: mprintf(rows[,cols,prows,pcols]) // only print up to prows and pcols
: mprintf(format[,prrows,procols])
: mprintf(.... rbeg,cbeg,prows,pcols) // only print this piece of matrix
VERBATIM
static double mprintf (void* vv) {
  int i, j, nx, rows, cols, rbeg, cbeg, prows, pcols, pfl; char *format;
  double *x;
  /* x will be printed out */
  nx = vector_instance_px(vv, &x);
  i=1; rbeg=cbeg=pfl=0; prows=rows=(int)ROWS; pcols=cols=(int)COLS;
  if (ifarg(i)) {
    if (!hoc_is_str_arg(i)) {ROWS=prows=rows=(int)*getarg(i++); COLS=pcols=cols=(int)*getarg(i++);}
    if (ifarg(i) && hoc_is_str_arg(i)) { pfl=1; format=gargstr(i++); }
    if (ifarg(i)) prows=(int)*getarg(i++);
    if (ifarg(i)) pcols=(int)*getarg(i++);
    if (ifarg(i)) {rbeg=prows; prows=(int)*getarg(i++);}
    if (ifarg(i)) {cbeg=pcols; pcols=(int)*getarg(i++);}
  }
  if (nx != rows*cols) { printf("Vector size mismatch"); hxe(); }
  if (prows>rows||pcols>cols){printf("Prsize mismatch: %d %d %d %d",prows,rows,pcols,cols);hxe();}
  for (i=rbeg;i<prows;i++) {
    for (j=cbeg;j<pcols;j++) {
      if (pfl) printf(format,x[i*cols+j]); else printf("%g\t",x[i*cols+j]);
    }
    printf("\n");
  }
  return nx;
}
ENDVERBATIM
 
:* vec.mshuffle(cols) shuffle within each row
VERBATIM
static double mshuffle (void* vv) {
  int i,j,k,n,nx,rows,cols; double *x, y[1], temp, augstep;
  nx=vector_instance_px(vv, &x);
  COLS=cols=(ifarg(1))?(int)*getarg(1):(int)COLS;
  ROWS=(double)nx/COLS;
  if (ROWS!=floor(ROWS)) {printf("matrix mshuffle ROW/COL ERR %g %g %d\n",ROWS,COLS,nx); hxe();}
  rows=ROWS;
  for (i=0;i<nx;i+=cols) dshuffle(x+i,cols);
  return (double)nx;
}
ENDVERBATIM

:* x.sector(mat,i,j,width,cols) // for a square
:* x.sector(mat,i,j,widx,widy,cols)
: pull out a piece of the matrix from i,y to i+widx,j+widy
VERBATIM
static double sector (void* vv) {
  int i, j, ii, jj, kk, nx, rows, cols, widx, widy, nm, err;
  double *x, *m;
  nx = vector_instance_px(vv, &x);
  nm = vector_arg_px(1, &m);
  i = (int)*getarg(2);
  j = (int)*getarg(3);
  widx = (int)*getarg(4);
  if (ifarg(6)) { widy = (int)*getarg(5); cols = (int)*getarg(6);
  } else {       widy=widx;               cols = (int)*getarg(5);  }
  err=0;
  ROWS=rows=nm/cols;
  if (nx != widx*widy) {printf("sector() ERRA: Vector size mismatch:%d %d\n",nx,widx*widy); hxe();}
  for (ii=0,kk=0;ii<widy && !err;ii++) for (jj=0;jj<widx && !err;jj++) {
    if ((i+ii)>rows || (j+jj)>cols) { 
      printf("WARN: fell off edge: %d %d %d %d\n",(i+ii),rows,(j+jj),cols); err=1; }
    x[kk++]=m[(i+ii)*cols+(j+jj)];
  }
  return x[i*cols+j];
}
ENDVERBATIM
 
:* vec.ppmrd(FILE)
VERBATIM 
static double ppmrd (void* vv) {
  int code, ii, type, num, maxsz, rows, cols, max, a,b,c;
  char tstr[256];
  double *x;
  FILE* f;

  num = vector_instance_px(vv, &x);
  maxsz=vector_buffer_size((IvocVect*)vv);
  f =     hoc_obj_file_arg(1);
  fgets(tstr,256,f);
  sscanf(tstr,"P%d",&type);
  fgets(tstr,256,f); // # CREATOR line; ? optional    
  fgets(tstr,256,f);
  sscanf(tstr,"%d %d",&cols,&rows);
  printf("Reading %d rows x %d cols\n",rows,cols);    
  fgets(tstr,256,f); // usually 255 for max
  sscanf(tstr,"%d",&max);
  if (maxsz<rows*cols) { printf("ppmrd vec too small %d %d.",maxsz,rows*cols); hxe(); }
  vector_resize((IvocVect*)vv, rows*cols);
  switch (type) {
    case 3:
    ii=0;
    while (fscanf(f,"%d %d %d",&a,&b,&c)==3) x[ii++]=(double)(a*65536+b*256+c);
    if (ii!=rows*cols) { printf("ppmrd only read %d of %d.",ii,rows*cols); hxe(); }
    break;
    case 6:
    for (ii=0;ii<rows*cols;ii++) {
      fread(&tstr,1,3,f);
      x[ii]=(double)tstr[0]*65536.+(double)tstr[1]*256.+(double)tstr[2];
    }
    break;
    default:      
    printf("ppmrd can't read type %d.\n",type); hxe();
    break;
  }
  return cols;
}
ENDVERBATIM

:* vec.ppmwr(FILE,cols)
VERBATIM 
static double ppmwr (void* vv) {	
  int code, ii, type, num, maxsz, rows, cols, max;
  unsigned char a[3], err;
  char tstr[256];
  double *x;
  FILE* f;

  num = vector_instance_px(vv, &x);
  f =     hoc_obj_file_arg(1);
  cols = (int)*getarg(2);
  if (ifarg(3)) type = (int)*getarg(3); else type=6;
  ROWS=rows=num/cols; err=0; 
  fprintf(f,"P%d\n",type);
  fprintf(f,"# CREATOR: NEURON ppmwr %s %s\n","$Date: 2010/05/10 15:14:46 $","$Revision: 1.64 $");
  fprintf(f,"%d %d\n",cols,rows);
  printf("Saving %d rows x %d cols\n",rows,cols);    
  max=255;
  fprintf(f,"%d\n",max);
  for (ii=0;ii<num && !err;ii++) {
    a[0]=(short)(x[ii]/65536); a[1]=(short)((x[ii]-(double)a[0]*65536.)/256.); 
    a[2]=(short)(x[ii]-(double)a[0]*65536.-(double)a[1]*256.); 
    switch (type) {
      case 6:
      fwrite(&a,1,3,f);
      break;
      case 3:
      fprintf(f,"%d %d %d ",a[0],a[1],a[2]);
      if (ii%6==0) fprintf(f,"\n");
      break;
      default:      
      printf("ppmrd can't read type %d.\n",type); hxe();
      err=1;
      break;
    }
  }
  fflush(f);
  return ii;
}
ENDVERBATIM

:* vec.rgb(R,G,B[,1])
: vec2rgb or rgb2vec, with flag goes rgb2vec
VERBATIM
static double rgb (void* vv) {
  int i, j, nx, nr, ng, nb, rows, cols, flag;
  double *x, *r, *g, *b;
  nx = vector_instance_px(vv, &x);
  nr = vector_arg_px(1, &r);
  ng = vector_arg_px(2, &g);
  nb = vector_arg_px(3, &b);
  if (nr!=nx || ng!=nx || nb!=nx) { printf("ERR rgb() size mismatch\n"); hxe(); }
  if (ifarg(4)) { // if flag move rgb to vec
    for (i=0;i<nx;i++) x[i]=r[i]*65536.+g[i]*256.+b[i];
  } else for (i=0;i<nx;i++) {
    r[i]=floor(x[i]/65536.); 
    g[i]=floor((x[i]-r[i]*65536.)/256.); 
    b[i]=x[i]-r[i]*65536.-g[i]*256.; 
  }
  return 0.0;
}
ENDVERBATIM

:* mget(i,j,cols)
VERBATIM
static double mget (void* vv) {
  int i, j, nx, rows, cols;
  double *x;
  nx = vector_instance_px(vv, &x);
  i = (int)*getarg(1);
  j = (int)*getarg(2);
  COLS=cols=(ifarg(3))?(int)*getarg(3):(int)COLS;
  if (i*cols+j >= nx) {
    hoc_execerror("Indices out of bounds", 0);
  }
  return x[i*cols+j];
}
ENDVERBATIM
 
:* mrow(mat,i,cols)
VERBATIM
static double mrow (void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  i = (int)*getarg(2);
  COLS=cols=(ifarg(3))?(int)*getarg(3):(int)COLS;
  ROWS=rows=ny/cols;
  if (i>=rows) { hoc_execerror("Indices out of bounds", 0); }
  if (cols!=nx) x=vector_newsize((IvocVect*)vv, nx=cols);
  for (j=0;j<nx;j++) { x[j] = y[i*cols+j]; }
  return nx;
}
ENDVERBATIM
 
:* mcol(mat,j,cols)
VERBATIM
static double mcol(void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  j = (int)*getarg(2);
  COLS=cols=(ifarg(3))?(int)*getarg(3):(int)COLS;
  if (j>=cols) { hoc_execerror("Indices out of bounds", 0); }
  ROWS=rows=ny/cols;
  if (rows!=nx) x=vector_newsize((IvocVect*)vv, nx=rows);
  for (i=0;i<nx;i++) { x[i] = y[i*cols+j]; }
  return nx;
}
ENDVERBATIM
 
:* msetrow(mat,i,cols)
VERBATIM
static double msetrow (void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  i = (int)*getarg(2);
  COLS=cols=(ifarg(3))?(int)*getarg(3):(int)COLS;
  if (cols!=nx || i>=ny/cols) {
    hoc_execerror("Indices out of bounds", 0);
  }
  for (j=0;j<nx;j++) { y[i*cols+j] = x[j]; }
  return nx;
}
ENDVERBATIM

:* msetcol(mat,j,cols)
VERBATIM
static double msetcol (void* vv) {
  int i, j, nx, ny, rows, cols;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  j = (int)*getarg(2);
  COLS=cols=(ifarg(3))?(int)*getarg(3):(int)COLS;
  if (cols!=ny/nx || j>=cols) {
    hoc_execerror("Indices out of bounds", 0);
  }
  for (i=0;i<nx;i++) { y[i*cols+j] = x[i]; }
  return nx;
}
ENDVERBATIM
 
:* mset (i,j,cols,val)
VERBATIM
static double mset(void* vv) {
  int i, j, nx, rows, cols;
  double *x, val;
  nx = vector_instance_px(vv, &x);
  i = (int)*getarg(1);
  j = (int)*getarg(2);
  COLS=cols=(ifarg(3))?(int)*getarg(3):(int)COLS;
  val=(ifarg(4))?*getarg(4):*getarg(3);
  if (i*cols+j >= nx) {
    hoc_execerror("Indices out of bounds", 0);
  }
  return (x[i*cols+j]=val);
}
ENDVERBATIM

VERBATIM 

typedef struct RUNT { // horizontal run
  int left,right,y,width;
} runt;

typedef struct RUNL { // run list
  runt* p;
  int sz;
  int bufsz;
  int* marked;
} runl;

typedef struct COMPT { // component
  int left,right,top,bot,pixels,nruns,w,h;
  int* prun;
} compt;

typedef struct COMPL { // component list
  compt* p;
  int sz;
  int bufsz;
  runl* pruns;
} compl_struct;

compl_struct* alloccompl (int bufsz, runl* pruns) {
  compl_struct* r;
  r = (compl_struct*)calloc(1,sizeof(compl_struct));
  r->sz=0;
  r->bufsz = bufsz;
  r->p = (compt*) calloc(bufsz,sizeof(compt));
  r->pruns = pruns;
  return r;
}

void addcomp (compl_struct* pc, int maxruns) {
  if(pc->sz + 1 >= pc->bufsz) {
    pc->bufsz *= 10;
    pc->p = (compt*) realloc(pc->p,pc->bufsz*sizeof(compt));
  }
  memset(&pc->p[pc->sz],0,sizeof(compt));
  pc->p[pc->sz].prun = (int*) calloc(maxruns,sizeof(int));
  pc->sz++;
}

void freecompl (compl_struct** cc) {
  compl_struct* c; int i;
  c=cc[0];
  for(i=0;i<c->sz;i++) free(c->p[i].prun);
  free(c->p);
  cc[0]=0x0;
}

void addruntocomp (compl_struct* pc, int cidx, int ridx ) {
  compt* c;
  c = &pc->p[cidx];
  c->prun[c->nruns] = ridx;
  if( pc->pruns->p[ridx].left < c->left ) c->left = pc->pruns->p[ridx].left;
  if( pc->pruns->p[ridx].right > c->right ) c->right = pc->pruns->p[ridx].right;
  if( pc->pruns->p[ridx].y < c->top ) c->top = pc->pruns->p[ridx].y;
  if( pc->pruns->p[ridx].y > c->bot ) c->bot = pc->pruns->p[ridx].y;
  c->w = c->right - c->left + 1;
  c->h = c->bot - c->top + 1;
  c->pixels += pc->pruns->p[ridx].width;
  c->nruns++;
}

runl* allocrunl (int bufsz) {
  runl* r;
  r = (runl*)calloc(1,sizeof(runl));
  r->sz=0;
  r->bufsz = bufsz;
  r->p = (runt*) calloc(bufsz,sizeof(runt));
  r->marked = (int*) calloc(bufsz,sizeof(int));
  return r;
}

void freerunl (runl** rr) {
  runl* r;
  r = rr[0];
  free(r->p);
  free(r);
  free(r->marked);
  rr[0]=0x0;
}

void addrun (runl* pr, int left, int right, int y) {
  int oldsz,i;
  oldsz = pr->bufsz;
  if(pr->sz + 1 >= oldsz) {
    pr->bufsz *= 10;
    pr->p = (runt*) realloc(pr->p,pr->bufsz*sizeof(runt));
    pr->marked = (int*) realloc(pr->marked,pr->bufsz*sizeof(int));
    for(i=oldsz;i<pr->bufsz;i++) pr->marked[i]=0;
  }
  pr->p[pr->sz].left=left;
  pr->p[pr->sz].right=right;
  pr->p[pr->sz].width=right-left+1;
  pr->p[pr->sz].y=y;
  pr->sz+=1;
}

runl* findruns (double* pin, int rows, int cols) {
  int i,x,y,right,left,bufsz;
  runl* pr;
  bufsz=1000;
  pr = allocrunl(bufsz);
  for(y=0;y<rows;y++) {
    for(x=0;x<cols;x++) {
      if(pin[y*cols+x]) {
        right=left=x;
        for(right=x;right<cols && pin[y*cols+right];right++);
        if(right==cols) right=cols-1;
        for(left=x;left>=0 && pin[y*cols+left];left--);
        if(left<0) left=0;
        addrun(pr,left,right,y);
      }
    }
  }
  return pr;
}

int vintrunts (runt* r1, runt* r2) { // return 1 if they overlap horizontally and are 1 pixel apart in vertical
  return !(r1->left>r2->right || r2->left>r1->right) && abs(r1->y-r2->y)==1;
}

void floodruns (compl_struct* pc, int ridx, int cidx) {
  int i;
  runl* pr;
  pr = pc->pruns;
  if(pr->marked[ridx]) return;
  addruntocomp(pc,cidx,ridx);  
  pr->marked[ridx]=1;
  for(i=ridx-1;i>=0;i--) {
    if( pr->p[ridx].y - pr->p[i].y > 1) break;
    if(vintrunts(&pr->p[i],&pr->p[ridx])) floodruns(pc,i,cidx);
  }
  for(i=ridx+1;i<pr->sz;i++) {
    if( pr->p[i].y - pr->p[ridx].y > 1) break;
    if(vintrunts(&pr->p[i],&pr->p[ridx])) floodruns(pc,i,cidx);
  }
}

compl_struct* findcomps (runl* pruns) {
  compl_struct* pc; int i,j,cidx; runt *r1,*r2;
  pc = alloccompl(1000, pruns);
  for(i=0;i<pruns->sz;i++) {
    if(pruns->marked[i]) continue;
    cidx = pruns->sz;
    addcomp(pc,pruns->sz);
    floodruns(pc,i,cidx);
  }
  return pc;
}

static double ccomps (void* vv) {
  int x,y,rows,cols,szin,szout;
  double *pin,*plab,val;
  runl* pr; compl_struct* pc;
  szin=vector_instance_px(vv,&pin);
  szout=vector_arg_px(1, &plab);
  COLS=cols=(int)*getarg(3);
  ROWS=rows=szin/cols;
  if(szin!=szout) {
    printf("ccomps ERRA: input must be same size as output (%d,%d)\n",szin,szout);
    return 0.0;
  }

  pr = findruns(pin,rows,cols);
  pc = findcomps(pr);

  

  freerunl(&pr);
  freecompl(&pc);

  return 1.0;
  
}
ENDVERBATIM
 
:* PROCEDURE install_matrix()
PROCEDURE install_matrix() {
  if (MATRIX_INSTALLED==0) { 
  MATRIX_INSTALLED=1
VERBATIM
  /* the list of additional methods */
  install_vector_method("outprod", outprod);
  install_vector_method("mmult", mmult);
  install_vector_method("spmult", spmult);
  install_vector_method("spget", spget);
  install_vector_method("mkspcp", mkspcp);
  install_vector_method("chkspcp", chkspcp);
  install_vector_method("spltp", spltp);
  install_vector_method("transpose", transpose);
  install_vector_method("revrows", revrows);
  install_vector_method("mprintf", mprintf);
  install_vector_method("mshuffle", mshuffle);
  install_vector_method("mget", mget);
  install_vector_method("mset", mset);
  install_vector_method("mrow", mrow);
  install_vector_method("mcol", mcol);
  install_vector_method("msetrow", msetrow);
  install_vector_method("msetcol", msetcol);
  install_vector_method("sector", sector);
  install_vector_method("ppmrd", ppmrd);
  install_vector_method("ppmwr", ppmwr);
  install_vector_method("rgb", rgb);
  install_vector_method("symmclean", symmclean);
  install_vector_method("ccomps",ccomps);
ENDVERBATIM
  } else { printf("%s\n","$Id: matrix.mod,v 1.64 2010/05/10 15:14:46 billl Exp $") }
}
