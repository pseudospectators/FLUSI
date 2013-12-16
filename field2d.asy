size(20cm,20cm);

real paletteheight=6cm;

// Take 2D cuts from images datasets in binary format.

// To transform flusi .h5 files to binary, use hdf2binary (which just
// wraps h5dump).

// Specify the contour lines to be drawn:
// asy -f pdf img.asy -u "con=new real[] {-2,2};"

// Specify the legend name via command-line arguments:
// asy -f pdf img.asy -u "legend=\"$v$\""
// (currently deprecated).


import graph;
import palette;
import utils;

// get data size
int nx=getint("nx");
int ny=getint("ny");
int nz=getint("nz");

string legend="";
real[] con={};
string xfilename="";
string yfilename="";
usersetting();

string xname, yname;
if(xfilename == "") xname=getstring("xfilename");
else xname=xfilename;
if(yfilename == "") yname=getstring("yfilename");
else yname=yfilename;
// Select direction of cross-section
string cutdir=getstring("cut direction: x,y,z");
int idir;
if(cutdir != "x" && cutdir != "y" && cutdir != "z") {
  write("Please select x, y, or z.");
  exit();
} else {
  if(cutdir== "x") idir=3;
  if(cutdir== "y") idir=2;
  if(cutdir== "z") idir=1;
}
int c=getint(cutdir+"-height");

real[][][] fx=readfile(nx,ny,nz,xname);
real[][][] fy=readfile(nx,ny,nz,yname);

// Optionally set the field to zero if the mask is not set to zero.
if(getstring("use mask") =="y") {
  string maskname=getstring("mask filename");
  real[][][] mask=readfile(nx,ny,nz,maskname);
  for(int i=0; i < nx; ++i) {
    for(int j=0; j < ny; ++j) {
      for(int k=0; k < nz; ++k) {
	if(mask[i][j][k] != 0.0)
	  fx[i][j][k] = 0.0;
      }
    }
  }
}

// Take a 2D cut of the file
real[][] fx2=cut2(fx,nx,ny,nz,c,idir);
real[][] fy2=cut2(fy,nx,ny,nz,c,idir);

// Get dimensions of image:
string sl1, sl2;
if(cutdir == "x") { sl1="yl"; sl2="zl";}
if(cutdir == "y") { sl1="xl"; sl2="zl";}
if(cutdir == "z") { sl1="xl"; sl2="yl";} // CHECKME
real l1=getreal(sl1);
real l2=getreal(sl2);

real d1=l1/nx;
real d2=l2/ny;

for(int i=0; i < nx; i += 2) {
  for(int j=0; j < ny; j += 2) {
    real ax=fx2[i][j];
    real ay=fy2[i][j];
    arrow((i*d1,j*d2),(ax,ay),sqrt(ax*ax+ay*ay));
  }
}


// Draw shape and remove wall region from image:
string shape=getstring("boundary shape");
if(shape == "circle") {
  pair O=(l1/2,l2/2);
  real ay=getreal("ay");
  if(ay != 0.0) {
    path wall=ellipse(O,1,1/sqrt(ay));
    draw(wall,red);
    //    clip(wall);
  }
}
