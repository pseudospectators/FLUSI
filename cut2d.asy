size(10cm,10cm);

real paletteheight=6cm;

// take 2D cuts from images produced by mhd3d.

// specify the legend name via command-line arguments:
// asy -f pdf img.asy -u "legend=\"$v$\""


import graph;
import palette;
import utils;

// get data size
int nx=getint("nx");
int ny=getint("ny");
int nz=getint("nz");

string legend="";
usersetting();

string name=getstring("filename");

// select direction of cross-section
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

real[][][] f=readfile(nx,ny,nz,name);

if(getstring("use mask") =="y") {
  string maskname=getstring("mask filename");
  real[][][] mask=readfile(nx,ny,nz,maskname);
  
  for(int i=0; i < nx; ++i) {
    for(int j=0; j < ny; ++j) {
      for(int k=0; k < nz; ++k) {
	if(mask[i][j][k] != 0.0)
	  f[i][j][k] = 0.0;
      }
    }
  }
}

real[][] f2=cut2(f,nx,ny,nz,c,idir);

// Get dimensions of image:
string sl1, sl2;
if(cutdir == "x") { sl1="yl"; sl2="zl";}
if(cutdir == "y") { sl1="xl"; sl2="zl";}
if(cutdir == "z") { sl1="xl"; sl2="yl";} // CHECKME
real l1=getreal(sl1);
real l2=getreal(sl2);

// Find bounds for palette:
real f2max=f2[0][0];
real f2min=f2[0][0];
for(int i=0; i < f2.length; ++i) {
  for(int j=0; j < f2[i].length; ++j) {
    if(f2[i][j] > f2max) f2max=f2[i][j];
    if(f2[i][j] < f2min) f2min=f2[i][j];
  }
}
real f2absmax=max(abs(f2max),abs(f2min));

// choose a palette
//pen[] Palette=BWRainbow();
pen[] Palette=paraview_cooltowarm;
// symmetric colour bar:
bounds range=image(f2,Range(-f2absmax,f2absmax),(0,0),(l1,l2),Palette);
// Full colour bar:
//bounds range=image(f2,Full,(0,0),(l1,l2),Palette);

// Draw shape and remove wall region from image:
string shape=getstring("boundary shape");
if(shape == "circle") {
  pair O=(l1/2,l2/2);
  real ay=getreal("ay");
  if(ay != 0.0) {
    path wall=ellipse(O,1,1/sqrt(ay));
    draw(wall);
    clip(wall);
  }
}
if(shape == "rectangle") {
  real w=getreal("width");
  real h=getreal("height");

  pair O=(l1/2,l2/2);
  pair pw=(0.5*w,0), ph=(0,0.5*h);
  path wall=O-pw-ph--O-pw+ph--O+pw+ph--O+pw-ph--cycle;
  draw(wall);
  clip(wall);
}

// Add the palette bar:
picture bar;
string barlegend="";
palette(bar,barlegend,range,(0,0),(0.5cm,paletteheight),Right,Palette,
        PaletteTicks(ptick=linewidth(0.5*linewidth())));
add(bar.fit(),point(E),30E);

label(legend+", "+"$i_"+cutdir+"$"+"="+string(c),point(N),N);
