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

// choose a palette
pen[] Palette=BWRainbow();

//image(f[z],(0,0),(1,1),Palette); // just add the image:

real lx=2.513274;
real ly=2.513274;
bounds range=image(f2,(0,0),(lx,ly),Palette); // add the image

pair O=(lx/2,ly/2);
real ay=0;
ay=getreal("ay");
if(ay != 0.0) {
  path wall=ellipse(O,1,1/sqrt(ay));
  //dot(O);
  // path ellipse(pair c, real a, real b)
  draw(wall);
  clip(wall);
}


// add the palette bar:
picture bar;
string barlegend="";
palette(bar,barlegend,range,(0,0),(0.5cm,paletteheight),Right,Palette,
        PaletteTicks(ptick=linewidth(0.5*linewidth())));
add(bar.fit(),point(E),30E);


/*
pair center=(0.5+1/nx,0.5+1/ny);
real rad=0.5/(0.4*pi);
draw(circle(center,rad));
*/

label(legend+", "+"$i_"+cutdir+"$"+"="+string(c),point(N),N);
