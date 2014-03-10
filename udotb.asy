size(10cm,10cm);

real paletteheight=6cm;

import graph;
import palette;
import utils;


// get data size
int nx=getint("nx");
int ny=getint("ny");
int nz=getint("nz");

real[][][] ux=readfile(nx,ny,nz,getstring("ux"));
real[][][] uy=readfile(nx,ny,nz,getstring("uy"));
real[][][] uz=readfile(nx,ny,nz,getstring("uz"));

real[][][] bx=readfile(nx,ny,nz,getstring("bx"));
real[][][] by=readfile(nx,ny,nz,getstring("by"));
real[][][] bz=readfile(nx,ny,nz,getstring("bz"));

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

real[][] ux2=cut2(ux,nx,ny,nz,c,idir);
real[][] uy2=cut2(uy,nx,ny,nz,c,idir);
real[][] uz2=cut2(uz,nx,ny,nz,c,idir);

real[][] bx2=cut2(bx,nx,ny,nz,c,idir);
real[][] by2=cut2(by,nx,ny,nz,c,idir);
real[][] bz2=cut2(bz,nx,ny,nz,c,idir);

int nx=ux2.length;
int ny=ux2[0].length;

real[][] udotb=new real[nx][ny];
for(int i=0; i < nx; ++i) {
  for(int j=0; j < ny; ++j) {
    udotb[i][j]=
      (ux2[i][j]*bx2[i][j] + uy2[i][j]*by2[i][j] + uz2[i][j]*bz2[i][j])/
      sqrt(
	   (ux2[i][j]*ux2[i][j] + uy2[i][j]*uy2[i][j] + uz2[i][j]*uz2[i][j])*
	   (bx2[i][j]*bx2[i][j] + by2[i][j]*by2[i][j] + bz2[i][j]*bz2[i][j])
	   +1e-8
	   );
  }
}

pen[] Palette=paraview_cooltowarm;
bounds range;
pair a=(0,0);
pair b=imagedims(cutdir);
range=image(udotb,Full,a,b,Palette);  // Full colour bar

// Draw shape and remove wall region from image:
string shape=getstring("boundary shape");
if(shape == "circle") clipellipse(b.x,b.y,currentpicture);
if(shape == "rectangle") cliprectangle(b.x,b.y,currentpicture);

// Add the palette bar:
picture bar;
string barlegend="";
palette(bar,barlegend,range,(0,0),(0.5cm,paletteheight),Right,Palette,
        PaletteTicks(ptick=linewidth(0.5*linewidth())));
add(bar.fit(),point(E),30E);
