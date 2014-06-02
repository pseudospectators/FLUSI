// compute A \codt B /(||A|| ||B|| + eps)
size(10cm,10cm);

real paletteheight=6cm;

import graph;
import palette;
import utils;

// get data size
int nx=getint("nx");
int ny=getint("ny");
int nz=getint("nz");

string filenameAx=getstring("Ax");
string filenameAy=getstring("Ay");
string filenameAz=getstring("Az");
string filenameBx=getstring("Bx");
string filenameBy=getstring("By");
string filenameBz=getstring("Bz");

real[] con={};

usersetting();

real[][][] Ax=readfile(nx,ny,nz,filenameAx);
real[][][] Ay=readfile(nx,ny,nz,filenameAy);
real[][][] Az=readfile(nx,ny,nz,filenameAz);
real[][][] Bx=readfile(nx,ny,nz,filenameBx);
real[][][] By=readfile(nx,ny,nz,filenameBy);
real[][][] Bz=readfile(nx,ny,nz,filenameBz);

real[][][] f=new real[nx][ny][nz];
for(int i=0; i < nx; ++i) {
  for(int j=0; j < ny; ++j) {
    for(int k=0; k < nz; ++k) {
      real ax=Ax[i][j][k];
      real ay=Ay[i][j][k];
      real az=Az[i][j][k];
      real bx=Bx[i][j][k];
      real by=By[i][j][k];
      real bz=Bz[i][j][k];
      f[i][j][k]=(ax*bx+ay*by+az*bz)
	/sqrt((ax*ax+ay*ay+az*az)*(bx*bx+by*by+bz*bz) +1e-10);
    }
  }
}


pen[] Palette=paraview_cooltowarm;

bool drawshape=false;

if(drawshape) {

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

  // Take a 2D cut of the file
  real[][] f2=cut2(f,nx,ny,nz,c,idir);
  
  real[][][] mask;
  real[][] mask2;
  bool boundsmask=getstring("use mask for bounds") =="y";
  if(boundsmask) {
    string maskname=getstring("mask filename");
    mask=readfile(nx,ny,nz,maskname);
    mask2=cut2(mask,nx,ny,nz,c,idir);
  }
  
  // Find bounds for palette:
  real f2max=-realMax;
  real f2min=realMax;
  for(int i=0; i < f2.length; ++i) {
    for(int j=0; j < f2[i].length; ++j) {
      //    if(!boundsmask ||mask2[i][i]==0.0) {
      if(f2[i][j] > f2max) f2max=f2[i][j];
      if(f2[i][j] < f2min) f2min=f2[i][j];
      //    }
    }
  }

  write("range: "+string(f2min)+" "+string(f2max));
  real f2absmax=max(abs(f2max),abs(f2min));

  pair a=(0,0);
  pair b=imagedims(cutdir);
  
  // Draw image and specify colour bar:
  bounds range;
  bool symbar=getstring("symmetric colour bar (y/n)") =="y";
  if(symbar)
    range=image(f2,Range(-1,1),a,b,Palette);
  else
    range=image(f2,Full,a,b,Palette);  // Full colour bar


  
  import contour;
  if(con.length > 0) {
    //draw(contour(f2,a,b,con));
    Label[] Labels=sequence(new Label(int i) {
	return Label(con[i] != 0 ? (string) con[i] : "",
		     Relative(unitrand()),(0,0)
		     ,UnFill(1bp)
		     );
      },con.length);
    draw(Labels,contour(f2,a,b,con));
  }


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

} else {
  import graph;
  import stats;
  size(400,200,IgnoreAspect);

  
  bool boundsmask=getstring("use mask for bounds") =="y";
  real[][][] mask;
  //real[][] mask2;
  
  if(boundsmask) {
    string maskname=getstring("mask filename");
    mask=readfile(nx,ny,nz,maskname);
    //mask2=cut2(mask,nx,ny,nz,c,idir);
  }
  
  
  real[] a={};
  for(int i=0; i < nx; ++i) {
    for(int j=0; j < ny; ++j) {
      for(int k=0; k < nz; ++k) {
	if(!boundsmask || mask[i][j][k]==0.0) {
	  a.push(f[i][j][k]);
	}
      }
    }
  }
  // Calculate "optimal" number of bins a la Shimazaki and Shinomoto.
  //int N=10*bins(a);
  int N=100;

  real mina=-1;
  real maxa=1;
  real dx=(maxa-mina)/N;
  real[] freq=frequency(a,mina,maxa,N);
  freq /= dx*sum(freq); // normalize
  real[] x=sequence(N);
  x *= dx;
  x -= 1.0;

  string filenameout=getstring("Save data to filename");
  if(filenameout != "" ) {
    file fout=output(filenameout);
    for(int i=0; i < N; ++i ) {
      write(fout,x[i],freq[i]);
      write(fout,endl);
      //      write(x[i],freq[i]);
    }
    close(fout);
  }
  
  draw(graph(x,freq),red);
    
  write("using "+string(N)+" bins");
  //  histogram(a,mina,maxa,N,normalize=true,lightblue,low=0,bars=false);
  //histogram(a,min(a),max(a),N,normalize=true,low=0,bars=false);
  
  xaxis("$A\cdot{}B/AB$",BottomTop,LeftTicks);
  yaxis("count",LeftRight,RightTicks(trailingzero));
  //xequals(-1/sqrt(2));
}
