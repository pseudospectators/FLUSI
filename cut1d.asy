size(200,150,IgnoreAspect);
import graph;
import palette;
import utils;

// used to over-ride the normal legend
// usage:
// asy cut -u "runlegs=\"asdf\""
// asy cut -u "quantity=\"asdf\""

string runlegs;
string quantity="";
usersetting();

// find the comma-separated strings to use in the legend
bool myleg=((runlegs== "") ? false: true);
bool flag=true;
int n=-1;
int lastpos=0;
string legends[];
if(myleg) {
  string runleg;
  while(flag) {
    ++n;
    int pos=find(runlegs,",",lastpos);
    if(lastpos == -1) {runleg=""; flag=false;}
    
    runleg=substr(runlegs,lastpos,pos-lastpos);

    lastpos=pos > 0 ? pos+1 : -1;
    if(flag) legends.push(runleg);
  }
}

// get data size
int nx=getint("nx");
int ny=getint("ny");
int nz=getint("nz");

// load file
string filename=getstring("filename");
real[][][] f=readfile(nx,ny,nz,filename);

// select direction of cut
string cutdir=getstring("cut direction: x,y,z");
if(cutdir != "x" && cutdir != "y" && cutdir != "z") {
  write("Please select x, y, or z.");
  exit();
}

// find the line through the data:
string cutinfo;
int c1=getint("height 1");
int c2=getint("height 2");
int idir;
real[] f1;
if(cutdir == "x") {
  f1=new real[nx];
  idir=3;
  cutinfo=", y="+(string)c1+", z="+(string)c2;
}
if(cutdir == "y") {
  f1=new real[ny];
  idir=2;
  cutinfo=", x="+(string)c1+", z="+(string)c2;
}
if(cutdir == "z") {
  f1=new real[nz];
  idir=1;
  cutinfo=", x="+(string)c1+", y="+(string)c2;
}
f1=cut1(f,nx,ny,nz,c1,c2,idir);

// draw the graph, axes, and attach the legend.
string legend=myleg ? legends[0] : texify(filename);
draw(graph(sequence(f1.length),f1),legend+cutinfo);

yaxis(quantity,LeftRight,LeftTicks);
xaxis("position",BottomTop,LeftTicks);
//attach(legend(),point(plain.E),20plain.E);

