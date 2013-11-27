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



string filenames=getstring("filenames");
bool flag=true;
int n=-1;
int lastpos;

while(flag) {
  ++n;
  string filename;
  int pos=find(filenames,",",lastpos);
  if(lastpos == -1) {filename=""; flag=false;}
  filename=substr(filenames,lastpos,pos-lastpos);
  
  if(flag) {
    write(filename);
    lastpos=pos > 0 ? pos+1 : -1;
    
    // get data size
    int nx=getint("nx");
    int ny=getint("ny");
    int nz=getint("nz");


    // load file
    real[][][] f=readfile(nx,ny,nz,filename);

    // select direction of cut
    string cutdir=getstring("cut direction: x,y,z");
    if(cutdir != "x" && cutdir != "y" && cutdir != "z") {
      write("Please select x, y, or z.");
      exit();
    }
    real L=getreal("length");
    real d;
    
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
      d=L/nx;
    }
    if(cutdir == "y") {
      f1=new real[ny];
      idir=2;
      cutinfo=", x="+(string)c1+", z="+(string)c2;
      d=L/ny;
    }
    if(cutdir == "z") {
      f1=new real[nz];
      idir=1;
      cutinfo=", x="+(string)c1+", y="+(string)c2;
      d=L/nz;
    }
    f1=cut1(f,nx,ny,nz,c1,c2,idir);

    // draw the graph, axes, and attach the legend.
    string legend=myleg ? legends[0] : texify(filename);
    real x[];
    for(int i=0; i < f1.length; ++i)
      x.push(i*d);
    //if(n == 1) f1=-f1;
    
    string leg=texify(filename);
    draw(graph(x,f1),Pen(n),leg+cutinfo);
  }
}
yaxis(quantity,LeftRight,LeftTicks);
xaxis("position",BottomTop,LeftTicks);

attach(legend(),point(plain.E),20plain.E);
