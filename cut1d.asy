size(200,150,IgnoreAspect);
import graph;
import palette;
import utils;

// used to over-ride the normal legend
// usage:
// asy cut -u "runlegs=\"asdf\""
// asy cut -u "quantity=\"asdf\""


string type="cut";
string xlabel, ylabel;
type=getstring("type: cut or spectrum");


if(type == "cut") {
  scale(Linear,Linear);
  if(xlabel == "") xlabel="position";
  if(ylabel == "") ylabel="value";
}
if(type == "spectrum") {
  scale(Log,Log);
  if(xlabel == "") xlabel="wavenumber";
  if(ylabel == "") ylabel="magnitude";
}

string runlegs;
string quantity="";
usersetting();

// find the comma-separated strings to use in the legend
bool myleg=((runlegs== "") ? false: true);
string[] legends=set_legends(runlegs);

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

    // Optionally set the field to zero if the mask is not set to zero.
    if(getstring("use mask") =="y") maskit(f,nx,ny,nz);
    
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
    int idir;
    string sh1, sh2;
    real[] f1;
    if(cutdir == "x") {
      f1=new real[nx];
      idir=1;
      //cutinfo=", y="+(string)c1+", z="+(string)c2;
      d=L/nx;
      sh2="y height";
      sh1="z height";

    }
    if(cutdir == "y") {
      f1=new real[ny];
      idir=2;
      //cutinfo=", x="+(string)c1+", z="+(string)c2;
      d=L/ny;
      sh2="z height";
      sh1="x height";
    }
    if(cutdir == "z") {
      f1=new real[nz];
      idir=3;
      //cutinfo=", x="+(string)c1+", y="+(string)c2;
      d=L/nz;
      sh2="x height";
      sh1="y height";
    }
    int c2=getint(sh2);
    int c1=getint(sh1);
    f1=cut1(f,nx,ny,nz,c1,c2,idir);

    // draw the graph, axes, and attach the legend.
    string legend=myleg ? legends[0] : texify(filename);
    real x[];
    //if(n == 1) f1=-f1;

    pen p=Pentype(n);
    if(n == 0) p+=longdashed;
    if(n == 2) p=darkgreen+solid;

    if(n == 0) p=blue+dashed;
    if(n == 1) p=darkgreen+solid;
    
    string legend=myleg ? legends[n] : texify(filename);
    if(type == "cut") {
      for(int i=0; i < f1.length; ++i) x.push(i*d);
      draw(graph(x,f1),p,legend+cutinfo);
    }
    
    if(type == "spectrum") {
      for(int i=0; i < f1.length; ++i)
	x.push(i);
      pair[] ff=new pair[f1.length];
      for(int i=0; i < ff.length; ++i) 	ff[i]=f1[i];
      pair[] fff=fft(ff);
      //      string leg=texify(filename);
      draw(graph(x,abs(fff),x > 0 & x<fff.length/3.0-1),p,legend+cutinfo);
     }
    
  }
}
yaxis(ylabel+quantity,LeftRight,LeftTicks);
xaxis(xlabel,BottomTop,LeftTicks);

//attach(legend(),point(plain.E),20plain.E);
