import graph;
import utils;

size(200,150,IgnoreAspect);
//scale(Linear,Linear);


// used to over-ride the normal legend
// usage:
// asy vt -u "runlegs=\"asdf\""

// Use usersetting() to get optional 
string runlegs;
real tmax=realMax;
string thescale="linlog";
usersetting();
bool myleg=((runlegs == "") ? false: true);
string[] legends=set_legends(runlegs);

string yvar=getstring("y variable: dissu, dissb, total");


int ypos=0;
if(yvar == "dissu") ypos=1;
if(yvar == "dissb") ypos=2;
if(yvar == "total") ypos=-1;

string datafile="diss.t";


if(ypos == 0) {
  write("Invalid choice for y variable.");
  exit();
}

scale(Linear,Log);
if(thescale=="linlin")
  scale(Linear,Linear);

string runs=getstring("runs");
string run;
int n=-1;
bool flag=true;
int lastpos;
while(flag) {
  ++n;
  int pos=find(runs,",",lastpos);
  if(lastpos == -1) {run=""; flag=false;}
  run=substr(runs,lastpos,pos-lastpos);
  if(flag) {
    write(run);
    lastpos=pos > 0 ? pos+1 : -1;

    // load all of the data for the run:
    string filename, tempstring;
   
    filename=run+"/"+datafile;
    file fin=input(filename).line();
    real[][] a=fin.dimension(0,0);
    a=transpose(a);
    
    // get time:
    real[] t=a[0];
    pen p=Pentype(n);
    // if(n == 0) p+=longdashed;
    // if(n == 2) p=darkgreen+solid;
    
    string legend=myleg ? legends[n] : texify(run);
    if(ypos > 0) {
      draw(graph(t,a[ypos],t < tmax),p,legend);
    } else {
      draw(graph(t,a[1]+a[2], t < tmax),p,legend);
    }
  }
}

// Optionally draw different data from files as well:
draw_another(myleg,legends,n);

// Draw axes
if(yvar == "dissu") 
  yvar="$\nu\int\left|\omega\right|^2\mathrm{d}x$";
if(yvar == "dissb") 
  yvar="$\eta\int\left|\j\right|^2\mathrm{d}x$";
if(yvar == "total") 
  yvar="$\int \nu\left|\omega\right|^2 + \eta\left|\j\right|^2 \mathrm{d}x$";
yaxis(yvar,LeftRight,LeftTicks);
xaxis("time",BottomTop,LeftTicks);
  
//attach(legend(),point(plain.E),20plain.E);


