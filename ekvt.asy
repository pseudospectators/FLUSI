import graph;
import utils;

size(300,200,IgnoreAspect);
scale(Linear,Log);
//scale(Linear,Linear);

// used to over-ride the normal legend
// Usage:
// asy ekvt.asy -u "runlegs=\"asdf\""

// To turn on the linear interpolation, one must specify bounds, as in
// asy ekvt.asy -u "startstops=new real[] {60,130,70,120}"

// To only plot up to time tmax=30:
// asy ekvt.asy -u "tmax=30"

// Use usersetting() to get optional 
string runlegs;
real tmax=realMax;
real[] startstops;
usersetting();

//startstops= new real[] {60,130};

bool myleg=((runlegs == "") ? false: true);
string[] legends=set_legends(runlegs);

string yvar=getstring("y variable: Ekin,Ekinx,Ekiny,Ekinz");

int ypos=0;
if(yvar == "Ekin") ypos=1;
if(yvar == "Ekinx") ypos=2;
if(yvar == "Ekiny") ypos=3;
if(yvar == "Ekinz") ypos=4;

string datafile="ek.t";

pair linear(real[] x, real[] y, real start, real stop)
{
  pair AB=(0,0);

  real meanx=0.0, meany=0.0;
  int N=0;
  for(int i=0; i < x.length; ++i) {
    if(x[i] > start && x[i] < stop) {
      meanx += x[i];
      meany += y[i];
      ++N;
    }
  }
  if(N == 0) return AB;
  meanx /= N;
  meany /= N;

  real a=0, b=0;
  for(int i=0; i < x.length; ++i) {
    if(x[i] > start && x[i] < stop) {
      a += (x[i]-meanx)*(y[i]-meany);
      b += (x[i]-meanx)*(x[i]-meanx);
    }
  }
  
  return (meany-(a/b)*meanx,a/b);
 }

if(ypos == 0) {
  write("Invalid choice for y variable.");
  exit();
}

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

    // Restarting simulations from save-points earlier than the last
    // integral output produces overlapping data.  If the time series
    // is not monotonic, the code below discards later terms in favour
    // of earlier terms.
    // FIXME: this should be done when data is generated, not in
    // post-processing.
    real lastmax=-realMax;
    for(int i=0; i < t.length; ++i) {
      if(t[i] > lastmax)
	lastmax=t[i];
      else {
	int j=i;
	while(j < t.length && t[j] <= lastmax) {
	  t.delete(j);
	  a[ypos].delete(j);
	}
      }
    }
    
    string legend=myleg ? legends[n] : texify(run);
    draw(graph(t,a[ypos],t<tmax),Pen(n),legend);

    if(startstops.length > 0) {
      real start=startstops[2*n];
      real stop=startstops[2*n+1];
      
      real[] logEkin=log(a[ypos]);
      
      pair AB=linear(t,logEkin,start,stop);
      
      draw(graph(t,4*exp(AB.x)*exp(AB.y*t),
		 abs(t-0.5*(start+stop))<0.5*(stop-start)),
	   Pen(n)+dashed,"$\propto \exp("+string(AB.y,4)+"\,t)$");
    }
  }
}

// Optionally draw different data from files as well:
draw_another(myleg,legends,n);

// Draw axes
yaxis(yvar,LeftRight,LeftTicks);
xaxis("time",BottomTop,LeftTicks);
  
attach(legend(),point(plain.E),20plain.E);


