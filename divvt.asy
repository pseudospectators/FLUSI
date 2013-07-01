import graph;
import utils;

size(200,150,IgnoreAspect);
scale(Linear,Linear);

// used to over-ride the normal legend
// usage:
// asy vt -u "runlegs=\"asdf\""

// Use usersetting() to get optional 
string runlegs;
usersetting();
bool myleg=((runlegs == "") ? false: true);
string[] legends=set_legends(runlegs);

string yvar=getstring("y variable: divu, divb");

int ypos=0;
if(yvar == "divu") ypos=1;
if(yvar == "divb") ypos=2;

string datafile="dvt";


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
    
    string legend=myleg ? legends[n] : texify(run);
    draw(graph(t,a[ypos]),Pen(n),legend);
  }
}

// Optionally draw different data from files as well:
draw_another(myleg,legends,n);

// Draw axes
yaxis(yvar,LeftRight,LeftTicks);
xaxis("time",BottomTop,LeftTicks);
  
attach(legend(),point(plain.E),20plain.E);


