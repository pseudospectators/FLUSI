import graph;

// Find the comma-separated strings to use in the legend
string[] set_legends(string runlegs)
{
  string[] legends;
  bool myleg=((runlegs== "") ? false: true);
  write(myleg);
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
  return legends;
}

// optionally draw different data from files as well:
void draw_another(bool myleg,string[] legends, int n)
{
  bool anotherlegend=true;
  bool anotherone=(getstring("compare with other file? y/n") == "y");
  while(anotherone) {
    string filename;
    filename=getstring("filename:");
    file fin=input(filename).line();
    real[][] a=fin.dimension(0,0);
    a=transpose(a);
    int n0=getint("column x");
    int n1=getint("column y");
    string legend=myleg ? legends[n] : texify(filename);
    if(anotherlegend) {
      draw(graph(a[n0],a[n1]),Pen(n),legend);
      anotherlegend=false;
    } else {
      draw(graph(a[n0],a[n1]),Pen(n));
    }
    anotherone=(getstring("another one? y/n") == "y");
  }
}

