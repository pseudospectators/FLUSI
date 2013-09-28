import graph;

// Find the comma-separated strings to use in the legend
string[] set_legends(string runlegs)
{
  string[] legends;
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


// Load file (note format change from FORTRAN to C)
real[][][] readfile(int nx, int ny, int nz, string name) {
  file fin=input(name,mode="binary").singlereal();
  //file fin=input(name,mode="binary");
  //file fin=input(name);
  real[][][] f=new real[nx][ny][nz];

  // There is some really weird ordering in .h5 files.
  for(int k=0; k < nz; ++k) {
    for(int j=0; j < ny; ++j) {
      for(int i=0; i < nx; ++i) {
	f[i][j][k]=fin.read(1);
      }
    }
  }
  return f;
}

// return a 2D cross-section of a 3D file in index dir.
real[][] cut2(real[][][] f,
    int nx, int ny, int nz, int c, int dir) {
  real[][] f2;
  if(dir == 1) {
    f2=new real[nx][ny];
    for(int i=0; i < nx; ++i) {
      for(int j=0; j < ny; ++j) {
	f2[i][j]=f[i][j][c];
      }
    }
  }
  if(dir == 2) {
    f2=new real[ny][nz];
    for(int j=0; j < ny; ++j) {
      for(int k=0; k < nz; ++k) {
	f2[j][k]=f[c][j][k];
      }
    }
  }
  if(dir == 3) {
    f2=new real[nx][nz];
    for(int i=0; i < nx; ++i) {
      for(int k=0; k < nz; ++k) {
	f2[i][k]=f[i][c][k];
      }
    }
  }
  return f2;
}

// return a 1D cross-section of a 3D file in index dir.
real[] cut1(real[][][] f, int nx, int ny, int nz, int c1, int c2, int dir) {
  real [] f1;
  if(dir == 1) {
    f1=new real[nx];
    for(int i=0; i < nx; ++i)
      f1[i]=f[i][c1][c2];
  }
  if(dir == 2) {
    f1=new real[ny];
  for(int j=0; j < ny; ++j)
    f1[j]=f[c1][j][c2];
  }
  if(dir == 3) {
    f1=new real[nz];
    for(int k=0; k < ny; ++k)
      f1[k]=f[c1][c2][k];
  }
  return f1;
}


// Given a function y(x), find the fit to a+by for x in (start,stop).
// The data is returned in a pair (a,b).
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
