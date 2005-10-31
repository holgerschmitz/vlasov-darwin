BEGIN{
  Lx = 512;
  Ly = 256;
  
  for (i=1; i<=sqrt(Lx*Ly); ++i)
  {
    ComSize = i;
    dim1 = int(sqrt( Lx*(ComSize+0.001) / Ly ));
    if (dim1<1) dim1=1;
    
    dim2 = int(ComSize/dim1);
    if (dim2<1)
    {
      dim2=1;
      dim1=ComSize;
    }
    
    if (ComSize==(dim1*dim2))
    {
      print ComSize, dim1,"x",dim2;
    }
  }
  
}
