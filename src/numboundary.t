
template<class BoundX, class BoundY>
void SymmetricBoundary<BoundX, BoundY>::apply(NumMatrix<double,2> &u) const
{   
  int lx = u.getLow(0);
  int mx = u.getHigh(0);
  int ly = u.getLow(1);
  int my = u.getHigh(1);

  BoundX bx(lx,mx);
  BoundY by(ly,my);

  const int fx = bx.factor();
  const int fy = by.factor();

  for(int i = lx; i <= mx; i++) {
    u(i,by.low())  = fy*u(i,by.lowSrc());
    u(i,by.high()) = fy*u(i,by.highSrc());
  }

  for(int j = ly; j <= my; j++) {
    u(bx.low(),j)  = fx*u(bx.lowSrc(),j);
    u(bx.high(),j) = fx*u(bx.highSrc(),j);
  }
}

template<class BoundLeft, class BoundRight, class BoundBottom, class BoundTop>
void MixedBoundary<BoundLeft, BoundRight, BoundBottom, BoundTop>
  ::apply(NumMatrix<double,2> &u) const
{   
  int lx = u.getLow(0);
  int mx = u.getHigh(0);
  int ly = u.getLow(1);
  int my = u.getHigh(1);

  BoundLeft   bl(lx,mx);
  BoundRight  br(lx,mx);
  BoundBottom bb(ly,my);
  BoundTop    bt(ly,my);

  const int fl = bl.factor();
  const int fr = br.factor();
  const int ft = bt.factor();
  const int fb = bb.factor();

  for(int i = lx; i <= mx; i++) {
    u(i,bb.low())  = fb*u(i,bb.lowSrc());
    u(i,bt.high()) = ft*u(i,bt.highSrc());
  }

  for(int j = ly; j <= my; j++) {
    u(bl.low(),j)  = fl*u(bl.lowSrc(),j);
    u(br.high(),j) = fr*u(br.highSrc(),j);
  }
}


template<class BoundLeft, class BoundRight, class BoundBottom, class BoundTop>
void MixedBoundaryWithOffset<BoundLeft, BoundRight, BoundBottom, BoundTop>
  ::apply(NumMatrix<double,2> &u) const
{   
  int lx = u.getLow(0);
  int mx = u.getHigh(0);
  int ly = u.getLow(1);
  int my = u.getHigh(1);

  BoundLeft   bl(lx,mx);
  BoundRight  br(lx,mx);
  BoundBottom bb(ly,my);
  BoundTop    bt(ly,my);

  const int fl = bl.factor();
  const int fr = br.factor();
  const int ft = bt.factor();
  const int fb = bb.factor();

  for(int i = lx; i <= mx; i++) {
    u(i,bb.low())  = fb*u(i,bb.lowSrc()) + ob;
    u(i,bt.high()) = ft*u(i,bt.highSrc()) + ot;
  }

  for(int j = ly; j <= my; j++) {
    u(bl.low(),j)  = fl*u(bl.lowSrc(),j) + ol;
    u(br.high(),j) = fr*u(br.highSrc(),j) + or;
  }
}

