
template<class BoundX, class BoundY>
void MixedBoundary<BoundX, BoundY>::apply(NumMatrix<double,2> &u) const
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

