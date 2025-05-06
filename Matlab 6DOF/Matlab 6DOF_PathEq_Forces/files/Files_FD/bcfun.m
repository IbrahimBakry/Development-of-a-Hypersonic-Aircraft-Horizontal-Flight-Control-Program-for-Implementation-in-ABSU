function res = bcfun(yleft,yright)
Initial_conitions;
yl = [Vtot0 gamma0 chi0 h0 phi0 theta0 mf0 psiV0 psigamma0...
    psichi0 psih0 psiphi0 psitheta0 psim0];
yr = [1830 0 0 30000 55.58/57.3 37.94/57.3 1.6329e+04 0 0 0 0 0 0 0];

yleft(:,1) = yl;
yright(:,end) = yr;

res = [yleft, yright];
