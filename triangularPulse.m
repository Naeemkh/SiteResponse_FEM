function [x,y]=triangularPulse(a,b,c,x)

y = ((a<x)&(x<=b)).*(x-a)./(b-a)+((b<x)&(x<c)).*(c-x)./(c-b)+(x<=0)*0+(x>=c)*0;


end