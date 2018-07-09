function [true_p,p,y,a] = qint(xm1,ym1,x0,y0,xp1,yp1)
% QINT: Quadratic interpolation of 3 uniformly spaced samples
%
%               [p,y,a] = qint(ym1,y0,yp1)
%
%       returns extremum-location p, height y, and half-curvature a
%       of a parabolic fit through three points.
%       The parabola is given by y(x) = a*(x-p)^2+b,
%       where y(-1)=ym1, y(0)=y0, y(1)=yp1.

% TAKES AS INPUTS: three EQUALLY SPACED x-values (xm1, x0, xp1) and the
% three corresponding y-values (ym1, y0, yp1)

    peak = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
    
   if nargout>1
   p = peak;
   end;
   if nargout>2
     y = y0 - 0.25*(ym1-yp1)*p;
   end;
   if nargout>3
     a = 0.5*(ym1 - 2*y0 + yp1);
   end;
   
   true_p = peak*(x0-xm1)+x0;


% Source: https://www.dsprelated.com/freebooks/sasp/Matlab_Finding_Interpolated_Spectral.html#sec:peakml
