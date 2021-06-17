function [ shift ] = simpleshiftfit(y1,y2,y1err,y2err,shift_guess)
%Author: Tim Gorman
%Date: 2017-06-08
%Title: simpleshiftfit
%This function shifts the 2nd of two arrays to the first based on a least
%squares fitting.  The shift is additionally weighted by the error bars of
%the data involved. Same domain is assumed and NaNs are ignored
    opts1=  optimset('display','off');
    if isempty(y1err)==1 && isempty(y2err)==1 % no error bars input
        fitfunc =@(shift)sqrt((y1 -(shift+y2)).^2);
    elseif isempty(y1err)==0 && isempty(y2err)==0 % both errors not empty
        fitfunc =@(shift)sqrt((y1-(shift+y2)).^2./(y1err.^2+y2err.^2));
    elseif isempty(y1err)==1 && isempty(y2err)==0 % first error is empty
        fitfunc =@(shift)sqrt((y1-(shift+y2)).^2./y2err.^2);
    elseif isempty(y1err)==0 && isempty(y2err)==1 % second error is empty
        fitfunc =@(shift)sqrt((y1-(shift+y2)).^2./y1err.^2);
    end
    
    shift = lsqnonlin(fitfunc,shift_guess,[],[],opts1);
    
end

