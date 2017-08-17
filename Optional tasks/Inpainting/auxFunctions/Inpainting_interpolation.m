function [u] = Inpainting_interpolation(f, mask)
    data = f;
    data(logical(mask)) = NaN;
    u = fixgaps(data);
end

% FIXGAPS Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN
% in the input time series (may be complex), but ignores
% trailing and leading NaN.
function y=fixgaps(x)
    % R. Pawlowicz 6/Nov/99
    y=x;
    bd=isnan(x);
    gd=find(~bd);
    bd([1:(min(gd)-1) (max(gd)+1):end])=0;
    y(bd)=interp1(gd,x(gd),find(bd));
end