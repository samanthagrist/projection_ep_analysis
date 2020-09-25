function Iscale = sgautoscale(I)
I_max = max(I(:));
Iscale = uint8(double(255*(double(I)./double(I_max))));
end