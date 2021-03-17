function r = getPrjSliReg(o, s)
%GETSLIREG Get the region ID of the slice asked
r = find(s <= cumsum(o.nPrjRegSli), 1);
end

