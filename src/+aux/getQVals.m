function [outArg] = getQVals(inpArg)
%GETQVALS Calculates q values

inpSize = size(inpArg);
outArg = zeros(inpSize);
inpNonZeroFlag = inpArg(:) > 0;
inpNonZero = inpArg(inpNonZeroFlag);
if sum(inpNonZeroFlag) > 1
  try
    [~, outArg(inpNonZeroFlag)] = mafdr(inpNonZero);
  catch
    outArg(inpNonZeroFlag) = inpNonZero;
  end
else
  outArg(inpNonZeroFlag) = inpNonZero;
end

end

