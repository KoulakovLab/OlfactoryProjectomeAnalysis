function out = uplus(inp)
%UPLUS Unary plus on MAPSEQDATA returns a hard-copy of the current object

out = mapseqData;
out.name = inp.name;
out.srcImg = inp.srcImg;
out.srcRegName = inp.srcRegName;
out.nSrcRegSli = inp.nSrcRegSli;
out.prjImg = inp.prjImg;
out.prjRegName = inp.prjRegName;
out.nPrjRegSli = inp.nPrjRegSli;
out.brcId = inp.brcId;
out.brcName = inp.brcName;
out.brcCor = inp.brcCor;
out.data = inp.data;
if ~isempty(inp.different_prjRegSum)
  out.prjRegSum = inp.prjRegSum;
end
if ~isempty(inp.different_srcRegSum)
  out.srcRegSum = inp.srcRegSum;
end

end

