function oarr = partition(obj, cond)
%PARTITION Partition the data set into multiple depending on condition
% criteria
%   INPUT:
%       - obj:  The input object
%       - cond: The condition to check towards. Can either be a label
%       vector; or a function handle.
%   OUTPUT:
%       - oarr: Each slice that corresponds to a different index
LAB = zeros(obj.nBrc, 1);

if isa(cond, 'function_handle')
    LAB = cond(obj);
elseif isvector(cond) && isnumeric(cond) && (length(cond) == o.nBrc)
    % Check if valid label vector
    if all(round(cond) == cond, 'all') && all(cond > 0, 'all')
        LAB = cond;
    end
else
    error('Invalid label specification.');
end

NAME = unique(LAB);
oarr(length(NAME), 1) = mapseqData();

for i = 1:length(NAME)
    THIS = (LAB == NAME(i));
    oarr(i).name = [obj.name, ' label: ', num2str(NAME(i))];
    oarr(i).srcImg = obj.srcImg(THIS, :);
    oarr(i).prjImg = obj.prjImg(THIS, :);
    oarr(i).brcId = obj.brcId(THIS, :);
    oarr(i).brcCor = obj.brcId(THIS, :);
    oarr(i).nSrcRegSli = obj.nSrcRegSli;
    oarr(i).srcRegName = obj.srcRegName;
    oarr(i).nPrjRegSli = obj.nPrjRegSli;
    oarr(i).prjRegName = obj.prjRegName;
    if ~isempty(obj.different_prjRegSum)
      oarr(i).prjRegSum = obj.prjRegSum(THIS, :);
    end
    if ~isempty(obj.different_srcRegSum)
      oarr(i).srcRegSum = obj.srcRegSum(THIS, :);
    end
    cur = obj.data;
    dnames = fields(cur);
    for k = (find(contains(dnames, 'brc', 'IgnoreCase', 1))')
        cur.(dnames{k}) = cur.(dnames{k})(THIS, :);
    end
    oarr(i).data = cur;
    
end

