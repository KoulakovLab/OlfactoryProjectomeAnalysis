function reorderBarcodes(obj, order)
%REORDERBARCODES Establish a new barcode ordering with the specified index.
%   The only difference between this and manually doing is it also orders
%   any 'brc' specifier in the data field

obj.srcImg = obj.srcImg(order, :);
obj.prjImg = obj.prjImg(order, :);
if ~isempty(obj.brcId)
  obj.brcId = obj.brcId(order, :);
end
if ~isempty(obj.brcCor)
  obj.brcCor = obj.brcCor(order, :);
end
if ~isempty(obj.different_srcRegSum)
  obj.srcRegSum = obj.srcRegSum(order, :);
end
if ~isempty(obj.different_prjRegSum)
  obj.prjRegSum = obj.prjRegSum(order, :);
end

names = fields(obj.data);

for i = 1:length(names)
    if length(names{i}) > 2
        if strcmp('brc',  names{i}(1:3))
            obj.data.(names{i}) = obj.data.(names{i})(order, :);
        end
    end
end

end

