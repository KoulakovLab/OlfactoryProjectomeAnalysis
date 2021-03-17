function output = getDataset(obj, slice)
%GETDATASET Get the requested dataset from the aligned objects.

if ischar(slice)
    slice = {slice};
end
N = length(slice);

if iscell(slice)
    sto = zeros(length(slice), 1);
    for i = 1:length(slice)
        ind = find(strcmp(obj.name, slice{i}), 1);
        if isempty(ind)
            error('No such dataset ''', slice, ''' in dataset');
        end
        sto(i) = ind;
    end
    slice = sto;
elseif isvector(slice) && all(round(slice) == slice, 'all') && all(slice > 0, 'all')
    if any(slice > obj.nData, 'all');
        error('No such dataset');
    end
else
    error('Invalid slice specification');
end

output(length(slice), 1) = mapseqData;
for i = 1:length(slice)
    output(i).name = obj.name{slice(i)};
    output(i).srcImg = obj.srcImg{slice(i)};
    output(i).srcRegName = obj.srcRegName;
    output(i).nSrcRegSli = obj.nSrcRegSli;
    output(i).prjImg = obj.prjImg{slice(i)};
    output(i).prjRegName = obj.prjRegName;
    output(i).nPrjRegSli = obj.nPrjRegSli;
    output(i).brcId = obj.brcId{slice(i)};
    output(i).brcCor = obj.brcCor{slice(i)};
    output(i).brcName = obj.brcName;
    if ~isempty(obj.different_prjRegSum)
      output(i).prjRegSum = obj.prjRegSum{slice(i)};
    end
    if ~isempty(obj.different_srcRegSum)
      output(i).srcRegSum = obj.srcRegSum{slice(i)};
    end
    % Grab classifier if exist
    if isfield(obj.data, 'classifier')
      output(i).data.classifier = obj.data.classifier;
    end
    % Grab TSNE coordinates if exist
    if isfield(obj.data, 'brcTsne')
      output(i).data.brc_tsne = obj.data.brcTsne{slice(i)};
    end
    % Grab the template stuff if exists
    if isfield(obj.data, 'temp')
      output(i).data.temp = obj.data.temp;
    end
end

end

