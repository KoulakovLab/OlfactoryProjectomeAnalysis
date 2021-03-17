function o = plus(i1, i2)
  %PLUS Overload the addition function to merge data
  o = mapseqData;
  
  % Merge source regions
  if isequal(i1.nSrcRegSli, i2.nSrcRegSli) && isequal(i1.srcRegName, i2.srcRegName)
    o.srcImg = vertcat(i1.srcImg, i2.srcImg);
    o.nSrcRegSli = i1.nSrcRegSli;
    o.srcRegName = i1.srcRegName;
    if (~isempty(i1.different_srcRegSum)) & (~isempty(i2.different_srcRegSum))
      o.srcRegSum = vertcat(i1.srcRegSum, i2.srcRegSum);
    end
  else
    error('Incompatible slice numbers');
  end
  
  % Merge projection regions
  if isequal(i1.nPrjRegSli, i2.nPrjRegSli) && isequal(i1.prjRegName, i2.prjRegName)
    o.prjImg = vertcat(i1.prjImg, i2.prjImg);
    o.nPrjRegSli = i1.nPrjRegSli;
    o.prjRegName = i1.prjRegName;
    if (~isempty(i1.different_prjRegSum)) & (~isempty(i2.different_prjRegSum))
      o.prjRegSum = vertcat(i1.prjRegSum, i2.prjRegSum);
    end
  else
    error('Incompatible slice numbers');
  end
  
  % Check if the brcName field is same
  if isequal(i1.brcName, i2.brcName) && (size(i1.brcId, 2) == size(i2.brcId, 2))
    o.brcName = i1.brcName;
    o.brcId = vertcat(i1.brcId, i2.brcId);
    if ~(isempty(i1.brcCor) || isempty(i2.brcCor))
      o.brcCor = vertcat(i1.brcCor, i2.brcCor);
    end
  else
    warning('Incompatible barcode labels; leaving labels empty');
  end
  
  % Check if data has compatible fields
  field_list = unique([fieldnames(i1.data); fieldnames(i2.data)]);
  for i = 1:length(field_list)
    curf = field_list{i};
    if isfield(i1.data, curf) && isfield(i2.data, curf)
      if any(strcmp(curf(1:3), {'brc', 'src', 'prj'}))
        % Check if the secondary dimension is compatible
        if size(i1.data.(curf), 2) == size(i2.data.(curf), 2)
          o.data.(curf) = vertcat(i1.data.(curf), i2.data.(curf));
        else
          warning(['Incompatible ', curf]);
        end
      else
        warning(['Incompatible data fields; not merging ', curf, '.']);
      end
    else
      warning(['Incompatible data fields; not merging ', curf, '.']);
    end
  end
  
  % Merge info
  o.name = char(join({i1.name, i2.name}, ' & '));
  
end