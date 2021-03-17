classdef mapseqDataAligned < handle
  %MAPSEQDATAALIGNED Storage and analysis of mapseq data after alignment
  properties
    % Global
    name
    % Source (injection site) region stuff
    srcRegName              % Region names
    nSrcRegSli              % Number of slices per region
    srcImg
    % Projections region stuff
    prjRegName
    nPrjRegSli
    prjImg
    % Barcode label name
    brcName
    brcId
    brcCor
    % Hold extra info
    data
  end
  properties(Dependent)
    % Variables that are dynamically calculated
    nData                   % Number of data points
    nBrc                    % Number of barcodes per point
    nSrcReg                 % Number of regions in the injection site
    nSrcSli                 % Number of slices in injection site
    srcRegInd               % Indices per region of each slice
    nPrjReg
    nPrjSli
    prjRegInd
    prjRegSum
    srcRegSum
  end
  properties(Hidden)
    different_prjRegSum
    different_srcRegSum
  end
  
  methods
    function obj = mapseqDataAligned(arr)
      %MAPSEQDATA Create instance from a mapseqData array
      if isa(arr, 'mapseqData') && isvector(arr)
        obj.srcRegName = arr(1).srcRegName;
        obj.nSrcRegSli = arr(1).nSrcRegSli;
        obj.prjRegName = arr(1).prjRegName;
        obj.nPrjRegSli = arr(1).nPrjRegSli;
        % Check if all are equal
        if ~all(arrayfun( @(x) ...
            isequal(x.srcRegName, obj.srcRegName) && ...
            isequal(x.nSrcRegSli, obj.nSrcRegSli) && ...
            isequal(x.prjRegName, obj.prjRegName) && ...
            isequal(x.nPrjRegSli, obj.nPrjRegSli), arr))
          error('Input mapseqData array is not aligned');
        end
        obj.srcImg = {arr.srcImg}';
        obj.prjImg = {arr.prjImg}';
        obj.brcId  = {arr.brcId }';
        obj.brcCor  = {arr.brcCor }';
        obj.name = {arr.name}';
        if all(arrayfun(@(x) ~isempty(x.different_prjRegSum), arr))
          obj.prjRegSum = {arr.prjRegSum}';
        end
        if all(arrayfun(@(x) ~isempty(x.different_srcRegSum), arr))
          obj.srcRegSum = {arr.srcRegSum}';
        end
      else
        error('Constructor need a mapseqData array');
      end
      obj.data = struct();
    end
    
    % GLOBAL
    function val = get.nData(obj)
      % Total number of barcodes
      val = length(obj.name);
    end
    function val = get.nBrc(obj)
      % Total number of barcodes
      val = cellfun(@(x) size(x, 1), obj.srcImg);
    end
    
    % SOURCE SITE
    function val = get.nSrcSli(obj)
      % Total number of slices
      val = sum(obj.nSrcRegSli);
    end
    function val = get.nSrcReg(obj)
      % Total number of regions
      val = length(obj.srcRegName);
    end
    function val = get.srcRegInd(obj)
      % List of indices to access a region with
      val = arrayfun(@(x, y) (1:x) + y - x, ...
        obj.nSrcRegSli, cumsum(obj.nSrcRegSli), ...
        'UniformOutput', 0);
    end
    
    % PROJECTION SITE
    function val = get.nPrjSli(obj)
      % Total number of slices
      val = sum(obj.nPrjRegSli);
    end
    function val = get.nPrjReg(obj)
      % Total number of regions
      val = length(obj.prjRegName);
    end
    function val = get.prjRegInd(obj)
      % List of indices to access a region with
      val = arrayfun(@(x, y) (1:x) + y - x, ...
        obj.nPrjRegSli, cumsum(obj.nPrjRegSli), ...
        'UniformOutput', 0);
    end
    
    % REGION SUMS
    function val = get.prjRegSum(obj)
      % Collapse regions into one slice
      if isempty(obj.different_prjRegSum)
        val = cell(obj.nData, 1);
        for d = 1:obj.nData
          val{d} = zeros(obj.nBrc(d), obj.nPrjReg);
          for r = 1:obj.nPrjReg
            val{d}(:, r) = sum(obj.prjImg{d}(:, obj.prjRegInd{r}), 2);
          end
        end
      else
        val = obj.different_prjRegSum;
      end
    end
    function set.prjRegSum(obj, val)
      obj.different_prjRegSum = val;
    end
    function val = get.srcRegSum(obj)
      % Collapse regions into one slice
      if isempty(obj.different_srcRegSum)
        val = cell(obj.nData, 1);
        for d = 1:obj.nData
          val{d} = zeros(obj.nBrc(d), obj.nSrcReg);
          for r = 1:obj.nPrjReg
            val{d}(:, r) = sum(obj.srcImg{d}(:, obj.srcRegInd{r}), 2);
          end
        end
      else
        val = obj.different_srcRegSum;
      end
    end
    function set.srcRegSum(obj, val)
      obj.different_srcRegSum = val;
    end
    
  end
end

