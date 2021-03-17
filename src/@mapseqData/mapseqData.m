classdef mapseqData < handle
  %MAPSEQDATA Storage and analysis of mapseq data
  %   The information on MAPSEQ results are saved in an instance of the
  %   mapseqData object, and analyzed using the methods described
  
  properties
    % The name(s) of the dataset
    name
    % Source (injection site) region stuff
    srcImg                  % Imaging data for the region
    srcRegName              % Region names
    nSrcRegSli              % Number of slices per region
    % Projections region stuff
    prjImg
    prjRegName
    nPrjRegSli
    % Barcode information
    brcId                   % Identification per barcode
    brcName
    brcCor                  % Some coordinate system for dim. red. barcodes
    % Just an extra struct array to hold information
    data
  end
  properties(Dependent)
    % Variables that are dynamically calculated
    nBrc                    % Number of barcodes
    nSrcReg                 % Number of regions in the injection site
    nSrcSli                 % Number of slices in injection site
    srcRegInd               % Indices per region of each slice
    nPrjReg
    nPrjSli
    prjRegInd
    % Other things
    brcSrcIpr
    brcPrjIpr
    brcSrcDist
    brcPrjDist
    brcSrcVis
    brcPrjVis
    % Region sums
    prjRegSum
    srcRegSum
  end
  properties(Hidden)
    different_prjRegSum
    different_srcRegSum
  end
  
  methods
    function obj = mapseqData(varargin)
      %MAPSEQDATA Construct an instance of this class
      obj.data = struct();
    end
    
    % GLOBAL
    function val = get.nBrc(obj)
      % Total number of barcodes
      val = size(obj.srcImg, 1);
    end
    
    % SOURCE SITE
    function val = get.nSrcSli(obj)
      % Total number of slices
      val = size(obj.srcImg, 2);
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
    function val = get.brcSrcIpr(obj)
      % IPR of barcodes in the source region
      val = aux.ipr(obj.srcImg);
    end
    function val = get.brcSrcDist(obj)
      % Avg. distance of barcodes in the projection region
      val = (obj.srcImg * ((1:obj.nSrcSli)')) ./ sum(obj.srcImg, 2);
    end
    function val = get.brcSrcVis(obj)
      % IPR of barcodes in the projection region
      [~, val] = max(obj.srcImg, [], 2);
    end
    function val = get.srcRegSum(obj)
      % Collapse regions into one slice
      if isempty(obj.different_srcRegSum)
        val = zeros(obj.nBrc, obj.nSrcReg);
        for r = 1:obj.nSrcReg
          val(:, r) = sum(obj.srcImg(:, obj.srcRegInd{r}), 2);
        end
      else
        val = obj.different_srcRegSum;
      end
    end
    function set.srcRegSum(obj, val)
      obj.different_srcRegSum = val;
    end
    
    % PROJECTION SITE
    function val = get.nPrjSli(obj)
      % Total number of regions
      val = size(obj.prjImg, 2);
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
    function val = get.brcPrjIpr(obj)
      % IPR of barcodes in the projection region
      val = aux.ipr(obj.prjImg);
    end
    function val = get.brcPrjDist(obj)
      % Avg. distance of barcodes in the projection region
      val = (obj.prjImg * ((1:obj.nPrjSli)')) ./ sum(obj.prjImg, 2);
    end
    function val = get.brcPrjVis(obj)
      % Very important slice of barcodes in the projection region
      [~, val] = max(obj.prjImg, [], 2);
    end
    function val = get.prjRegSum(obj)
      % Collapse regions into one slice
%       if isempty(obj.different_prjRegSum)
        val = zeros(obj.nBrc, obj.nPrjReg);
        for r = 1:obj.nPrjReg
          val(:, r) = sum(obj.prjImg(:, obj.prjRegInd{r}), 2);
        end
%       else
%         val = obj.different_prjRegSum;
%       end
    end
    function set.prjRegSum(obj, val)
      obj.different_prjRegSum = val;
    end
  end
end

