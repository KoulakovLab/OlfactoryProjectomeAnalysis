function CLS = runClass(RES, CON)
  %RUNCLASS Do classification based on Voronoi control points
  [~, CLS] = min(pdist2(RES, CON), [], 2);
end

