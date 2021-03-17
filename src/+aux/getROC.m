function [FPR, TPR, AUC, THR] = getROC(SCORE, REAL)
  %GETROC Create an ROC curve from classification scores
  %   INPUT:
  %     SCORE: A matrix of (sample, class); with classification score for data
  %     REAL:  A binary matrix of (sample, class); that shows if the class is
  %       a real positive for that sample
  %   OUTPUT:
  %     CURVE: An array of (threshold, classifier, [fpr, tpr]) that
  %       contains the false positive rate and true positive rate for
  %       multiple thresholds.
  %     AUC: Area-Under-Curve, (class, [pos, neg]) for the ROC curve in question
  %     THRESH: The threshold values for the classification; where a score
  %       is determined as a positive if the score is >= the threshold
  %   Unlike matlab's ROC function; this function samles such that every
  %   transition is calculates
  
  % Get size of what we are dealing with
  N = size(SCORE, 1);
  C = size(SCORE, 2);
  
  % For each score; do individual ranking
  [rank_score, rank_order] = sort(SCORE, 1, 'descend');
  
  % For each score; do individual true/false stats
  rank_truth = REAL(sub2ind([N, C], rank_order, (1:C) .* ones(N, 1)));
  
  % Count true/false; including the case where threshold is Inf
  thrCls_tp = [zeros(1, C); cumsum( rank_truth, 1)];
  thrCls_fp = [zeros(1, C); cumsum(~rank_truth, 1)];
  thrCls_tn = [cumsum(~rank_truth, 1, 'reverse'); zeros(1, C)];
  thrCls_fn = [cumsum( rank_truth, 1, 'reverse'); zeros(1, C)];
  
  % Calculate the statistical variables
  thrCls_fpr = thrCls_fp ./ (thrCls_fp + thrCls_tn);
  thrCls_tpr = thrCls_tp ./ (thrCls_tp + thrCls_fn);
%   thrCls_fnr = 1 - thrCls_tpr;
%   thrCls_tnr = 1 - thrCls_fpr;
  
  % Output results for a curve plot
  FPR = thrCls_fpr;
  TPR = thrCls_tpr;
  
  % Record the thresholds on which the stats were taken
  THR = [inf(1, C); rank_score];
  
  % Do AUC calculation; each threshold adds area del_x * avg_y
  thrCls_delX = (thrCls_fpr(2:end, :) - thrCls_fpr(1:(end-1), :));
  thrCls_delY = (thrCls_tpr(2:end, :) + thrCls_tpr(1:(end-1), :)) .* .5;
  thrCls_delA = thrCls_delX .* thrCls_delY;
  AUC = sum(thrCls_delA, 1);
  
end

