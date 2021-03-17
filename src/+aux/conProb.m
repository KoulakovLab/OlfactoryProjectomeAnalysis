function [p_sta_con, p_sta_con_alt] = conProb(sta, con, wmet, elim)
%BAYES_X_GIVEN_Y Calculates conditional probability; given a measurement in
%the condition; what are the probabilities of finding in state
% INPUTS:
% |_sta:
% |     The measurements of the event that we want the probability
% |     distribution of
% |_con:
% |     The measurements of the event that we expect to be the given
% |     condition
% \_ The inputs are expected to be indicative of some event counts; with
%     the following format(sample, outcomes)
%    Both sets must have the same number of measurement samples
% OUTPUTS:
% |_p_sta_con:
% |     Conditional probability. Gives the probability that if the event in
% |     set "con" is measured; the probability of measuring the specific
% |     event in "sta".
% |     Has dimension (outcome of sta, outcome of con)
% |_p_sta_con_alt:
% |     Alternate measurements of the aforementioned values for sanity
% |    check
% \_  If there are unmeasured events; it's normalized to 0.
if ~exist('wmet', 'var')
  wmet = 'uniform';
end
if ~exist('elim', 'var')
  elim = 0;
end

% Preliminaries
removeNan = @(x) fillmissing(x, 'constant', 0);
if size(sta, 1) ~= size(con, 1)
    error('Sample sizes must match');
end
% Remove disproportionately bright samples
if elim > 0
  [~, sto1] = sort(sum(sta, 2), 'descend');
  [~, sto2] = sort(sum(con, 2), 'descend');
  sto = unique([sto1(1:elim); sto2(1:elim)]);
  sta(sto, :) = [];
  con(sto, :) = [];
end

% Axiom of choice (or whatever;) assume we pick samples equally
switch wmet
  case {'uniform', 'equal'}
    ps_sample = ones(size(sta, 1), 1) / size(con, 1);
    pc_sample = ones(size(con, 1), 1) / size(con, 1);
  case {'brightness'}
    ps_sample = sum(sta, 2) ./ sum(sta, 'all');
    pc_sample = sum(con, 2) ./ sum(con, 'all');
  otherwise
    error('Invalid weighting method');
end

% Get the P(con/sta; sample) the measurement prob given the sample
p_sta_sample = removeNan(sta ./ sum(sta, 2))';
p_con_sample = removeNan(con ./ sum(con, 2))';

% Utilize chain rule; P(A) = sum_b P(A|B) P(B)
% Get the P(con/sta);
p_sta = p_sta_sample * ps_sample;
p_con = p_con_sample * pc_sample;

% Utilize bayes rule: P(A|B) P(B) = P(B|A) P(A)
% or P(A|B) = P(A) * P(B|A) / P(B)
% Get P(sample; con/sta)
p_sample_sta = removeNan(ps_sample .* ((p_sta_sample ./ p_sta)'));
p_sample_con = removeNan(pc_sample .* ((p_con_sample ./ p_con)'));

% Use another chain rule; P(A|B) = sum_c P(A|C) P(C|B)
% Calculate the probability for P(sta_con)
p_sta_con = p_sta_sample * p_sample_con;
p_con_sta = p_con_sample * p_sample_sta;

% Use the bayes rule to calculate the alternatives for sanity checks
p_sta_con_alt = removeNan(p_sta .* (p_con_sta') ./ (p_con'));

end

