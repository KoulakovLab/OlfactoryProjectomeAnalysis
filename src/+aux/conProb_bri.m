function p_sta_con = conProb_bri(sta, con)
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
MASK = 5;

% Preliminaries
removeNan = @(x) fillmissing(x, 'constant', 0);
if size(sta, 1) ~= size(con, 1)
    error('Sample sizes must match');
end

% Axiom of choice (or whatever;) assume we pick samples equally
% p_sample = ones(size(sta, 1), 1) / size(con, 1);
% p_sample = sum([con, sta], 2) ./ sum([con, sta], 'all');
pc_sample = sum(con, 2) ./ sum(con, 'all');
ps_sample = sum(sta, 2) ./ sum(sta, 'all');

% Eliminate barcodes
[~, sto_crank] = sort(pc_sample, 'descend');
[~, sto_srank] = sort(ps_sample, 'descend');
sto_elim = unique([sto_crank(1:MASK); sto_srank(1:MASK)]);
sta(sto_elim, :) = [];
con(sto_elim, :) = [];
pc_sample = sum(con, 2) ./ sum(con, 'all');
clear('ps_sample');

% Get the P(con/sta; sample) the measurement prob given the sample
p_sta_sample = removeNan(sta ./ sum(sta, 2))';
p_con_sample = removeNan(con ./ sum(con, 2))';

% Utilize chain rule; P(A) = sum_b P(A|B) P(B)
% Get the P(con/sta);
% p_sta = p_sta_sample * ps_sample;
p_con = p_con_sample * pc_sample;

% Utilize bayes rule: P(A|B) P(B) = P(B|A) P(A)
% or P(A|B) = P(A) * P(B|A) / P(B)
% Get P(sample; con/sta)
% p_sample_sta = removeNan(ps_sample .* ((p_sta_sample ./ p_sta)'));
p_sample_con = removeNan(pc_sample .* ((p_con_sample ./ p_con)'));

% Use another chain rule; P(A|B) = sum_c P(A|C) P(C|B)
% Calculate the probability for P(sta_con)
p_sta_con = p_sta_sample * p_sample_con;
% p_con_sta = p_con_sample * p_sample_sta;

% Use the bayes rule to calculate the alternatives for sanity checks
% p_sta_con_alt = removeNan(p_sta .* (p_con_sta') ./ (p_con'));

end

