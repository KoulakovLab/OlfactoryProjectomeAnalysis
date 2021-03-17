function progressbar(cur, tot)
%PROGRESSBAR Print a commandline progress bar.
%   If cur is 1 (or 0); prints a full new line
%   If cur is not equal to tot; updates the previous line
%   If cur is equal to tot; also prints a newline
persistent temp
persistent inum
ICON = '|/-\';
LEN = 50;
PRC = round(100 * cur / tot);
NDOT = round(LEN * cur / tot);
FSTR = ['[', repmat('.', 1, NDOT), repmat(' ', 1, LEN-NDOT), ']', ...
    '%s(%3d%%)\n'];
ELEN = LEN + 10;
DELAY = .1;

if cur <= 1
    % Start recording time
    temp = tic;
    inum = 1;
    fprintf(['\n', FSTR], ICON(mod(inum - 1, length(ICON)) + 1), PRC);
    drawnow('update');
elseif cur == tot
    temp = [];
    inum = 0;
    fprintf([repmat('\b', 1, ELEN), FSTR], ' ', PRC);
    drawnow('update');
    % Stop recording time
elseif toc(temp) >= DELAY
    temp = tic;
    inum = inum + 1;
    fprintf([repmat('\b', 1, ELEN), FSTR], ...
        ICON(mod(inum - 1, length(ICON)) + 1), PRC);
    drawnow('update');
end


end

