function [cl, val, abu, IND, ind, scl] = watershedc(C, h)
% WATERSHEDC Clustering algorithm using watershed and forest-fire algorithm
N = size(C,1);
cl = zeros(N,1);

[r,c]=find(C);              % rows and columns
v = C(find(C(:)));          % values of correlation

if h>0
    maxv = max(v);
    mask = (v<(maxv-h));
    v = v.*mask+~mask*(maxv-h);     % apply marker
    
    ind = find(~mask);
    Cm = sparse(r(ind),c(ind),v(ind), size(C,1), size(C,1));
    
    cl = ff(Cm);
    
%     figure
%     plot(cl, '.-')
end

ind = r>c;                  % lower diagonal only
r = r(ind);
c = c(ind);
v = v(ind);

[v, ind] = sort(v, 'descend');
c = c(ind);
r = r(ind);

next_cluster_number=max(cl)+1;

for i=1:length(c)
    c1 = c(i);
    c2 = r(i);
    if (cl(c1)==0)&(cl(c2)==0)
        cl(c1) = next_cluster_number;
        cl(c2) = next_cluster_number;
        next_cluster_number = next_cluster_number+1;
    end
    if (cl(c1)==0)&(cl(c2)~=0)
        cl(c1) = cl(c2);
    end
    if (cl(c1)~=0)&(cl(c2)==0)
        cl(c2) = cl(c1);
    end    
end

[val,abu] = values2(cl);
[abu,ind] = sort(abu,'descend');
val = val(ind);

count=1;
ind = zeros(size(cl));
scl = ind;
IND = cell([length(val),1]);

for i=1:length(val)
    iii = find(cl==val(i));
    ind(count:(count+length(iii)-1)) = iii;
    scl(count:(count+length(iii)-1)) = val(i);
    count = count+length(iii);
    IND{i}=iii;
end

    function cl = ff(C)
    %   function cl = ff(C)
    %	Forest-fire algorithm with the matrix of links
    %	vv is the vector of values in matrix L
    %	c is the vector of cluster numbers
        C0 = (C>0);
        NN = size(C,1);
        cl=zeros(NN,1);                  % List of clusters
        cc=1;                           % Current cluster number
        n=0;
        for cn=1:NN
            fraction_done = cn/NN;
            %disp(char(repmat(8,1,50)))
            if cl(cn)==0				% Not classified yet
                % Start new cluster
                vv = sparse(cn, 1, 1, NN, 1);
                fv = cn;
                dn=1;
                nstep=0;
                while(dn>0)
                    nstep=nstep+1;
                    nv = C0*vv+vv;
                    fnv = find(nv);
                    dn = length(fnv)-length(fv);
                    vv=nv;
                    fv = fnv;
                end
                if nstep>1 % this is to exclude cluster with only one
                    % element
                    indd = find(vv);
                    cl(indd)=cc;
                    cc = cc+1;
                end
                fraction_done = cn/NN;
            end
        end
        return
    end

    function [val, abu] = values2(y)
        yy=sort(y(:));
        dy=diff(double(yy));
        indd=find([1; dy; 1]);
        val=yy(indd(1:(length(indd)-1)));
        abu=diff(indd);
    end
end