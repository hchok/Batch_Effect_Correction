function [K,normals,optim_vertices] = compute_exposing_hyperplanes(scaledF)
%compute convex hull of scaled data (it was pre-checked for being full rank)
K = convhulln(scaledF);

%identify vertices
vertices = unique(K(:));

%identify optimal vertices (partial order)
optim_vertices = [];
d = size(scaledF,2);
for i = 1:length(vertices)
    v = scaledF(vertices(i),:);
    tmp = scaledF(vertices,:);
    tmp(vertices(i),:) = [];
    diff = repmat(v,[size(tmp,1) 1]) - tmp;
    suboptimal_flag = sum(sum(diff <= 0,2) == d) > 0;
    if ~suboptimal_flag
        optim_vertices = [optim_vertices vertices(i)];
    end
end

%identify normals to all convex hull faces (convention: n*x_optim >= n*x_nonoptim)
normals = nan(size(K,1),d);
for i = 1:size(K,1)
    face = scaledF(K(i,:),:);
    v0 = face(1,:);
    tmp = face;
    tmp(1,:) = [];
    face_span = tmp - repmat(v0,[size(tmp,1) 1]);
    n = null(face_span);
    if size(n,2) == 1
        n = n/norm(n,1);
        offset = round(v0*n,10);
        other_vertices = find(~ismember(vertices,K(i,:)));
        non_optim_offset = round(max(scaledF(other_vertices,:)*n),10);
        if offset < non_optim_offset
            n = -n;
            non_optim_offset = -non_optim_offset;
            offset = -offset;
        end
        normals(i,:) = n';
    end
end
ndx = isnan(sum(normals,2));
normals(ndx,:) = [];
K(ndx,:) = [];
end











