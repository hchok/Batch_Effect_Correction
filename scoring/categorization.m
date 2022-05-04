function [centroid,constituents,constituents_category,vols] = categorization(varargin)
%process input arguments
input_file = varargin{1};
aggregator = varargin{2};
if length(varargin) > 2
    lb = varargin{3};
    ub = varargin{4};
    %check whether or not specified scoring space is empty
    flag = check_feasibility(lb,ub);
    feasibility_flag = flag == 1;
else
    lb = [];
    ub = [];
    feasibility_flag = true;
end

if feasibility_flag
    %read in data file
     T = readcell(input_file);

    %set aggregator
    if strcmpi(aggregator,'min')
        op = @min;
    elseif strcmpi(aggregator,'med')
        op = @median;
    end

    %process data into a sample-by-feature matrix
    num_pts = size(T,1) - 1;
    F = nan(num_pts,7);%presumed data dimension is 7
    for i = 1:num_pts
        datum1 = str2num(T{i+1,2});
        F(i,1) = op(datum1(1:2:end));%dim 1 @ sample i
        F(i,2) = op(datum1(2:2:end));%dim 2 @ sample i
        datum2 = str2num(T{i+1,3});
        F(i,3) = op(datum2(1:2:end));%dim 3 @ sample i
        F(i,4) = op(datum2(2:2:end));%dim 4 @ sample i
        datum3 = str2num(T{i+1,4});
        F(i,5) = op(datum3);%dim 5 @ sample i 
        datum4 = str2num(T{i+1,5});
        F(i,6) = op(datum4);%dim 6 @ sample i
        datum5 = str2num(T{i+1,6});
        F(i,7) = op(datum5);%dim 7 @ sample i    
    end

    %scale data
    data = (F - min(F)).*repmat((max(F) - min(F)).^(-1),[size(F,1) 1]);

    %identify problem dimension
    d = size(data,2);

    %define positive orthant constraints
    positive_orthant_constraints_A = -eye(d);
    positive_orthant_constraints_b = zeros(d,1);

    %define box constraints (if existent)
    if ~isempty(lb)
        box_constraints_A = [+eye(d);-eye(d)];
        box_constraints_b = [+ub';-lb'];
    end

    %define simplex hyperplane        
    simplex_hyperplane_Aeq = ones(1,d);
    simplex_hyperplane_beq = 1;

    %compute all vertex-exposing hyperplanes
    [K,normals,optim_vertices] = compute_exposing_hyperplanes(data);

    %compute all (constrained) polyhedral cones exposing optimum vertices
    vols = nan(1,length(optim_vertices));
    constituents = [];
    for category = optim_vertices
        %identify all hyperplanes exposing optimum vertex associated with above category
        hyperplane_ndx = transpose(find(sum(ismember(K,category),2) > 0));

        %identify normals of all hyperplanes above
        N = normals(hyperplane_ndx,:);

        %compute polyhedral cone defined over above normals
        [exposing_cone_A,exposing_cone_b,exposing_cone_Aeq,exposing_cone_beq] = vert2lcon([N;zeros(1,d)]);
        exposing_cone_A(abs(exposing_cone_b) > 1e-3,:) = [];%remove top hyperplane
        exposing_cone_b(abs(exposing_cone_b) > 1e-3,:) = [];%remove top hyperplane

        %define intersection of exposing cone, standard simplex, and (possibly) box constraints
        if ~isempty(lb)
            constrained_exposing_cone_A = [exposing_cone_A;positive_orthant_constraints_A;box_constraints_A];
            constrained_exposing_cone_b = [exposing_cone_b;positive_orthant_constraints_b;box_constraints_b];
        else
            constrained_exposing_cone_A = [exposing_cone_A;positive_orthant_constraints_A];
            constrained_exposing_cone_b = [exposing_cone_b;positive_orthant_constraints_b];        
        end

        %compute category constituents, volumes, and centroids
        V{category} = lcon2vert(constrained_exposing_cone_A,constrained_exposing_cone_b,[exposing_cone_Aeq;simplex_hyperplane_Aeq],[exposing_cone_beq;simplex_hyperplane_beq]);
        new_constituents{category} = round(V{category},4);
        try
            [~,vols(category)] = convhulln([V{category};zeros(1,d)]);
            cone = [V{category};zeros(1,d)];
            T = delaunayn(cone);
            centroid{category} = 0;
            for i = 1:size(T,1)
                [~,vol] = convhulln(cone(T(i,:),:));
                local_centroid = mean(cone(T(i,:),:));
                local_centroid = local_centroid/norm(local_centroid,1);
                centroid{category} = centroid{category} + local_centroid*vol/vols(category);
            end
        catch %if convhulln/delaunayn computation fails, category has no volume
            vols(category) = 0;
            centroid{category} = [];
        end
        constituents = [constituents; new_constituents{category}];
    end
    vols = vols/sum(vols);
    constituents = unique(constituents,'rows');
    for i = 1:length(new_constituents)
        if ~isempty(new_constituents{i})
            constituents_category{i} = ismember(constituents,new_constituents{i},'rows');
        else
            constituents_category{i} = [];
        end
    end
else
    fprintf('Specified box constraints are not feasible!\n');
    centroid = [];constituents = [];constituents_category = [];vols = [];
end
end
%%
function [flag] = check_feasibility(lb,ub)
cond1 = min((ub - lb)) >= 0;
cond2 = min(lb) >= 0;
cond3 = max(ub) <= 1;
cond4 = sum(lb) <= 1;
cond5 = sum(ub) >= 1;
flag = logical(cond1*cond2*cond3*cond4*cond5);
end