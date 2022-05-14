function [ u_out, r_found ] = my_point_location_problem( mpsol, theta )

% mpsol = solution.mpsol1;
% theta = ones(6,1);
P = mpsol{1}.Partition;
Nr = length(mpsol{1}.Fi);

u_out = [];
r_found = [];

r = 1;
while ( isempty(r_found) == true ) && r <= Nr
    flagIsInside = P(r).contains(theta);
    if flagIsInside == 1
        r_found = r;
    end
    r = r + 1;
end

u_out = mpsol{1}.Fi{r_found} * theta + mpsol{1}.Gi{r_found};

end

