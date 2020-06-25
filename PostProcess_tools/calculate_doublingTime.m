%Calculating doubling times within SuperSegger

function output = intGet3D( clist, g3D_ind )

output = [];
if iscell( clist )
    nc = numel(clist);
    
    for ii = 1:nc
        output = cat(1, output, intGet3D( clist{ii}, g3D_ind ));
    end
else
    output = clist.data3D(:,g3D_ind,:);
end

end

%extract Time in frames and Length and on that data fit exponentials 
%and do the rest of the calculations and finally add that result to the
%clist