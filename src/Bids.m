function bids = Bids(i,ne,nd)
    % Boundary ids of subdomain i
    if i==1
        bids = ne;
    elseif i==nd
        bids = 1;
    else
        bids = [1,ne];
    end
end
