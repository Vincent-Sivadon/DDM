function iids = Iids(bids,ne)
    iids = 1:ne;
    iids(bids) = [];
end
