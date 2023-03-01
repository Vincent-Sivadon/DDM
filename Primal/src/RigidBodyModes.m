function Rbs = RigidBodyModes(Sps,ni,nd)
    Rbs = zeros(2,2,nd);

    for i=1:nd
        Rbs(:,:,i) = pinv(Sps(:,:,i));
    end
end
