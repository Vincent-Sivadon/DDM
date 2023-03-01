function [x,xs] = NonOverlappingFields(ne,nd,ni)
    % xs represent the field at interfaces
    % x is nd vectors of ne elements
    % x(:,:,i) is the subdomain's local vector
    x  = zeros(ne,1,nd);
    xs = zeros(ni,1);
end
