function [Kii,Kib,Kbi,Kbb] = GetSubK(Ks,i,ne,nd)
   bids = Bids(i,ne,nd);
   bl = length(bids);

   ir = 1:ne-bl;
   br = ne-bl+1:ne;

   K = Ks(:,:,i);

   Kii = K(ir,ir);
   Kib = K(ir,br);
   Kbi = K(br,ir);
   Kbb = K(br,br);
end
