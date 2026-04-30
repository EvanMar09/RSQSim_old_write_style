function tau = initialTau(x,z, tauFname)
% this is ridiculous in that it loads from the file
% each time it is called

tmp = load(tauFname);
xx = tmp(:,1);
zz = tmp(:,2);
ttau = tmp(:,3);

dist = sqrt((x-xx).^2 + (z-zz).^2);

tau = mean(ttau(find(dist < 1.1*min(dist))));

end

