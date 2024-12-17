

pt = 1e-6;

x0 = [1, .4, 1.5, 0, 6.7]';

[h0,dh_dx] = h_cart(x0,true);
[h2,~] = h_cart(x0 + [0; pt; 0; 0; 0],false);
[h3,~] = h_cart(x0 + [0; 0; pt; 0; 0],false);


dh_dx2 = (h2 - h0)/pt
dh_dx3 = (h3 - h0)/pt

dh_dx

t = []; % s

x0 = rand(5, 1);
v0 = [0;0];


[fscript0,dfscript_dx,dfscript_dvtil] = ...
    fscript_cart(t,x0,[],v0,true);

A = nan(length(x0));
for i = 1:length(x0)

    dx = zeros(size(x0));
    dx(i) = pt;

    [fscript,~,~] = ...
        fscript_cart(t,x0 + dx,[],[0;0],false);

    A(:, i) = (fscript - fscript0)/pt;

end

D = nan(length(x0), length(v0));
for i = 1:length(v0)

    dv = zeros(size(v0));
    dv(i) = pt;

    [fscript,~,~] = ...
        fscript_cart(t,x0,[],v0 + dv,false);

    D(:, i) = (fscript - fscript0)/pt;

end









