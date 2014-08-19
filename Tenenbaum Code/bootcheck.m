function r=bootcheck(gamma,model,data)

y = model;

inds = find(1-isnan(y));

B = zeros(1000,1);

if ((length(inds)>0)&(var(y(inds))>eps))
  for j = 1:1000
    samp = unidrnd(length(inds),length(inds),1);
    g = fminsearch('powertransform',1,[],y(inds(samp)),data(inds(samp)));
    B(j) = -powertransform(g,y(inds(samp)),data(inds(samp)));
  end
  [val ind] = sort(B);
  r = [val(26) val(975) std(B)];
else
  r = 0;
end

