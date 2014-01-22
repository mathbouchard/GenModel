nc = 60; 
nr = 10;

A = sparse(nr,nc);

c = zeros(1,nc);
ct = zeros(1,nc);
for i=1:nr
    for j=1:nc
	if(rand() > 0.9)
	    A(i,j) = 1.0;
	end
    end
end

lb = zeros(1,nc);
ub = zeros(1,nc);
ub = ub+1.0;

rt = zeros(nr,1);
rhs = ones(nr,1);

for j=1:nr
    rt(j) = 'G';
end

for i=1:nc
    c(i) = i;
    ct(i) = 'C';
end

max = 0

[sol rc dual slack cstat rstat B] = gensolve(A, c, lb, ub, ct, rt, rhs, max);
sol
I = -eye(nr);
AA=full([A(:,cstat==1) I(:,rstat==1)]);
full(B*AA)
full(B*rhs)
