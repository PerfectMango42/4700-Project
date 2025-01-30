function w = wspace(t,nt);

if (nargin<2)
    nt = length(t);
    dt = t(2) - t(1);
    t = t(nt) - t(1) + dt;
end

if (nargin == 2)
    dt = t/nt;
end

w = 2*pi*(0:nt-1)'/t;
kv = find(w >= pi/dt);
w(kv) = w(kv) - 2*pi/dt;