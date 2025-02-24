function w = wspace(t,nt); % time and number of time steps
% Generate a vector of angular frequencies

if (nargin<2) % if only one parameter was passed to wspace
    % Ensure the time represents the entire duration of nt time steps
    nt = length(t);
    dt = t(2) - t(1);
    t = t(nt) - t(1) + dt;
end

if (nargin == 2) % if both parameters were passed to wspace
    % calculate dt as the total time divided by numbre of time steps
    dt = t/nt;
end

w = 2*pi*(0:nt-1)'/t; % angular frequency vector w
kv = find(w >= pi/dt); % find the frequencies of w greater than pi/dt
w(kv) = w(kv) - 2*pi/dt; % wraps the higher frequencies back below pi/dt