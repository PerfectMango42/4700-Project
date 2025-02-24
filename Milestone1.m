set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'Defaultaxeslinewidth', 2)

set(0, 'DefaultFigureWindowStyle', 'docked')

c_c = 299792458;    % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;    % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;   % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;   % s^2/Fm
c_q = 1.60217653e-19;       % charge of electron
c_hb = 1.05457266913e-34;   % Dirac constant
c_h = c_hb*2*pi;            % 2 pi Dirac constant
g_fwhm = 3.53e+012/100;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.001;

RL = 0;                  % Left reflectivity coefficient
RR = 0;                  % Right reflectivity coefficient

% creating structure InputParasL and assigning values in the structure
% But InputParasR is just a regular scalar value
InputParasL.E0 = 1e15;   % Amplitude of electric field
InputParasL.we = 0;     % Frequency offset
InputParasL.t0 = 2e-12; % Time offset of Gaussian wave
InputParasL.wg = 5e-13; % Standard deviation of the wave
InputParasL.phi = 0;    % Starting phase of the wave
InputParasR = 0;        % No wave starting from the right

beta_r = 0;            % Real part of the detuning 
beta_i = 0;            % Imaginary part of the detuning

kappa0 = 100;           % Grating value
kappaStart = 1/3;
kappaStop = 2/3;

n_g = 3.5;              % index of refraction
vg = c_c/n_g*1e2;       % TWM cm/s group velocity
Lambda = 1550e-9;       % wavelength in nm 

plotN = 10;             % value to know which N you should plot the points for

L = 1000e-6*1e2;        % cm
XL = [0,L];             % X axis range in a matrix
YL = [-1*InputParasL.E0,1*InputParasL.E0];% Y axis range in a matrix

Nz = 200;               % total grid steps in the graph
dz = L/(Nz-1);          % spacial step size along the length (L)
dt = dz/vg;             % time step with the corresponding spacial step size
fsync = dt*vg/dz;       % Always equals 1, syncronizing normalized factor

Nt = floor(2*Nz);       % number of time steps (discrete number time points in simulation)
tmax = Nt*dt;           % total simulation time
t_L = dt*Nz;            % time to travel length

% nan(not a number) to be filled in later
z = linspace(0,L,Nz).'; % Nz points, Nz-1 segments (column vector of length Nz, each entry from 0 - L)
time = nan(1, Nt);      % Row vector of size 1 - Nt
InputL = nan(1,Nt);     % Row vector of size 1 - Nt
InputR = nan(1,Nt);     % Row vector of size 1 - Nt
OutputL = nan(1,Nt);    % Row vector of size 1 - Nt
OutputR = nan(1,Nt);    % Row vector of size 1 - Nt

% change the grating to not be a constant value, and instead slowely
% increase further into each grating


kappa = ones(size(z));
kappa(z<L*kappaStart) = 0; % no kappa for the first 3rd of waveguide
kappa(z>L*kappaStop) = 0; % no kappa for the last 3rd of waveguide

% Initializing Electric fields
Ef = zeros(size(z));    % z has Nz elements 
Er = zeros(size(z));    % z has Nz elements
Efprev = zeros(size(z));

% Initializing polarization values
Pf = zeros(size(z));
Pr = zeros(size(z));

% Setting previous values
Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

% Bit code generation
binary_input = [1 0 1 0 0 1 0 0 0 1 1 1]; % change sequence for grating 

% length of the input
length_input = length(binary_input);

% spacing of each grating on waveguide
l = L*kappaStart / length_input;

% each grating location is from L/3 to L/3 + l up to 2L/3
% loop through the binary input and setup the gratings for 1 values
% The kappa for each grating increases further down the waveguide,
% additionally each kappa within each grating greadually increases
for i = 1:length_input

    % each part of the waveguide where the grating could go
    start_position = L*kappaStart + l*(i-1);
    end_position = L*kappaStart + l*i;

    % the space along the waveguide to set the kappa value
    index_kappa = (z >= start_position & z < end_position);
    
    % minimum and maximum kappa values
    max_kappa = 125;
    min_kappa = 100;

    % Create an exponential growth trend over the waveguide
    kappa_trend = min_kappa + (max_kappa - min_kappa) * (exp(linspace(0, 1, length_input)) - 1) / (exp(1) - 1);

    % Doesnt matter about 1 becuase kappa is already full of 1
    if binary_input(i) == 0
        % put no grating here
        kappa(index_kappa) = 0;
    else

        % smooth transition within the grating
        z_smooth = z(index_kappa);  % Extract the relevant z values

         % Ensure we don't go out of bounds
        if i < length_input
            kappa_end = kappa_trend(i+1);
        else
            % Keep constant for the last grating
            kappa_end = kappa_trend(i); 
        end

        % gradual increase within the grating
        kappa_values = linspace(kappa_trend(i), kappa_end, length(z_smooth)); 
        
        % Assign gradually increasing kappa values within the grating
        kappa(index_kappa) = kappa_values;
    end
end


% SourceFct creates a function handle (allow you to pass functions as
% arguments to other functions, store them...)
Ef1 = @SourceFct;
ErN = @SourceFct;

t = 0;                  % Starting time for the simulation
time(1) = t;            % already initialized as nan matrix, assigns first value to t

% nan values that get assigned values from the SourceFct with the InputParas structures, and the current time
InputL(1) = Ef1(t, InputParasL); 
InputR(1) = ErN(t, InputParasR);

% nan values that get assigned the boundaries of forward and reverse electric fields 
OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);

% Assigning the electric field input values at the boundaries opposite to
% the outputs
Ef(1) = InputL(1);
Er(Nz) = InputR(1);   

% Create a figure that contains three sublots, each display different data
% about the electric fields inputs and outputs
% Forward propagating electric field at the left side
figure('name', 'Fields')
subplot(3,1,1)
plot(z*10000,real(Ef),'r');
hold off
xlabel('z(\mum)')
ylabel('E_f')
% Reverse propagating electric field at the right side
subplot(3,1,2)
plot(z*10000,real(Er),'b');
xlabel('z(\mum)')
ylabel('E_r')
hold off
% Inputs and Outputs over time in picoseconds
subplot(3,1,3)
plot(time*1e12, real(InputL), 'r'); hold on
plot(time*1e12, real(OutputR), 'r--'); 
plot(time*1e12, real(InputR), 'b'); hold on
plot(time*1e12, imag(OutputL), 'b--');
xlabel('time(ps)')
ylabel('E')
hold off

beta = ones(size(z))*(beta_r+1i*beta_i);    % Creates a matrix of 1's of specified size times the detuning term
exp_det = exp(-1i*dz*beta);                 % Creates the exponential gain/loss according to detuning term (array of exponential terms)

for i = 2:Nt        % Iterate from 2 to the number of time steps
    t = dt*(i-1);   % Determine next time according to spacial step size and current iteration
    time(i) = t;    % Increment time
    
    % nan values that get assigned values from the SourceFct with the InputParas structures, and the current time
    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, InputParasR);

    Ef(1) = InputL(i) + RL*Er(1);       % Adding reflectivity coefficients
    Er(Nz) = InputR(i) + RR*Ef(Nz);     % Adding reflectivity coefficients

    % Start with zero polarization
    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    % Calculate polarization coefficients
    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    %determine the new forward and reverse fields based on the polarization
    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));

    % Updates the current Ef and Er over the spatial grid, ensuring to
    % normalize (fsync scaling) and accounts for the exponential gain/loss according to detuning parameters
    % does element wise multiplication of the corresponding exponential
    if i >= 3
        % Efprev has been populated
        Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i*dz*kappa(2:Nz).*Er(2:Nz);
        Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Efprev(1:Nz-1);
    else
        % Efprev has not been populated
        Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i*dz*kappa(2:Nz).*Er(2:Nz);   
        Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Ef(1:Nz-1);
    end 

    Efprev(2:Nz) = Ef(2:Nz);    % populating Efprev

    % nan values that get assigned the boundaries of forward and reverse electric fields
    OutputR(i) = Ef(Nz)*(1-RR);     % Adding the loss from the mirrors reflectivity
    OutputL(i) = Er(1)*(1-RL);      % Adding the loss from the mirrors reflectivity

    % Create the plots that visualize the forward and reverse propagating
    % electric fields
    if mod(i,plotN) == 0            % updates every plotN iterations
        % Real and imaginary parts of forward propagating wave
        subplot(2,2,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        % Real and imaginary parts of the reverse propagating wave
        subplot(2,2,2)
        plot(z*10000,real(Er),'b'); hold on
        plot(z*10000,imag(Er),'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re','\Im')
        hold off
        % Input and Output signals over time
        subplot(2,2,3);
        plot(time*1e12, real(InputL), 'r'); hold on
        plot(time*1e12, real(OutputR), 'g'); 
        plot(time*1e12, real(InputR), 'b');
        plot(time*1e12, imag(OutputL), 'm');
        xlim([0,Nt*dt*1e12])            % Sets total simulation time (number of time steps * length of each step)
        ylim(YL)                        % Ensures plot is big enough for electric field amplitude
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off

        subplot(2,2,4);
        plot (z*1000, kappa, 'r');
        xlabel('z(\mum)')
        ylabel('kappa')
        hold off
        pause(0.01)                     % Short delay in iterations of for loop (sleep())
    end
    % resets the previous values
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
end

% Create the Fourier Transform calculation 
fftOutputR = fftshift(fft(OutputR));
fftOutputL = fftshift(fft(OutputL));
omega = fftshift(wspace(time));

fftInputL = fftshift(fft(InputL));

% Plot the Fourier Transform results as frequency
% Plot of right output vs time
figure('name', 'Fields')
subplot(4,1,1)
plot(time*1e12, real(OutputR),'m'); hold on
plot(time*1e12, real(OutputL),'b');
xlabel('time(ps)')
ylabel('Output')
legend('Right', 'Left')
hold off

% Plot of absolute value of fourier transform vs omega
subplot(4,1,2)
plot(omega,abs(fftOutputR),'y'); hold on
plot(omega,abs(fftOutputL),'r');
xlabel('omega')
ylabel('FFT |E|')
legend('Right', 'Left')
hold off
% Inputs and Outputs over time in picoseconds
subplot(4,1,3)
plot(omega,unwrap(angle(fftOutputR)),'g'); hold on
plot(omega,unwrap(angle(fftOutputL)),'c');
xlabel('omega')
ylabel('FFT Angle')
legend('Right', 'Left')
hold off

% Magnitude of input vs output in frequency domain
subplot(4,1,4)
plot(omega,abs(fftOutputR),'g'); hold on
plot(omega,abs(fftInputL),'b');
xlabel('omega')
ylabel('FFT |E|')
legend('Output', 'Input')
hold off
