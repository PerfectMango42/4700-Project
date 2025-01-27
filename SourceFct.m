% Compute electric field signal based on time and input parameters 
function E = SourceFct(t, InputParas)   % Defines function name for other MATLAB file to call

if isfield(InputParas, 'rep')           % if input parameters contains 'rep'
    n = floor(t/InputParas.rep);        % How many repatitions have occured
    t = t-n*InputParas.rep;             % Update time based on current repetition cycle
end

if ~isstruct(InputParas)
    E = InputParas;                     % Handle if incorrect input parameters format
else
    % Compute the electric field with the input parameters values
    E = InputParas.E0*exp(-(t-InputParas.t0)^2/InputParas.wg^2)*exp(1i*(InputParas.we*t + InputParas.phi));
end

end
