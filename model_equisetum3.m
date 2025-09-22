function out = model_equisetum3(par, data)
    % Unpack parameters from vector
    ak      = par(1);        % fractionation factor for diffusion of water vapour through air
    akh     = ak^(2/3);      % kinetic fractionation factor for transport across boundary layer
    qk      = par(2);        % 17O exponent (diffusion of water vapour through air)
    ea      = par(3);        % vapor pressure (air) 
    es      = par(4);        % vapor pressure at leaf surface
    y       = par(5);        % transpiration tuning exponent
    T       = par(6);        % Temperature (celsius) 
    d18Omw  = par(7);        % d18O of source water (VSMOW scale)
   % m       = par(8);        % slope, bitches!
    n       = max(1, round(par(8)));         % slope, bitches!

    % Constants
    R18_SMOW = 0.0020052;
    R17_SMOW = 0.0003799;
    q17 = 0.529;
    % Number of stem segments
    
    %n = length(data.xdata);
    delta18 = zeros(n,1);
    delta17 = zeros(n,1);
    D17     = zeros(n,1);
        
    % --- T-dependent terms ---
    TK = 273.1 + T; % celsius to Kelvin
    aeq = exp((1.137*10^6/(TK^2)-4.156*10^2/TK-2.0667)/1000);
    % aeq = 
    esat =  1 / 1000 * exp(34.494-4924.99/(T+237.1))/(T+105)^(1.57); %saturation vapour pressure of air (kPa)
    ei = esat;  % intercellular vapor pressure is fully saturated
    aeq17 = aeq^q17; % equilibrium water-vapour fractionation for 17O 
    ak17 = ak^(qk); % fractionation during water diffusion in air for 17O
    akh17 = akh^(qk); % effective fractionation during transport through boundary layer

    %initialise stem water
    Rx = R18_SMOW * (1 + d18Omw / 1000 ); 
    Rx17 = R17_SMOW * ( Rx / R18_SMOW )^q17;
    Ra = Rx / aeq;
    Ra17 = Rx17 / aeq17; 
    
    % Loop over segments
    for i = 1:n
        n_remaining = n - i;
        x_i = (1 / (1 + n_remaining))^y;
    
        % --- 1. Chain-of-lakes mixing plus Flanagan leaf model ---
        beta = ei * ak - es * ak - ea * akh + es * akh;
        Rseg = ( Rx + ea * Ra * x_i  / beta ) / ( 1 - x_i + ei * x_i / ( beta * aeq ) ); 

        % --- 2. Enabled for 17O ---
        beta17 = ei * ak17 - es * ak17 - ea * akh17 + es * akh17;
        R17seg =  ( Rx17 + ea * Ra17 * x_i  / beta17 ) / ( 1 - x_i + ei * x_i / ( beta17 * aeq17 ) ); 
        %R17seg = R17_SMOW * ( Rseg / R18_SMOW )^qk; % alternative
        %implementation
        
        % --- 3. Calculate delta values ---
        delta18(i) = 1000 * (Rseg / R18_SMOW - 1);
        delta17(i) = 1000 * (R17seg / R17_SMOW - 1);

        % --- 4. Calculate D'17O ---
        d18_ln = 1000 * reallog( Rseg / R18_SMOW );
        d17_ln = 1000 * reallog( R17seg / R17_SMOW );
        D17(i) = 1000 * (d17_ln - 0.528 * d18_ln);  
        
        Rx = Rseg;     % infer Rl, treat as new Rx
        Rx17 = R17seg;   % same for 17O
    end

    % Return model predictions for ?18O and D'17O
    out = [delta18, D17];
end
