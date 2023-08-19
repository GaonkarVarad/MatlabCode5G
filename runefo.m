function [res, par, sta, sys] = runefo(par, sta, sys)
% DESCRIPTION
% [res, par, sta, sys] = runefo(par, sta, sys)
% The basic rune dynamic function for OFDMA.
% INPUT  (all inputs are optional)
%  par --               BASIC SIMULATION PARAMETERS
%  par.frametime --     time interval in the simulation loop
%  par.nframes--        number of iterations in the main loop
%  par.slottime --      time per slot in seconds
%  par.nslots --        number of slots per frame
%  par.seed --          seed to all random sequencies in the simulation
%  par.cellradius --    cell radius [m]
%  par.km --            km^2+lm^2+km*lm => the number of sites
%  par.lm --            related to km above
%  par.sps --           number of sectors per site
%  par.gainconst --     gain at 1 meter distance [dB]
%  par.alpha --         distance attenuation coefficient
%  par.noise --         thermal noise floor [dBm]
%  par.sigma --         standard deviation for the lognormal fading [dB]
%  par.raa --           lognormal correlation down link (typical 0.5)
%  par.corrdist --      lognormal fading correlation distance [m]
%  par.offtraf --       average number of offered calls to a cell [Erlang/cell]
%  par.mht --           mean holding time [seconds]
%  par.amean --         average acceleration [m/s/s]
%  par.vmean --         average speed  [m/s]
%  par.pmax --          max transmit link power [dBm]
%  par.rbermax --       raw bit error level over which a call is dropped
%  par.usefastf --      make use of fast fading
%
%  par.packetlambdad -- probabilty of packet arrival downlink
%  par.packetsized --   average number of bits per packet downlink
%  par.packetlambdau -- probabilty of packet arrival uplink
%  par.packetsizeu      average number of bits per packet uplink
%  par.tau --           average time delay in the channel
%  par.freq --          carrier frequency
%  par.chbw --          chunk bandwitdh 
%  par.nchunks --       number of chunks over the badwidth
%  par.chunkduration -- chunk duration [s] Same as par.slottime


%  sta --               STATE variable that is changed at every iteration
%  sta.seed --          the seed before next iteration
%  sta.time --          sample time
%  sta.m --             mobile identity number
%  sta.mtop --          max mobile number used sofar
%  sta.xym --           position in complex form [m]
%  sta.xyv --           speed imag <=> north real <=> east [m/s]
%  sta.gmb --           the gain matrix (slow part)
%  sta.rgumbc --        rayleigh fading gain for the uplink
%  sta.rgdmbc --        rayleigh fading gain for the downlink
%
%  sta.pdmc --          transmitted power down link  [dBm]
%  sta.cdmc --          carrier downlink [dBm]
%  sta.idmc --          interference down link [dBm]
%  sta.sirdmc --        signal to interference ratio down link [dB]
%  sta.pumc --          transmitted power uplink  [dBm]
%  sta.cumc --          carrier uplink [dBm]
%  sta.iumc --          interference uplink [dBm]
%  sta.sirumc --        signal to interference ratio uplink [dB]
%
%  sys --               SYSTEM variables
%  sys.xyb --           base positions [m]
%  sys.fib --           cell center vector [m]
%  sys.rhombvec --      system vectors [m]
%  sys.raylmap --       rayleigh map [dB]
%  sys.raylmapvec --    rayleigh map vectors [m]
%  sys.lognmap --       lognormal map [dB]
%  sys.lognmapvec --    lognormal map vectors [m]
% OUTPUT
%  res --               cellarray of sta collected for each iteration of the simulation
%  par --               same as input otherwise created within the function
%  sta --               same structure as input otherwise created within the function
%  sys --               same as input otherwise created within the function
% SEE ALSO
% TRY
%  [res, par, sta, sys] = runefo

% by Magnus Almgren 080221

% set simulation parameter par if not present as an input
setifnotexistorempty('par',setparo) % default parameter setting
oseed = setseed(par.seed); % Set the seed in random generators

% Create the sys variable if not present as an input.
if ~exist('sys','var') | isempty(sys)
    % generate base station position and directions
    clear sys
    [sys.xyb, sys.fib, sys.rhombvec] = crecells(par.cellradius,par.sps,par.km,par.lm);

    % wrap coordinates so the fit will be better with positions of users when
    % plotted. This is only a cosmetic operation since wrapping is done anyway
    % when gain is calculated.
    sys.xyb = wrapinto(sys.xyb,sys.rhombvec,'hex');

    % Create a lognormal map. The lognormal map is dependent on the seed.
    oseed1 = setseed(par.seed);  % Set seed of pseudo random generator for the map.
    [sys.lognmap, sys.lognmapvec] = crelognmap(sys.xyb, sys.rhombvec, par.corrdist);

    % Create a rayleigh fading map
    freqs = flatten(par.freq + par.chbw*(1:par.nchunks),3); % make vektor point into 3rd dimension
    [sys.raylmap, sys.raylmapvec] = creraylmap(sys.rhombvec,freqs,par.tau);
    sys.xyboffsrayl = irandn(size(sys.xyb))*sys.rhombvec(1); % offset for fast fading

    setseed(oseed1);  % Restore seed to original value.
end

% Init of state variable
% The variables below are altered in the for loop and saved after each iteration
if ~exist('sta','var') | isempty(sta)
    clear sta
    e = zeros(0,1);
    sta.seed = par.seed; sta.time = 0; sta.m = e; sta.mtop  = 0; sta.xym = e; sta.xyv = e; sta.gmb = e;
    sta.nbitsdm = e; sta.nbitsum = e;
end

setseed(sta.seed); % Set seed in random generators.

% The simulation loop. One iteration corresponds to one frame
for iframe = 1:par.nframes
    dt = par.frametime-par.nslots*par.slottime;
    sta.time = sta.time+dt; % Timestamp

    % Terminate some calls due to hang up.
    terones = (rand(size(sta.xym)) < 1-exp(-mdiv(par.frametime,par.mht))); % Toss one biased Coin per active user

    % Make a realisation of new users.
    % At first iteration full traffic is generated
    incr = 1 - (sta.time~=dt).*exp(-mdiv(par.frametime,par.mht));
    % Average number of users to create. Poisson arrivals.
    nmob = poisson(par.offtraf .* size(sys.xyb,2).* incr);
    keep = ~terones;  % Clean out terminated users

    % Several of these updates are not necessary but is done in order to ease debugging
    update = inline('cate(1,partof(a,1,keep),nan(nmob,1))');
    updatez = inline('cate(1,partof(a,1,keep),zeros(nmob,1))');
    sta.xym   = update(sta.xym,keep,nmob);
    sta.xyv   = update(sta.xyv,keep,nmob);
    sta.m     = update(sta.m  ,keep,nmob);
    sta.gmb   = update(sta.gmb,keep,nmob);

    sta.nbitsdm  = updatez(sta.nbitsdm, keep,nmob);
    sta.nbitsum  = updatez(sta.nbitsum, keep,nmob);
    sta.mtop = sta.mtop+nmob; % highest mobile id so far

    % Move old users and initiate a position to new users.
    [sta.xym,sta.xyv] = mobmove(sta.xym,sta.xyv,par.amean,par.vmean,dt,sys.rhombvec);

    % Calculate the gain matrix. Size is terminals  by bases
    sta.gmb = pathgain(sta.xym, sys.xyb, sys.fib, par.sps, sys.rhombvec, ...
        par.gainconst, par.alpha, par.sigma, par.raa, ...
        sys.lognmap, sys.lognmapvec);

    for islot = 1:par.nslots
        dt = par.slottime;
        sta.time = sta.time + dt;

        % Move users one slot forward
        [sta.xym,sta.xyv] = mobmove(sta.xym,sta.xyv,par.amean,par.vmean,dt,sys.rhombvec);
        % Generate packets
        ispacketdm = poisson(adjsiza(par.packetlambdad,sta.nbitsdm));
        sta.nbitsdm = sta.nbitsdm + round(expdistr(par.packetsized.*ispacketdm));

        ispacketum = poisson(adjsiza(par.packetlambdau,sta.nbitsum));
        sta.nbitsum = sta.nbitsum + round(expdistr(par.packetsizeu.*ispacketum));

        % Compute flat fading gain
        sta.rgumbc = useraylmap(mplus(+sys.xyboffsrayl,-sta.xym), sys.raylmap, sys.raylmapvec);
        sta.rgdmbc = useraylmap(mplus(-sys.xyboffsrayl,-sta.xym), sys.raylmap, sys.raylmapvec);

        % Radio resource allocation
        isdesdmb = mand(ismax(sta.gmb,2),sta.nbitsdm>0); % Desired base station. Routing downlink
        isdesumb = mand(ismax(sta.gmb,2),sta.nbitsum>0); % Desired base station. Routing uplink
        sta.isdeslinkdmb = mand(ismax(mprod(isdesdmb,sta.nbitsdm),1),isdesdmb); % Select user. Scheduling bug
        sta.isdeslinkumb = mand(ismax(mprod(isdesumb,sta.nbitsum),1),isdesumb); % Select user.

        sta.pdm = par.pmax+lin2db(sum(sta.isdeslinkdmb,2)>0);
        sta.pum = par.pmax+lin2db(sum(sta.isdeslinkumb,2)>0);

        % Perform the transmission from the sending to the receiving side, all in dB or dBm.
        chdmbc = mplus(sta.gmb,lin2db(abs(sta.rgdmbc).^2));
        chdmuc = mplus(sta.gmb,lin2db(abs(sta.rgumbc).^2));
        [sta.cdmc, sta.idmc, sta.sirdmc] = transmit(sta.isdeslinkdmb, sta.pdm, chdmbc, par.noise, 2);
        [sta.cumc, sta.iumc, sta.sirumc] = transmit(sta.isdeslinkumb, sta.pum, chdmuc, par.noise, 1);
        % Calculate throughput 
        capdmc = par.chbw*par.chunkduration*log2(1+db2lin(sta.sirdmc));
        capumc = par.chbw*par.chunkduration*log2(1+db2lin(sta.sirumc));
        sta.nbitsdm = max(0,sta.nbitsdm - sum(capdmc,3)); % decrease user queue
        sta.nbitsum = max(0,sta.nbitsum - sum(capumc,3));

        sta.seed = setseed; % Save seed to faciliate simulation from this state
        res(islot,iframe)= sta; % Save result from this iteration
    end
    
    sta.seed = setseed; % Save seed to faciliate simulation from this state
end
setseed(oseed); % Restore seed to original value
1;

