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
%  par.RBbw --          RB bandwitdh 
%  par.nRB --           number of RBs over the badwidth
%  par.RBduration --    RB duration [s] Same as par.slottime


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
    
    %%plothex(sys.xyb,sys.fib); %% plot the cells

    % Create a lognormal map. The lognormal map is dependent on the seed.
    oseed1 = setseed(par.seed);  % Set seed of pseudo random generator for the map.
    [sys.lognmap, sys.lognmapvec] = crelognmap(sys.xyb, sys.rhombvec, par.corrdist);

    % Create a rayleigh fading map
    freqs = flatten(par.freq + par.RBbw*(1:par.nRB),3); % make vektor point into 3rd dimension
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
for iframe = 1:par.nframes % Outer loop is for FRAMES
    dt = par.frametime-par.nslots*par.slottime;
    sta.time = sta.time+dt; % Timestamp

    % Terminate some calls due to hang up.
    terones = (rand(size(sta.xym)) < 1-exp(-mdiv(par.frametime,par.mht))); % Toss one biased Coin per active user

    % Make a realisation of new users.
    % At first iteration full traffic is generated
    
    incr = 1 - (sta.time~=dt).*exp(-mdiv(par.frametime,par.mht)); %% TODO: How to set par.mht ??
   
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
    
    for islot = 1:par.nslots % Inner loop is for SLOTS
        
        % new time = old time + slot time
        sta.time = sta.time + par.slottime; 

        % Move users one slot forward
        [sta.xym,sta.xyv] = mobmove(sta.xym,sta.xyv,par.amean,par.vmean,dt,sys.rhombvec);
       
        % Generate packets for DL
        ispacketdm = poisson(adjsiza(par.packetlambdad,sta.nbitsdm));
        sta.nbitsdm = sta.nbitsdm + round(expdistr(par.packetsized.*ispacketdm));
        
        %Generate packets for UL
        ispacketum = poisson(adjsiza(par.packetlambdau,sta.nbitsum));
        sta.nbitsum = sta.nbitsum + round(expdistr(par.packetsizeu.*ispacketum));

        % Compute flat fading gain
        sta.rgdmbc = useraylmap(mplus(-sys.xyboffsrayl,-sta.xym), sys.raylmap, sys.raylmapvec); 
        sta.rgumbc = useraylmap(mplus(+sys.xyboffsrayl,-sta.xym), sys.raylmap, sys.raylmapvec);
       
        %%%%%%%%%%%%  Assign users to BS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   gmb = N_ue * N_bs 
        %   step1: choose the best BS for all the user ismax(sta.gmb,2) , return 1
        %   step2: if it has something in the queue, return 1
        %   Returns: For e.g: 
        %   isdesdmb:   1 0 0  means User 1 has to receive something from BS 1
        %               0 0 0  means User 2 has nothing to receive (queue is empty)
        %               0 0 1  means User 3 has to receive something from BS 3
        
        isdesdmb = mand(ismax(sta.gmb,2),sta.nbitsdm>0); % Desired base station. Routing downlink
        isdesumb = mand(ismax(sta.gmb,2),sta.nbitsum>0); % Desired base station. Routing uplink

        %calculate sinr for users that has to transmit
%         chdmbc0=mplus(sta.gmb,lin2db(abs(sta.rgdmbc).^2));
%         pw = par.pmax+lin2db(sum(isdesdmb,2)>0); 
%         [cdmc0, idmc0, sirdmc0] = transmit(isdesdmb, pw, chdmbc0, par.noise, 2);
%         capdmc0 = par.RBbw*par.nRB*log2(1+db2lin(sirdmc0)); 
%         %capdmc0=sum(capdmc0,3);
%         sched_matrix= capdmc0<0;
%         sched_matrix=sched_matrix+0;
        
        
             %%%%%%%%%%%%%%% Assign RBs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  step 1: Multiply isdesdmb by nbitsdm
        %  step 2: Find the user with Maximum bits to transmit Dim 1: ismax(mprod(isdesdmb,sta.nbitsdm),1)
        %  step 3: Multiply with isdesdmb 
        %  Returns: isdeslinkdmb: 
        %           1 0 0 User 1 is scheduled
        %           0 0 0
        %           0 1 0 User 3 is scheduled 
        
        
%         if(iframe==1&&islot==1)
%            %b=sta.nbitsdm<0;
%            %b=b+0;
%            tx_bits = sta.nbitsdm<0;
%            avg_rate=sta.nbitsdm<0;
%            avg_rate=avg_rate+1; % just to ignore the infinite value in metric 
%            %avg_rate=mplus(0.999.*avg_rate,mprod(capdmc0,alloc_var).*0.001);
%         end
%         
%         avg_rate=avg_rate+tx_bits;
%         metric_pf=mdiv(capdmc0,avg_rate/(iframe*islot));
        %metric_pf=mdiv(sta.nbitsdm+1,avg_rate/(iframe*islot));
       
%         tot_bs = size(sta.gmb,2);
%         for bs = 1:tot_bs
%             tot_users_to_be_scheduled = sum(isdesdmb(:,bs),1);
%             users_to_be_scheduled = find(isdesdmb(:,bs)==1);
%             if(tot_users_to_be_scheduled==1)
%                    sched_matrix(users_to_be_scheduled,:)=1;
%             else %% divide the RB
%                     for user= 1:tot_users_to_be_scheduled
%                         [M,I] = max(capdmc0(users_to_be_scheduled(user),:));
%                         sched_matrix(users_to_be_scheduled(user),1,I)=1;
%                     end
%                 
%                 end
%              
%             end
%           sched_matrix=sum(sched_matrix,3)>0;
% 
%         sta.isdeslinkdmb= mprod(isdesdmb,sched_matrix);
        sta.isdeslinkdmb = mand(ismax(mprod(isdesdmb,sta.nbitsdm),1),isdesdmb); % Select user.
        
        %sta.isdeslinkdmb = mand(ismax(mprod(isdesdmb,metric_pf),1),isdesdmb); % Proportional fair
        %sta.isdeslinkdmb = mand(ismax(mprod(isdesdmb,sta.nbitsdm),1),isdesdmb);% Max data 
        %sta.isdeslinkdmb = mand(ismax(mprod(isdesdmb,db2lin(sta.gmb)),1),isdesdmb); % MAx gain
        sta.isdeslinkumb = mand(ismax(mprod(isdesumb,sta.nbitsum),1),isdesumb); % Select user.
        
        
        
        %%%%%%%%%% just for debug purpose to find scheduling percentage
        %a=find(sum(sta.isdeslinkdmb,2)==1);
        %for i = 1:length(a)
        %b(a(i),1)=b(a(i),1)+1;
        %end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % set the transmitting power only for those scheduled to transmit
        sta.pdm = par.pmax+lin2db(sum(sta.isdeslinkdmb,2)>0);  %% returns -inf in nothing to transmit else returns 33
        sta.pum = par.pmax+lin2db(sum(sta.isdeslinkumb,2)>0);  

        % Perform the transmission from the sending to the receiving side, all in dB or dBm.
        chdmbc = mplus(sta.gmb,lin2db(abs(sta.rgdmbc).^2)); %% adds the fastfading gain to the channel N_UE*N_BS*N_RBs
        chdmuc = mplus(sta.gmb,lin2db(abs(sta.rgumbc).^2));
        
        %chdmbc= mplus(chdmbc,sched_matrix);
        
        
        
        %%%%%% Transmission and Reception %%%%%%%%%%%%%%%%%%%%%%%
       
        [sta.cdmc, sta.idmc, sta.sirdmc] = transmit(sta.isdeslinkdmb, sta.pdm, chdmbc, par.noise, 2);
        [sta.cumc, sta.iumc, sta.sirumc] = transmit(sta.isdeslinkumb, sta.pum, chdmuc, par.noise, 1);
        
        % Calculate throughput 
        capdmc = par.RBbw*par.RBduration*log2(1+db2lin(sta.sirdmc)); % capdmc FORMAT: N_UE * 1 * N_RBS
        capumc = par.RBbw*par.RBduration*log2(1+db2lin(sta.sirumc)); % capumc FORMAT: N_UE * 1 * N_RBs
        %tx_bits = mplus(tx_bits,sum(capdmc,3));
        %sta.tp= sum(capdmc,3);
        sta.nbitsdm = max(0,sta.nbitsdm - sum(capdmc,3)); % DEQUEUE [Residual Bits] = [Bits in queue] - [Bits transmitted]
        sta.nbitsum = max(0,sta.nbitsum - sum(capumc,3));

        sta.seed = setseed; % Save seed to faciliate simulation from this state
        res(islot,iframe)= sta; % Save result from this slot
    end
    
    sta.seed = setseed; % Save seed to faciliate simulation from this state
end
setseed(oseed); % Restore seed to original value
%bar(b); % to see Scheduling percentage
%ylim([0 1]);
%xlabel('USER ID #');
%ylabel('total slots scheduled');
%title('Scheduled Opportunity')
1;

