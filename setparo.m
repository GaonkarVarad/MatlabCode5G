function par = setparo
% DESCRIPTION
% par = setparo
%  Sets the structured parameter par that is the input to runefo.
%  This file is intended to use as a template.
% OUTPUT      
%  par --  input parameter to rune
% SEE ALSO
%  runefc, runefg, runeft
% TRY
%  setparo

% by Magnus Almgren 080221

% Set up the default for the gain part
par=setparg;

% Simulation specific paramters (Classical LTE: frametime = 0.010, slottime = 0.001, nslots = 10)
par.frametime = 1e-2; % time interval between simulated frames. Could be nslots*slottime 
par.nframes = 1;       % number of iterations in the main loop 
par.slottime = 1e-3;   % time per slot 
par.nslots  = 10;      % number of slots per frame
par.seed = 1;          % seed to all random sequencies in the simulation
		       
% Transmission specific parameters 
par.homargin = 3;      % gain margin between two bases used at Hand Off
par.pmax = 33;         % power set to each new link
par.usefastf = 1;      % make use of fast fading

% specifics for OFDM
par.packetlambdad = 0.1; % probabilty of packet arrival downlink
par.packetsized = 100;   % average number of bits per packet downlink
par.packetlambdau = 0.1; % probabilty of packet arrival uplink
par.packetsizeu = 100;   % average number of bits per packet uplink
par.tau = 1e-7;          % average time delay in the channel
par.freq = 2e9;          % carrier frequency
par.RBbw = 2e5;          % RB bandwitdh 
par.nRB = 3;         % number of RBs over the bandwidth
par.RBduration = 1e-3;% RB duration [s] Same as par.slottime

par.product = 300; % 300 kbps in bps
par.packetlambdau_values = 0.01:0.01:1; % range of packetlambdau values to test
par.optimum_packetlambdau =[];
par.optimum_packetsizeu = [];
par.optimum_delay_probability = inf;
par.optimum_mean_delay = inf;
par.optimum_packetlambdau = 0;
par.optimum_packetsizeu = 0;
par.minimum_delay_probability = 0.05;
par.maximum_mean_delay = 10;


