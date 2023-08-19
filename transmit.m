function [cdl, idl, sirdl]=transmit(isdeslink, pdb, gdb, noisedb, dt)
% DESCRIPTION
%  all variables in dB
% INPUT
%  isdeslink -- binary matrix indicating desired links
%  pdb       -- transmitted power
%  gdb       -- gain matrix
%  noisedb   -- noise floor [dBm]
%  dt        -- dimension spanning transmitters
% OUTPUT
%  cdl       -- received carrier power down link [dBm]
%  idl       -- interference power down link  [dBm]
%  sirdl     -- signal to interference ratio down link [dB]
% TRY
%

% by Magnus Almgren 080121

[p,g,noise] = db2lin(pdb,gdb,noisedb); % convert to linear measures
if dt == 2 % Down link. Move power over to second dimension
    p = sum(mprod(p,isdeslink,1)); % calculate the base station power 
                                   % in this case we have one user assigned
                                   % to the BS, so total power is used to
                                   % transmit to the allocated user
                                   % this works when BS transmits with different power 
                                   % to multiple users 
end

% find the interference + noise( I+N )  in DL
% step 1: consider all the users not scheduled ~isdeslink
% step 2: multiply by TX power and gain 
% step 3: Sum of all powers from each BS
% step 4: Add noise
%
% find received signal ( S )
% step 1: consider all the user scheduled isdeslink
% step 2: multiply by TX power and gain 
%
% find the ratio of above two S/(I+N) i.e., SINR
%
printerf = sum(mprod(p,g,~isdeslink),dt) + noise; % non desired power at the receiver
prdes    = sum(mprod(p,g, isdeslink),dt); % desired power at the receiver

% for the UL:  we only consider all the non-desiredpower multiplied by the
% isdeslink i.e.,  isdeslinkum

if dt == 1 % Up link. Move over to first dimension
    printerf = sum(mprod(printerf,isdeslink),2);
    prdes =    sum(mprod(prdes,   isdeslink),2);
end

snr = mdiv(prdes,printerf);
snr(mand(printerf==0,prdes==0)) = 0; % make 0/0 go from Nan to 0

[cdl, idl, sirdl] = lin2db(prdes,printerf,snr); % convert back to dB

