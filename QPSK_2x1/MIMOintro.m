clear all;  close all;
dbstop if error;
%%
% We start by defining some common simulation parameters
frmLen = 100;       % frame length
numPackets = 1000;  % number of packets
EbNo = 0:2:12;      % Eb/No varying to 20 dB
maxNumErrs = 300;       % maximum number of errors
maxNumPackets = 3000;   % maximum number of packets
pLen = 8;               % number of pilot symbols per frame
W = hadamard(pLen);

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
ostbcEnc = comm.OSTBCEncoder;
ostbcComb = comm.OSTBCCombiner;

awgn1Rx = comm.AWGNChannel(...
    'NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
    'SignalPower', 1);
awgn2Rx = clone(awgn1Rx);

% Create comm.ErrorRate calculator System objects to evaluate BER.
errorCalc1 = comm.ErrorRate;
errorCalc2 = comm.ErrorRate;
errorCalc3 = comm.ErrorRate;

%% QPSK
P = 4;				% modulation order
qpskMod = comm.QPSKModulator;
qpskDemod = comm.QPSKDemodulator('OutputDataType','double');

% Set up a figure for visualizing BER results
fig = figure; 
grid on;
ax = fig.CurrentAxes;
hold(ax,'on');

ax.YScale = 'log';
xlim(ax,[EbNo(1), EbNo(end)]);
ylim(ax,[1e-4 1]);
xlabel(ax,'Eb/No (dB)');
ylabel(ax,'BER'); 
fig.NumberTitle = 'off';
fig.Name = 'Orthogonal Space-Time Block Coding';
fig.Renderer = 'zbuffer';
title(ax,'Alamouti-coded 2x2 System');
set(fig,'DefaultLegendAutoUpdate','off');
fig.Position = figposition([41 50 25 30]);
%% Alamoti 2x2 
N = 2;              % maximum number of Tx antennas
M = 2;              % maximum number of Rx antennas

pilots = W(:, 1:N);     % orthogonal set per transmit antenna

chan = comm.MIMOChannel( ...
    'MaximumDopplerShift', 0, ...
    'SpatialCorrelationSpecification', 'None', ...
    'NumTransmitAntennas', N, ...
    'NumReceiveAntennas', M, ...
    'PathGainsOutputPort', true);

% object to M that is 2
release(ostbcEnc);
release(ostbcComb);
ostbcEnc.NumTransmitAntennas = N;
ostbcComb.NumReceiveAntennas = M;

% Release the hAWGN2Rx System object
release(awgn2Rx);

% Set the global random stream for repeatability
s = rng(55408);

% Pre-allocate variables for speed
HEst = zeros(frmLen, N, M);
ber_Estimate = zeros(3,length(EbNo));

% Loop over several EbNo points
for idx = 1:length(EbNo)
    reset(errorCalc1);
    reset(errorCalc2);
    awgn2Rx.EbNo = EbNo(idx); 

    % Loop till the number of errors exceed 'maxNumErrs'
    % or the maximum number of packets have been simulated
    while (ber_Estimate(2,idx) < maxNumErrs) && ...
          (ber_Estimate(3,idx)/frmLen < maxNumPackets)
        % Generate data vector per frame 
        data = randi([0 M-1], frmLen, 1);
        
        % Modulate data
        modData = qpskMod(data);           
        
        % Alamouti Space-Time Block Encoder
        encData = ostbcEnc(modData);

        % Prepend pilot symbols for each frame
        txSig = [pilots; encData];

        % Pass through the 2x2 channel        
        reset(chan);
        [chanOut, H] = chan(txSig);

        % Add AWGN
        rxSig = awgn2Rx(chanOut);

        % Channel Estimation
        %   For each link => N*M estimates
        HEst(1,:,:) = pilots(:,:).' * rxSig(1:pLen, :) / pLen;
        %   assume held constant for the whole frame
        HEst = HEst(ones(frmLen, 1), :, :);

        % Combiner using estimated channel
        decDataEst = ostbcComb(rxSig(pLen+1:end,:), HEst);
        
        % ML Detector (minimum Euclidean distance)
        demodEst   = qpskDemod(decDataEst);      % estimated  
        
        % Calculate and update BER for current EbNo value
        %   for estimated channel
        ber_Estimate(:,idx) = errorCalc1(data, demodEst);
    
    end % end of FOR loop for numPackets

    % Plot results
%     semilogy(ax,EbNo(1:idx), ber_Estimate(1,1:idx), 'ro');
%     drawnow;
end  % end of for loop for EbNo

% Perform curve fitting and replot the results
fitBEREst   = berfit(EbNo, ber_Estimate(1,:));
semilogy(ax,EbNo, fitBEREst, 'k.-');

% Restore default stream
rng(s)

%% Alamoti 2x1 
N = 2;              % maximum number of Tx antennas
M = 1;              % maximum number of Rx antennas

pilots = W(:, 1:N);     % orthogonal set per transmit antenna

chan = comm.MIMOChannel( ...
    'MaximumDopplerShift', 0, ...
    'SpatialCorrelationSpecification', 'None', ...
    'NumTransmitAntennas', N, ...
    'NumReceiveAntennas', M, ...
    'PathGainsOutputPort', true);

% object to M that is 2
release(ostbcEnc);
release(ostbcComb);
ostbcEnc.NumTransmitAntennas = N;
ostbcComb.NumReceiveAntennas = M;

% Release the hAWGN2Rx System object
release(awgn2Rx);

% Set the global random stream for repeatability
s = rng(55408);

% Pre-allocate variables for speed
HEst = zeros(frmLen, N, M);
ber_Estimate = zeros(3,length(EbNo));

% Loop over several EbNo points
for idx = 1:length(EbNo)
    reset(errorCalc1);
    reset(errorCalc2);
    awgn2Rx.EbNo = EbNo(idx); 

    % Loop till the number of errors exceed 'maxNumErrs'
    % or the maximum number of packets have been simulated
    while (ber_Estimate(2,idx) < maxNumErrs) && ...
          (ber_Estimate(3,idx)/frmLen < maxNumPackets)
        % Generate data vector per frame 
        data = randi([0 M-1], frmLen, 1);
        
        % Modulate data
        modData = qpskMod(data);           
        
        % Alamouti Space-Time Block Encoder
        encData = ostbcEnc(modData);

        % Prepend pilot symbols for each frame
        txSig = [pilots; encData];

        % Pass through the 2x2 channel        
        reset(chan);
        [chanOut, H] = chan(txSig);

        % Add AWGN
        rxSig = awgn2Rx(chanOut);

        % Channel Estimation
        %   For each link => N*M estimates
        HEst(1,:,:) = pilots(:,:).' * rxSig(1:pLen, :) / pLen;
        %   assume held constant for the whole frame
        HEst = HEst(ones(frmLen, 1), :, :);

        % Combiner using estimated channel
        decDataEst = ostbcComb(rxSig(pLen+1:end,:), HEst);
        
        % ML Detector (minimum Euclidean distance)
        demodEst   = qpskDemod(decDataEst);      % estimated  
        
        % Calculate and update BER for current EbNo value
        %   for estimated channel
        ber_Estimate(:,idx) = errorCalc1(data, demodEst);
    
    end % end of FOR loop for numPackets

    % Plot results
%     semilogy(ax,EbNo(1:idx), ber_Estimate(1,1:idx), 'k.');
%     drawnow;
end  % end of for loop for EbNo

% Perform curve fitting and replot the results
fitBEREst   = berfit(EbNo, ber_Estimate(1,:));
semilogy(ax,EbNo, fitBEREst, 'bd-');
% hold(ax,'off');

% Restore default stream
rng(s)

%% AWGN Channel

% Release the hAWGN2Rx System object
release(awgn2Rx);

% Set the global random stream for repeatability
s = rng(55408);

% Pre-allocate variables for speed
ber_Estimate = zeros(3,length(EbNo));

% Loop over several EbNo points
for idx = 1:length(EbNo)
    reset(errorCalc1);
    reset(errorCalc2);
    awgn2Rx.EbNo = EbNo(idx); 

    % Loop till the number of errors exceed 'maxNumErrs'
    % or the maximum number of packets have been simulated
    while (ber_Estimate(2,idx) < maxNumErrs) && ...
          (ber_Estimate(3,idx)/frmLen < maxNumPackets)
        % Generate data vector per frame 
        data = randi([0 M-1], frmLen, 1);
        
        % Modulate data
        modData = qpskMod(data);           

        % Add AWGN
        rxSig = awgn2Rx(modData);
        
        % ML Detector (minimum Euclidean distance)
        demodEst   = qpskDemod(rxSig);      % estimated  
        
        % Calculate and update BER for current EbNo value
        %   for estimated channel
        ber_Estimate(:,idx) = errorCalc1(data, demodEst);
    
    end % end of FOR loop for numPackets

    % Plot results
%     semilogy(ax,EbNo(1:idx), ber_Estimate(1,1:idx), 'k.');
%     drawnow;
end  % end of for loop for EbNo

% Perform curve fitting and replot the results
fitBEREst   = berfit(EbNo, ber_Estimate(1,:));
semilogy(ax,EbNo, fitBEREst, 'r*-');

% Restore default stream
rng(s)