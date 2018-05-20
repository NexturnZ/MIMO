close all; clear all;
dbstop if error;
%% Initialization

prmQPSKTxRx = commqpsktxrx_init; % QPSK system parameters 
prmQPSKTxRx.EbNo = 3;
useScopes = false; % true if scopes are to be used
printReceivedData = true; %true if the received data is to be printed
compileIt = false; % true if code is to be compiled
useCodegen = false; % true to run the generated mex file
%% Execution and Results

if compileIt
    codegen -report runQPSKSystemUnderTest.m -args {coder.Constant(prmQPSKTxRx),coder.Constant(useScopes),coder.Constant(printReceivedData)} %#ok
end
if useCodegen
    BER = runQPSKSystemUnderTest_mex(prmQPSKTxRx, useScopes, printReceivedData);  
else
    BER = runQPSKSystemUnderTest(prmQPSKTxRx, useScopes, printReceivedData);
end
fprintf('Error rate = %f.\n',BER(1));
fprintf('Number of detected errors = %d.\n',BER(2));
fprintf('Total number of compared samples = %d.\n',BER(3));
displayEndOfDemoMessage(mfilename)
