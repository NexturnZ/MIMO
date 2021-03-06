classdef QPSKTx < matlab.System  
%#codegen
% Generates the QPSK signal to be transmitted
    
%   Copyright 2012-2016 The MathWorks, Inc.
    
    properties (Nontunable)
        UpsamplingFactor = 4;
        MessageLength = 105;
        DataLength = 174;
        TransmitterFilterCoefficients = 1;
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];
        NumTxAntenna = 2;
        pLen = 8;
    end
    
    properties (Access=private)
        pBitGenerator
        pQPSKModulator 
        pTransmitterFilter
        pMIMOEncoder
        pBarkerCode = [+1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1]; % Bipolar Barker Code
        pModulatedHeader
        pPilots
    end
    
    methods
        function obj = QPSKTx(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj)
            obj.pBitGenerator = QPSKBitsGen(...
                'MessageLength', obj.MessageLength, ...
                'BernoulliLength', obj.DataLength-obj.MessageLength, ...
                'ScramblerBase', obj.ScramblerBase, ...
                'ScramblerPolynomial', obj.ScramblerPolynomial, ...
                'ScramblerInitialConditions', obj.ScramblerInitialConditions);
            obj.pQPSKModulator  = comm.QPSKModulator('BitInput',true, ...
                'PhaseOffset', pi/4);
            obj.pTransmitterFilter = dsp.FIRInterpolator(obj.UpsamplingFactor, ...
                obj.TransmitterFilterCoefficients);
            
            % MIMO
            obj.pMIMOEncoder = comm.OSTBCEncoder;
            obj.pMIMOEncoder.NumTransmitAntennas = obj.NumTxAntenna;
           
            obj.pModulatedHeader = sqrt(2)/2 * (-1-1i) * obj.pBarkerCode;
            
            W = hadamard(obj.pLen);
            obj.pPilots = W(:, 1:obj.NumTxAntenna);
        end
        
        function transmittedSignal = stepImpl(obj)
           
            [transmittedData, ~] = obj.pBitGenerator();                % Generates the data to be transmitted           
            modulatedData = obj.pQPSKModulator(transmittedData);       % Modulates the bits into QPSK symbols
            MIMOData = obj.pMIMOEncoder(modulatedData);
            
            transmittedSignal = [obj.pPilots; MIMOData];
        end
        
        function resetImpl(obj)
            reset(obj.pBitGenerator);
            reset(obj.pQPSKModulator );
            reset(obj.pTransmitterFilter);
        end
        
        function releaseImpl(obj)
            release(obj.pBitGenerator);
            release(obj.pQPSKModulator );
            release(obj.pTransmitterFilter);
        end
        
        function N = getNumInputsImpl(~)
            N = 0;
        end
    end
end

