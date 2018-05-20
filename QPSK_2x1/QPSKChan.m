classdef QPSKChan < matlab.System 
%#codegen
    
%   Copyright 2012-2016 The MathWorks, Inc.
   
    
    properties (Nontunable)
        DelayType = 'Triangle';
        RaisedCosineFilterSpan = 10;
        PhaseOffset = 47;
        SignalPower = 0.25;
        FrameSize = 100;
        UpsamplingFactor = 4;
        EbNo = 7;
        BitsPerSymbol = 2;
        FrequencyOffset = 5000;
        SampleRate = 200000;
    end
    
    properties (Access=private)
        pPhaseFreqOffset
        pVariableTimeDelay
        pAWGNChannel
        pMIMOChannel
    end
    
    properties (Constant, Access=private)
        pDelayStepSize = 0.05;
        pDelayMaximum = 8;
        pDelayMinimum = 0.1;
    end
    
    methods
        function obj = QPSKChan(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~, ~)
            obj.pPhaseFreqOffset = comm.PhaseFrequencyOffset(...
                'PhaseOffset', obj.PhaseOffset, ...
                'FrequencyOffset', obj.FrequencyOffset, ...
                'SampleRate',obj.SampleRate);
            obj.pVariableTimeDelay = dsp.VariableFractionalDelay(...
                'MaximumDelay', obj.FrameSize*obj.UpsamplingFactor);
            obj.pAWGNChannel = comm.AWGNChannel('EbNo', obj.EbNo, ...
                'BitsPerSymbol', obj.BitsPerSymbol, ...
                'SignalPower', obj.SignalPower, ...
                'SamplesPerSymbol', obj.UpsamplingFactor);
            obj.pMIMOChannel = comm.MIMOChannel( ...
                'MaximumDopplerShift', 0, ...
                'SpatialCorrelationSpecification', 'None', ...
                'NumTransmitAntennas', 2, ...
                'NumReceiveAntennas', 2, ...
                'PathGainsOutputPort', true);
         end
        
        
        function [corruptSignal, H] = stepImpl(obj, TxSignal, count)
            
            
            % MIMO channel
            [MIMOsig,H] = obj.pMIMOChannel(TxSignal);
            corruptSignal = obj.pAWGNChannel(MIMOsig);
            
        end
        
        function resetImpl(obj)
            reset(obj.pPhaseFreqOffset);
            reset(obj.pVariableTimeDelay);            
            reset(obj.pAWGNChannel);
        end
        
        function releaseImpl(obj)
            release(obj.pPhaseFreqOffset);
            release(obj.pVariableTimeDelay);            
            release(obj.pAWGNChannel);            
        end
        
        function N = getNumInputsImpl(~)
            N = 2; 
        end
    end
end

