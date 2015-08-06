function [ output ] = normax( input )
%normax( input ) return the normalized form of signal according to its max
%   input - signal
%   output - normalized signal

if any(input)
    output = input/max(abs(input(:)));
else 
    output = input;
end

end



