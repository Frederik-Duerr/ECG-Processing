function [filtered_normalized_meanRemoved_signal] = normalisation_function(filtered_signal)
%% Applying the Normalisation to the Data Set
% Loading the filtered data
% Formular of from the Paper: x(n) = x(n)/ sum[n = 0 - N-1](x^2)`; N =
% length of the Signal. Afterwards the mean was removed



 filtered_normalized_meanRemoved_signal = zeros(size(filtered_signal));

for i = 1:size(filtered_signal,2)
    filtered_normalized_signal(:,i) = filtered_signal(:,i)/sum(filtered_signal(:,i).^2); % calculates the energy of the Signal in column i
    % The division makes Energy = 1
    
    % Remove the mean from the normalized signal = removes the DC-Offset
    filtered_normalized_meanRemoved_signal(:,i) = filtered_normalized_signal(:,i) - mean(filtered_normalized_signal(:,i));
    
end

end