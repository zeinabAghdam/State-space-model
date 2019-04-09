function [centeredData] = locallyCenterData(vData, window)
    assert(window >= 2);
    halfWindow = floor(window/2);
    
    outI = 1;
    centeredData = zeros(1, length(vData) - halfWindow*2);
    
    for i = (halfWindow+1):(length(vData) - halfWindow)
        windowBegin = i - halfWindow;
        windowEnd   = i + halfWindow;
        
        windowMean = circ_mean(vData(windowBegin:windowEnd)');
        
        centeredData(outI) = vData(i) - windowMean;
        centeredData(outI) = mod(centeredData(outI) + pi, 2*pi) - pi;
        outI = outI + 1;
    end
end