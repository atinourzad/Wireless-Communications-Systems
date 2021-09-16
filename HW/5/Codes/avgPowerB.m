function [c,ceq] = avgPowerB(x)
    c = [];
    
    deltaG = 0.0001;
    
    g1 = x(1):deltaG:x(2);
    p1 = sum((10*log10(2)-g1) .* (1/sqrt(2*pi*64)) .* exp(-((g1-10).^2)/128).* deltaG);
    
    g2 = x(2):deltaG:x(3);
    p2 = sum((10*log10(6)-g2) .* (1/sqrt(2*pi*64)) .* exp(-((g2-10).^2)/128).* deltaG);
    
    g3 = x(3):deltaG:x(4);
    p3 = sum((10*log10(28)-g3) .* (1/sqrt(2*pi*64)) .* exp(-((g3-10).^2)/128).* deltaG);
    
    g4 = x(4):deltaG:x(4)+10;
    p4 = sum((10*log10(120)-g4) .* (1/sqrt(2*pi*64)) .* exp(-((g4-10).^2)/128) .* deltaG);
    
    ceq = p1 + p2 + p3 + p4;
end


