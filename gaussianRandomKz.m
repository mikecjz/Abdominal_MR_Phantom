function kz = gaussianRandomKz(npar,nt,sigma,varargin)
if ~isempty(varargin)
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'seed'
                seed = varargin{i+1};
                
            otherwise % skip it
                error('Unsuportted argument')
        end
    end
    
    rng(seed)
    parSelect = round(randn(1,nt)*sigma+(npar/2+1));
    
    outofboundidx = find((parSelect<1|parSelect>npar));
    
    preLength = nt;
    additionalNumber = length(outofboundidx);
    
    while additionalNumber>0
        rng(seed)
        additionalSlect = round(randn(1,preLength+additionalNumber)*sigma+(npar/2+1));
        additionalSlect=additionalSlect(1,preLength+1:preLength+additionalNumber);
        
        preLength = preLength+additionalNumber;
        
        parSelect(1,outofboundidx) = additionalSlect;
        
        outofboundidx = find((parSelect<1|parSelect>npar));
        additionalNumber = length(outofboundidx);
        
    end
else
    parSelect = zeros(1,nt);
    for i = 1:nt
        newrand = round(randn(1)*sigma+(npar/2+1));
        while(newrand<1||newrand>npar)
            newrand = round(randn(1)*sigma+(npar/2+1));
        end
        
        parSelect(i) = newrand;
    end
end

kz = parSelect;


end