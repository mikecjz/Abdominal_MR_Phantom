function time = TRNumToTime(TRNum) %in msec. For SLIDER
    val = floor((TRNum-1)./400);
    
    remainder = mod(TRNum-1,400)+1;
    time = 0;
    if remainder>=1 && remainder<=300
        time = remainder*6+3266*val+100*(val+1);
    elseif remainder>=301 && remainder<=325
        time = remainder*6+3266*val+100*(val+1)+700;
    elseif remainder>=326
        time = remainder*6+3266*val+100*(val+1)+766;
    end
    
    
end