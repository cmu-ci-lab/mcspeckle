function p = evalampfunc_amplitude(ampfunc,angs,dim)
    % provided amplitude function rotation amplitude
    N = length(ampfunc.evalAmp);
    
    if(dim == 2)
        p = ampfunc.evalAmp(round(angs*(N-1)/(2*pi))+1);
    else
        p = ampfunc.evalAmp(round(angs*(N-1)/(pi))+1);
    end
end