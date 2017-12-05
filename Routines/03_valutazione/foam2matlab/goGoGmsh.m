    function result=goGoGmsh(command,gcom,pref)
        if nargin == 2
            pref = '';
        end
        %fullcom = sprintf('~/Documents/gmsh %s',command);
        fullcom = sprintf('%s %s %s',pref,gcom,command);
        
        result=unix(['export LD_LIBRARY_PATH=""; . /usr/lib/x86_64-linux-gnu/; ' fullcom]);
    
    end