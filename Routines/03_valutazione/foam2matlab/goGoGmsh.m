    function result=goGoGmsh(command,gcom)
        
        %fullcom = sprintf('~/Documents/gmsh %s',command);
        fullcom = sprintf('%s %s',gcom,command);
        
        result=unix(['export LD_LIBRARY_PATH=""; . /usr/lib/x86_64-linux-gnu/; ' fullcom]);
    
    end