    function result=goGoGmsh(command)
        %fullcom = sprintf('~/Documents/gmsh %s',command);
        fullcom = sprintf('gmsh %s',command);
        result=unix(['export LD_LIBRARY_PATH=""; . /usr/lib/x86_64-linux-gnu/; ' fullcom]);
    end