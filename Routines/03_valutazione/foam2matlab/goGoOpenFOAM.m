function [result,log]=goGoOpenFOAM(command)

[result,log]=unix(['export LD_LIBRARY_PATH=""; . /opt/openfoam5/etc/bashrc; ' command]);

end