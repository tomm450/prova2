function result=goGoOpenFOAM(command)
result=unix(['export LD_LIBRARY_PATH=""; . /opt/openfoam5/etc/bashrc; ' command]);
end