function nautilus(link)
% apre cartella da linea di comando

[~,~] = system(sprintf...
    ('export LD_LIBRARY_PATH=""; ./;  nautilus --display=DISPLAY %s >/dev/null'...
    ,link));

end

