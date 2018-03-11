function [style] = plot_style(i)

%            b     blue          .     point              -     solid
%            g     green         o     circle             :     dotted
%            r     red           x     x-mark             -.    dashdot 
%            c     cyan          +     plus               --    dashed   
%            m     magenta       *     star             (none)  no line
%            y     yellow        s     square
%            k     black         d     diamond
%            w     white         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram


c = 'bgrcmk';
s = 'ox+*sdv^<>ph';

i_max = size(c,2)*size(s,2);

if i > i_max
    error('casi maggiori delle possibili rappresentazioni \n')
end

si = ceil(i/size(c,2));
ci = i- (si-1)*size(c,2);

style = [c(ci),s(si),'-'];

end