function [] = close_all_but_gui()
%close_all_but_gui: This function closes all figures except the gui
%   Detailed explanation goes here
  fig_h = permute( findobj( 0, 'Type', 'Figure' ), [2,1] );
        for fh = fig_h
            uih = findobj( fh, 'Type', 'uicontrol' );
            if isempty( uih )
                delete( fh );
            end
        end

end

