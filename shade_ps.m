function shade_ps(T0,T0_ps,ylim)
h1 = area([0 0 T0_ps T0_ps], [ylim(1) ylim(2) ylim(2) ylim(1)], ...   % draw rectangles
                ylim(1), ...                % fill from bottom of plot instead of from 0
                'FaceColor', 'k', ...       % black
                'FaceAlpha', 0.12, ...      % transparent
                'LineStyle', 'none', ...    % no border
                'ShowBaseLine', 'off', ...  % no baseline
                'HandleVisibility', 'off'); % prevent this object from being added to legends
uistack(h1, 'bottom');           % put the rectangles behind all other plot elements

h2 = area([T0 T0 T0+T0_ps T0+T0_ps], [ylim(1) ylim(2) ylim(2) ylim(1)], ...   % draw rectangles
                ylim(1), ...                % fill from bottom of plot instead of from 0
                'FaceColor', 'k', ...       % black
                'FaceAlpha', 0.12, ...      % transparent
                'LineStyle', 'none', ...    % no border
                'ShowBaseLine', 'off', ...  % no baseline
                'HandleVisibility', 'off'); % prevent this object from being added to legends
uistack(h2, 'bottom');           % put the rectangles behind all other plot elements

h3 = area([2*T0 2*T0 2*T0+T0_ps 2*T0+T0_ps], [ylim(1) ylim(2) ylim(2) ylim(1)], ...   % draw rectangles
                ylim(1), ...                % fill from bottom of plot instead of from 0
                'FaceColor', 'k', ...       % black
                'FaceAlpha', 0.12, ...      % transparent
                'LineStyle', 'none', ...    % no border
                'ShowBaseLine', 'off', ...  % no baseline
                'HandleVisibility', 'off'); % prevent this object from being added to legends
uistack(h3, 'bottom');           % put the rectangles behind all other plot elements

h4 = area([3*T0 3*T0 3*T0+T0_ps 3*T0+T0_ps], [ylim(1) ylim(2) ylim(2) ylim(1)], ...   % draw rectangles
                ylim(1), ...                % fill from bottom of plot instead of from 0
                'FaceColor', 'k', ...       % black
                'FaceAlpha', 0.12, ...      % transparent
                'LineStyle', 'none', ...    % no border
                'ShowBaseLine', 'off', ...  % no baseline
                'HandleVisibility', 'off'); % prevent this object from being added to legends
uistack(h4, 'bottom');           % put the rectangles behind all other plot elements


end