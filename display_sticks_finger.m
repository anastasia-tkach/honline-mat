function [] = display_kalman_2D(segments, data_points, model_points, iter)

%if iter == 1 || iter == 2
    %if iter == 1, shift = 0.185; else shift = 0.507; end
    %figure('units', 'normalized', 'outerposition', [shift, 0.4, 0.33, 0.47]); axis off; axis equal; hold on;
%else
    %clf; hold on; axis off; axis equal;
%end
%set(gcf,'color','w');
%set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');

mylines(data_points, model_points, [0.85, 0.85, 0.85]);
mypoints({segments{1}.global(1:2, 4)}, [0.0, 0.3, 0.7], 35);
for i = 1:length(segments)
    if segments{i}.parent_id == 0, continue; end
    a = segments{i}.global(1:2, 4);
    b = segments{segments{i}.parent_id}.global(1:2, 4);
    myline(a, b, [1.0, 0.6, 0.5], 5);
    mypoints({a}, [1.0, 0.6, 0.5], 30);
    mypoints({b}, [1.0, 0.6, 0.5], 30);
end
mypoints(data_points, [0.3, 0.6, 0.8], 20);
xlim([-10, 10]); ylim([-5, 10]);
%xlim([-9, 2]); ylim([-1, 10]);
drawnow;