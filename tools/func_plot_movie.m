function [] = func_plot_movie(start_step, movie_step, time_all, q1, q2, properties, line_prop)

l1=properties(3);
l2=properties(4);

value_x_test=l1*cos(q1)+l2*cos(q1+q2);
value_y_test=l1*sin(q1)+l2*sin(q1+q2);

value_x1=l1*cos(q1);
value_y1=l1*sin(q1);

if strcmp(line_prop, 'solid') == 1
    line_property = 1;
elseif strcmp(line_prop, 'dotted') == 1
    line_property = 2;
else
    disp('error: please input the line property!')
end

filename='./results/double_arm.gif';
figure()
for i=start_step:movie_step:time_all
    clf;
    hold on
    
    trace=1000;
    
    if i > trace * movie_step
        if line_property == 1
            plot(value_x_test(i-trace*movie_step:i), value_y_test(i-trace*movie_step:i), 'b', 'LineWidth', 1);
        elseif line_property == 2
            plot(value_x_test(i-trace*movie_step:i), value_y_test(i-trace*movie_step:i), 'b--', 'LineWidth', 1);
        end
        
    else
        if line_property == 1
            plot(value_x_test(1:i), value_y_test(1:i), 'b', 'LineWidth', 1);
        elseif line_property == 2
            plot(value_x_test(1:i), value_y_test(1:i), 'b--', 'LineWidth', 1);
        end
    end
    
    line([0,value_x1(i)], [0, value_y1(i)], 'LineWidth', 3);
    line([value_x1(i),value_x_test(i)], [value_y1(i),value_y_test(i)], 'LineWidth', 3)
    line([0, 0], [-1, 1], 'Color', 'black', 'LineStyle', '--')
    line([-1, 1], [0, 0], 'Color', 'black', 'LineStyle', '--')
    xlabel('x')
    ylabel('y')
    xlim([-1,1])
    ylim([-1,1])
    title(['time = ' ,num2str(i*0.01),'s'])
    drawnow
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if i==start_step
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end

end

