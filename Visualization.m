function Visualization(db, r, p, kappa, c, del, time, FrameSkip, h, eta1fit, eta2fit, eta1mod, eta2mod, eta1modapprox, eta2modapprox, eta1free, eta2free, ...
                       c_fit, l_fit, cmod, c2,lmod,l2)

    [a, M] = size(r);
    N = (M-1)/2;
    t_num=time/(FrameSkip*h)+1;
    timegrid = (0:t_num-1) * FrameSkip * h;

    frameRate = 10; 
    writerObj = VideoWriter(fullfile('Videos', 'output.mp4'), 'MPEG-4');
    writerObj.FrameRate = frameRate;
    open(writerObj);

    % Create an off-screen (invisible) figure
    fig = figure('Visible','off', 'Position',[100,100,1120,840]);

    for i = 1:a
        %--- Subplot #1
        ax1 = subplot(3,2,1, 'Parent', fig);
        plot(ax1, (-N:N), r(i,:), 'b.'); 
        hold(ax1, 'on');
        plot(ax1, (-N:N), profile(db, 'r', c_fit(i), l_fit(i), (-N:N)),'r');
        title(ax1, '$r(j)$');
        ylim(ax1, [min(r(:))-del/20, max(r(:))+del/20]);
        xlim(ax1, [-(N/10)+timegrid(i)*c, ...
                     (N/10)+timegrid(i)*c]);
        hold(ax1, 'off');

        %--- Subplot #2
        ax2 = subplot(3,2,2, 'Parent', fig);
        plot(ax2, (-N:N), p(i,:), 'b.');
        hold(ax2, 'on');
        plot(ax2, (-N:N), profile(db, 'p', c_fit(i), l_fit(i), (-N:N)),'r');
        title(ax2, '$p(j)$');
        ylim(ax2, [min(p(:))-del/20, max(p(:))+del/20]);
        xlim(ax2, [-(N/10)+timegrid(i)*c, ...
                     (N/10)+timegrid(i)*c]);
        hold(ax2, 'off');

        %--- Subplot #3
        ax3 = subplot(3,2,3, 'Parent', fig);
        plot(ax3, (-N:N), eta1mod(i,:), 'r.');
        title(ax3, '$\eta_1(j)$');
        hold(ax3, 'on');
        plot(ax3, (-N:N), eta1fit(i,:), 'b.');
        hold(ax3, 'on');
        plot(ax3, (-N:N), eta1modapprox(i,:), 'g.');
        % hold(ax3, 'on');
        % plot(ax3, (-N:N), eta1free(i,:),'m.');
        ylim(ax3, [min(eta1fit(:))+min(eta1mod(:)), max(eta1fit(:))+max(eta1mod(:))]);
        xlim(ax3, [-(N/10)+timegrid(i)*c, ...
                     (N/10)+timegrid(i)*c]);
        legend(ax3, 'mod', 'fit','2');

        %--- Subplot #4
        ax4 = subplot(3,2,4, 'Parent', fig);
        plot(ax4, (-N:N), eta2mod(i,:), 'r.');
        title(ax4, '$\eta_2(j)$');
        hold(ax4, 'on');
        plot(ax4, (-N:N), eta2fit(i,:), 'b.');
        hold(ax4, 'on');
        plot(ax4, (-N:N), eta2modapprox(i,:), 'g.');
        % hold(ax4, 'on');
        % plot(ax4, (-N:N), eta2free(i,:), 'm.');
        ylim(ax4, [min(eta2fit(:))+min(eta2mod(:)), max(eta2fit(:))+max(eta2mod(:))]);
        xlim(ax4, [-(N/10)+timegrid(i)*c, ...
                     (N/10)+timegrid(i)*c]);
        legend(ax4, 'mod', 'fit','2');

        %--- Subplot #5
        ax5 = subplot(3,2,5, 'Parent', fig);
        plot(ax5, timegrid(1:i), c_fit(1:i),'b');
        hold(ax5, 'on');
        plot(ax5, timegrid(1:i), 1+cmod(1:i)/24,'r');
        hold(ax5, 'on');
        plot(ax5, timegrid(1:i), 1+c2(1:i)/24,'g');
        legend(ax5, 'fit', 'mod','2');
        xlabel(ax5, '$t$');
        title('$c$')
        %ylim(ax5, [0.7*c_fit(1), 1.3*c_fit(1)]);
        xlim(ax5, [0, time]);
        hold(ax5, 'off');


                %--- Subplot #6
        ax6 = subplot(3,2,6, 'Parent', fig);
        plot(ax6, timegrid(1:i),l_fit(1:i)-timegrid(1:i)*c,'b');
        hold(ax6, 'on');
        plot(ax6, timegrid(1:i),lmod(1:i)-timegrid(1:i)*c,'r');
        hold(ax6, 'on');
        plot(ax6, timegrid(1:i),l2(1:i)-timegrid(1:i)*c,'g');
        legend(ax6, 'fit', 'mod','2');
        xlabel(ax6, '$t$');
        title('$l-c_*$')
        %ylim(ax6, [0.7*l_values(1), 1.3*c_values(1)]);
        xlim(ax6, [0, time]);
        hold(ax6, 'off');

       

        drawnow;                 % Ensure the figure is updated
        frame = getframe(fig);   % Capture the *invisible* figure
        writeVideo(writerObj, frame);

        clf(fig);  % Clear the figure for the next iteration
    end

    close(writerObj);
    close(fig);

end
