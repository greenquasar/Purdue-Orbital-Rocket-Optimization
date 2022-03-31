function RocketDrawer(stage_length, stage_width, inner_width)
    
    yaxis = [0 stage_length stage_length 0 0];
    xaxis = [0 0 stage_width stage_width 0];
    i = 1;
    for theta = 0:0.01:2*pi
        webX(i) = stage_width*2 + inner_width/2 * cos(theta);
        webY(i) = stage_width/2 + inner_width/2 * sin(theta);
        widthX(i) = stage_width*2 + stage_width/2 * cos(theta);
        widthY(i) = stage_width/2 + stage_width/2 * sin(theta);
        i = i + 1;
    end
    figure(2)
    plot(xaxis, yaxis);
    hold on
    plot(webX, webY);
    plot(widthX, widthY);
    legend("Body", "Web", "Body");
    axis padded
    axis equal
    movegui('northeast');
    hold off
end