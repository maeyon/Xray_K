% 16/11/3 L計算
% 32x32マス、(0,0)-(32,32)
% Xs(t):t枚目のX線源のX座標
% Ys(t):t枚目のX線源のY座標
% Xd(t,i):t枚目i番目の検出器のX座標
% Yd(t,i):t枚目i番目の検出器のY座標
% X線源とi番目の検出器を結ぶ直線は下記の連立方程式より解ける
% Ys = aXs + b
% Ydi = aXdi + b
% よって、a = (Ys-Ydi)/(Xs-Xdi), b = Ys-(Ys-Ydi)/(Xs0-Xdi)*Xs
%
%%%%%%
%%%%%%
% 事前指定
i = 1;
XMAX = 32;
YMAX = 32;
DMAX = 32;

Xs0 = XMAX/2;   %回転角ゼロのX線源座標
Ys0 = -5;       %回転角ゼロのX線源座標
TMAX = 30;     %投影枚数
XCENTER = XMAX/2;      %回転中心
YCENTER = YMAX/2;      %回転中心

%回転角ゼロの検出器座標
Xd0 = XMAX/2;
Yd0 = YMAX-Ys0;

%回転角を考慮したt枚目の投影のX線源の座標(ラジアン注意)
t = 1:TMAX;
Xs(t) = (Xs0-XCENTER)*cos(t*2*pi/TMAX)-(Ys0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER;
Xs = transpose(repmat(Xs,DMAX,1));
Ys(t) = (Xs0-XCENTER)*sin(t*2*pi/TMAX)+(Ys0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER;
Ys = transpose(repmat(Ys,DMAX,1));


%回転角を考慮したt枚目,i番目の投影の検出器の座標(ラジアン注意)
t = 1:TMAX;
    for i = 1:DMAX
        Xd(t,i) = (Xd0-XCENTER)*cos(t*2*pi/TMAX)-(Yd0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER-(DMAX/2-i+0.5)*cos(t*2*pi/TMAX);
        Yd(t,i) = (Xd0-XCENTER)*sin(t*2*pi/TMAX)+(Yd0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER+(DMAX/2-i+0.5)*sin(t*2*pi/TMAX);
    end

%%%%%%
%%%%%%
%1x1の枠の中を線が通るか？
%各辺ごとにif文をつくる
%傾きで分岐

x = 0;
y = 0;
j = 1;

flg_in_up = 0;
flg_in_left = 0;
flg_in_down = 0;
flg_in_right = 0;

flg_out_up = 0;
flg_out_left = 0;
flg_out_down = 0;
flg_out_right = 0;

a = (Ys-Yd)./(Xs-Xd);
b = Ys-(Ys-Yd)./(Xs-Xd).*Xs;

for t = 1:TMAX
    for i = 1:DMAX
        for x = 0:XMAX-1
            for y = 0:YMAX-1
                
                %交点の座標
                xc_up = (y+1-b(t,i))/a(t,i);
                yc_left = a(t,i)*x+b(t,i);
                xc_down = (y-b(t,i))/a(t,i);
                yc_right = a(t,i)*(x+1)+b(t,i);

                %入口
                if Ys(t,i) <= YCENTER %X線源が検出器より下にあるとき
                    if a(t,i) <= 0 %傾きが負のとき、下辺/右辺から入って上辺/左辺から出る
                        if (xc_down >= x) && (xc_down <= x+1)
                            flg_in_down = 1;
                            flg_in_right = 0;
                        elseif (yc_right >= y) && (yc_right <= y+1)
                            flg_in_down = 0;
                            flg_in_right = 1;
                        else
                            flg_in_down = 0;
                            flg_in_right = 0;
                        end
                    else %傾きが正のとき左辺/下辺から入って上辺/右辺から出る
                        if (xc_down >= x) && (xc_down <= x+1)
                            flg_in_down = 1;
                            flg_in_left = 0;
                        elseif (yc_left >= y) && (yc_left <= y+1)
                            flg_in_down = 0;
                            flg_in_left = 1;
                        else
                            flg_in_down = 0;
                            flg_in_left = 0;
                        end
                    end
                else %X線源が検出器より上にあるとき
                    if a(t,i) <= 0 %傾きが負のとき、上辺/左辺から入って下辺/右辺から出る
                        if (xc_up >= x) && (xc_up <= x+1)
                            flg_in_up = 1;
                            flg_in_left = 0;
                        elseif (yc_left >= y) && (yc_left <= y+1)
                            flg_in_up = 0;
                            flg_in_left = 1;
                        else
                            flg_in_up = 0;
                            flg_in_left = 0;
                        end
                    else %傾きが正のとき上辺/右辺から入って下辺/左辺から出る
                        if (xc_up >= x) && (xc_up <= x+1)
                            flg_in_up = 1;
                            flg_in_right = 0;
                        elseif (yc_right >= y) && (yc_right <= y+1)
                            flg_in_up = 0;
                            flg_in_right = 1;
                        else
                            flg_in_up = 0;
                            flg_in_right = 0;
                        end
                    end
                end
                
                %出口
                if Ys(t,i) <= YCENTER %X線源が検出器より下にあるとき
                    if a(t,i) <= 0 % 傾きが負のとき、下辺/右辺から入って上辺/左辺から出る
                        if (xc_up >= x) && (xc_up <= x+1)
                            flg_out_up = 1;
                            flg_out_left = 0;
                        elseif (yc_left >= y) && (yc_left <= y+1)
                            flg_out_up = 0;
                            flg_out_left = 1;
                        else
                            flg_out_up = 0;
                            flg_out_left = 0;
                        end
                    else %傾きが正のとき左辺/下辺から入って上辺/右辺から出る
                        if (xc_up >= x) && (xc_up <= x+1)
                            flg_out_up = 1;
                            flg_out_right = 0;
                        elseif (yc_right >= y) && (yc_right <= y+1)
                            flg_out_up = 0;
                            flg_out_right = 1;
                        else
                            flg_out_up = 0;
                            flg_out_right = 0;
                        end
                    end
                else %X線源が検出器より上にあるとき
                    if a(t,i) <= 0 %傾きが負のとき、上辺/左辺から入って下辺/右辺から出る
                        if (xc_down >= x) && (xc_down <= x+1)
                            flg_out_down = 1;
                            flg_out_right = 0;
                        elseif (yc_right >= y) && (yc_right <= y+1)
                            flg_out_down = 0;
                            flg_out_right = 1;
                        else
                            flg_out_down = 0;
                            flg_out_right = 0;
                        end
                    else %傾きが正のとき上辺/右辺から入って下辺/左辺から出る
                        if (xc_down >= x) && (xc_down <= x+1)
                            flg_out_down = 1;
                            flg_out_left = 0;
                        elseif (yc_left >= y) && (yc_left <= y+1)
                            flg_out_down = 0;
                            flg_out_left = 1;
                        else
                            flg_out_down = 0;
                            flg_out_left = 0;
                        end
                    end
                end

                %長さの計算
                %入り口と出口によって場合分け
                
                %下辺->上辺
                if (flg_in_down == 1) && (flg_out_up == 1)
                    l(j) = sqrt((xc_down-xc_up)^2 + 1);
                    j = j + 1;
                %下辺->左辺
                elseif (flg_in_down == 1) && (flg_out_left == 1)
                    l(j) = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %下辺->右辺
                elseif (flg_in_down == 1) && (flg_out_right == 1)
                    l(j) = sqrt((yc_right-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %右辺->上辺
                elseif (flg_in_right == 1) && (flg_out_up == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                    j = j + 1;
                %右辺->左辺
                elseif (flg_in_right == 1) && (flg_out_left == 1)
                    l(j) = sqrt((yc_right-yc_left)^2 + 1);
                    j = j + 1;
                %右辺->下辺
                elseif (flg_in_right == 1) && (flg_out_down == 1)
                    l(j) = sqrt((yc_right-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %左辺->上辺
                elseif (flg_in_left == 1) && (flg_out_up == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_left-(y+1))^2);
                    j = j + 1;
                %左辺->右辺
                elseif (flg_in_left == 1) && (flg_out_right == 1)
                    l(j) = sqrt((yc_left-yc_right)^2 + 1);
                    j = j + 1;
                %左辺->下辺
                elseif (flg_in_left == 1) && (flg_out_down == 1)
                    l(j) = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %上辺->下辺
                elseif (flg_in_up == 1) && (flg_out_down == 1)
                    l(j) = sqrt((xc_down-xc_up)^2 + 1);
                    j = j + 1;
                %上辺->左辺
                elseif (flg_in_up == 1) && (flg_out_left == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_left-(y+1))^2);
                    j = j + 1;
                %上辺->右辺
                elseif (flg_in_up == 1) && (flg_out_right == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                    j = j + 1;
                else
                    l(j) = 0;
                    j = j + 1;
                end
            end
        end
    end
end

L = reshape(l,[32768,30]);