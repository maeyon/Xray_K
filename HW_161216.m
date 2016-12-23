% 16/12/14 L計算
% [16/12/9の宿題をもとに変更点]
% (済)(1)l(j)->lに変更
% (2)傾き無限大の対策
% (3)投影像の表示
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
XMAX = 3;  %X画素数
YMAX = 3;  %Y画素数
DMAX = 2; %検出器の数
DWIDTH = 1;

Xs0 = XMAX/2;   %回転角ゼロのX線源のX座標
Ys0 = -5;       %回転角ゼロのX線源のY座標
TMAX = 3;     %投影枚数
XCENTER = XMAX/2;      %回転中心
YCENTER = YMAX/2;      %回転中心

%Xd0 = XMAX/2;   %回転角ゼロの検出器のX座標
Xd0 = [0:XMAX/(DMAX-1):XMAX];   %回転角ゼロの検出器のX座標
Yd0 = (YMAX-Ys0)*ones(1,DMAX);  %回転角ゼロの検出器のY座標
Genten_Hosei_D = [Xd0 - XCENTER; Yd0 - YCENTER];    %回転角ゼロの検出器の座標(原点に移動)

%回転角を考慮したt枚目の投影のX線源の座標(ラジアン注意)
clear Xs Ys Xd Yd
for t = 1:TMAX
    % 回転行列R　の　作成
    Rt = [cos((t-1)*2*pi/TMAX) -sin((t-1)*2*pi/TMAX);sin((t-1)*2*pi/TMAX) cos((t-1)*2*pi/TMAX)];
    %　原点補正後の検出器座標の回転
    temp = Rt*Genten_Hosei_D;
    %　もとの座標位置に復元
    Xd(t,:) = temp(1,:) + XCENTER;
    Yd(t,:) = temp(2,:) + YCENTER;
    
    %　原点補正後のX線源座標の回転
    temp = Rt*[Xs0-XCENTER;Ys0-YCENTER];
    % もとの座標位置に復元
    Xs(t) = temp(1) + XCENTER; 
    Ys(t) = temp(2) + YCENTER;
end

Xs = repmat(Xs(:),1,DMAX);
Ys = repmat(Ys(:),1,DMAX);

% tic;
% % 回転行列R　の　作成
% Rall = [cos([1:TMAX]'*2*pi/TMAX) -sin([1:TMAX]'*2*pi/TMAX);sin([1:TMAX]'*2*pi/TMAX) cos([1:TMAX]'*2*pi/TMAX)];
% %　原点補正後の座標の回転
% temp = Rall*Genten_Hosei_D;
% %　原点補正後の座標を回転したものをもとの座標位置に復元
% Xd2 = temp(1:TMAX,:) + XCENTER; 
% Yd2 = temp(TMAX+1:end,:) + YCENTER;
% toc;
% 
% t = 1:TMAX;
% Xs(t) = (Xs0-XCENTER)*cos(t*2*pi/TMAX)-(Ys0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER;
% Xs = transpose(repmat(Xs,DMAX,1));
% Ys(t) = (Xs0-XCENTER)*sin(t*2*pi/TMAX)+(Ys0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER;
% Ys = transpose(repmat(Ys,DMAX,1));
% 
% 
% 
% %回転角を考慮したt枚目,i番目の投影の検出器の座標(ラジアン注意)
% %ここ修正
% t = 1:TMAX;
%     for i = 1:DMAX
%         Xd(t,i) = (Xd0-XCENTER)*cos(t*2*pi/TMAX)-(Yd0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER-(DMAX/2-i+0.5)*cos(t*2*pi/TMAX);
%         Yd(t,i) = (Xd0-XCENTER)*sin(t*2*pi/TMAX)+(Yd0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER+(DMAX/2-i+0.5)*sin(t*2*pi/TMAX);
%     end


    
%%%%%%
%%%%%%
%1x1の枠の中を線が通るか？
%各辺ごとにif文をつくる
%傾きで分岐

x = 0;
y = 0;
j = 1;

a = (Ys-Yd)./(Xs-Xd);   %傾き
b = Ys-(Ys-Yd)./(Xs-Xd).*Xs;    %切片

c = (Xs-Xd)./(Ys-Yd);
d = Xs-(Xs-Xd)./(Ys-Yd).*Ys;

clear TDs Ps ls
k=1;
for t = 1:TMAX
    for i = 1:DMAX
        for x = 0:XMAX-1
            for y = 0:YMAX-1
                flg_in_up = 0;
                flg_in_left = 0;
                flg_in_down = 0;
                flg_in_right = 0;
                
                flg_out_up = 0;
                flg_out_left = 0;
                flg_out_down = 0;
                flg_out_right = 0;

                %交点の座標
                if abs(a(t,i)) > 1
                    xc_up = c(t,i)*(y+1)+d(t,i);
                    yc_left = (x-d(t,i))/c(t,i);
                    xc_down = c(t,i)*y+d(t,i);
                    yc_right = (x+1-d(t,i))/c(t,i);
                else
                    xc_up = (y+1-b(t,i))/a(t,i);
                    yc_left = a(t,i)*x+b(t,i);
                    xc_down = (y-b(t,i))/a(t,i);
                    yc_right = a(t,i)*(x+1)+b(t,i);
                end

                %入口
                if Ys(t,i) <= YCENTER %X線源が検出器より下にあるとき
                    if a(t,i) <= 0 %傾きが負のとき、下辺もしくは右辺から入る
                        if (xc_down >= x) && (xc_down < x+1)
                            flg_in_down = 1;
                            flg_in_right = 0;
                        elseif (yc_right >= y) && (yc_right < y+1)
                            flg_in_down = 0;
                            flg_in_right = 1;
                        else
                            flg_in_down = 0;
                            flg_in_right = 0;
                        end
                    else %傾きが正のとき左辺もしくは下辺から入る
                        if (xc_down >= x) && (xc_down < x+1)
                            flg_in_down = 1;
                            flg_in_left = 0;
                        elseif (yc_left >= y) && (yc_left < y+1)
                            flg_in_down = 0;
                            flg_in_left = 1;
                        else
                            flg_in_down = 0;
                            flg_in_left = 0;
                        end
                    end
                else %X線源が検出器より上にあるとき
                    if a(t,i) <= 0 %傾きが負のとき、上辺/左辺から入る
                        if (xc_up >= x) && (xc_up < x+1)
                            flg_in_up = 1;
                            flg_in_left = 0;
                        elseif (yc_left >= y) && (yc_left < y+1)
                            flg_in_up = 0;
                            flg_in_left = 1;
                        else
                            flg_in_up = 0;
                            flg_in_left = 0;
                        end
                    else %傾きが正のとき上辺/右辺から入る
                        if (xc_up >= x) && (xc_up < x+1)
                            flg_in_up = 1;
                            flg_in_right = 0;
                        elseif (yc_right >= y) && (yc_right < y+1)
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
                    if a(t,i) <= 0 % 傾きが負のとき上辺/左辺から出る
                        if (xc_up >= x) && (xc_up < x+1)
                            flg_out_up = 1;
                            flg_out_left = 0;
                        elseif (yc_left >= y) && (yc_left < y+1)
                            flg_out_up = 0;
                            flg_out_left = 1;
                        else
                            flg_out_up = 0;
                            flg_out_left = 0;
                        end
                    else %傾きが正のとき上辺/右辺から出る
                        if (xc_up >= x) && (xc_up < x+1)
                            flg_out_up = 1;
                            flg_out_right = 0;
                        elseif (yc_right >= y) && (yc_right < y+1)
                            flg_out_up = 0;
                            flg_out_right = 1;
                        else
                            flg_out_up = 0;
                            flg_out_right = 0;
                        end
                    end
                else %X線源が検出器より上にあるとき
                    if a(t,i) <= 0 %傾きが負のとき下辺/右辺から出る
                        if (xc_down >= x) && (xc_down < x+1)
                            flg_out_down = 1;
                            flg_out_right = 0;
                        elseif (yc_right >= y) && (yc_right < y+1)
                            flg_out_down = 0;
                            flg_out_right = 1;
                        else
                            flg_out_down = 0;
                            flg_out_right = 0;
                        end
                    else %傾きが正のとき下辺/左辺から出る
                        if (xc_down >= x) && (xc_down < x+1)
                            flg_out_down = 1;
                            flg_out_left = 0;
                        elseif (yc_left >= y) && (yc_left < y+1)
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
                    l = sqrt((xc_down-xc_up)^2 + 1);
                %下辺->左辺
                elseif (flg_in_down == 1) && (flg_out_left == 1)
                    l = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                %下辺->右辺
                elseif (flg_in_down == 1) && (flg_out_right == 1)
                    l = sqrt((yc_right-y)^2 + (xc_down-(x+1))^2);
                %右辺->上辺
                elseif (flg_in_right == 1) && (flg_out_up == 1)
                    l = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                %右辺->左辺
                elseif (flg_in_right == 1) && (flg_out_left == 1)
                    l = sqrt((yc_right-yc_left)^2 + 1);
                %右辺->下辺
                elseif (flg_in_right == 1) && (flg_out_down == 1)
                    l = sqrt((yc_right-y)^2 + (xc_down-(x+1))^2);
                %左辺->上辺
                elseif (flg_in_left == 1) && (flg_out_up == 1)
                    l = sqrt((xc_up-x)^2+(yc_left-(y+1))^2);
                %左辺->右辺
                elseif (flg_in_left == 1) && (flg_out_right == 1)
                    l = sqrt((yc_left-yc_right)^2 + 1);
                %左辺->下辺
                elseif (flg_in_left == 1) && (flg_out_down == 1)
                    l = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                %上辺->下辺
                elseif (flg_in_up == 1) && (flg_out_down == 1)
                    l = sqrt((xc_down-xc_up)^2 + 1);
                %上辺->左辺
                elseif (flg_in_up == 1) && (flg_out_left == 1)
                    l = sqrt((xc_up-x)^2+(yc_left-(y+1))^2);
                %上辺->右辺
                elseif (flg_in_up == 1) && (flg_out_right == 1)
                    l = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                else
                    l = 0;
                end
                
                j=sub2ind([YMAX,XMAX], y+1, x+1);
                td=sub2ind([DMAX,TMAX], i, t);
                if l~=0 %l(j-1)~=0 % 非ゼロ要素があるならスパース列に追加
                    TDs(k) = td;%(t-1)*DMAX+i;
                    Ps(k) = j;%j-1
                    ls(k) = l;%l(j-1)
                    k=k+1;
                end
            end
        end
    end
end
Lsparse = sparse(TDs,Ps,ls,TMAX*DMAX,XMAX*YMAX);

RGB = imread('C:\Users\katsube\Documents\MATLAB\M.png');
gray = rgb2gray(RGB);
GRAY_D = double(gray); %倍精度変換

imshow(gray)    %元画像の表示


%L = reshape(l,[XMAX*YMAX*DMAX,TMAX]);