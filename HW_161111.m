% 日本語を書いてみるテスト 16/11/11 L計算(修正後)
% 32x32マス、(0,0)-(32,32)
% Xs:X線源のX座標
% Ys:X線源のY座標
% Xdi:i番目の検出器のX座標
% Ydi:i番目の検出器のY座標
% X線源とi番目の検出器を結ぶ直線は下記の連立方程式より解ける
%(y=ax+bとするとaが大きくなりすぎて不安定になるため、x=ay+bのかたちに置き換える)
% Xs = aYs + b
% Xdi = aYdi + b
% よって、a = (Xs-Xdi)/(Ys-Ydi), b = Xs-(Xs-Xdi)/(Ys-Ydi)*Ys
% (y=ax+bとする場合は上記のXとYを置き換えればよい)
%%%%%%
%%%%%%
% 事前指定
i = 1;
XMAX = 32;
YMAX = 32;
DMAX = 32;
Xs = XMAX/2;
Ys = -5;
i = 1:DMAX;
Xd(i) = i - 0.5;
Yd = YMAX-Ys;

%%%%%%
%%%%%%
%1x1の枠の中を線が通るか？
%各辺ごとにif文をつくる
%傾きで分岐

%交点の座標を求める
%(1)傾き
a(i) = (Xs-Xd(i))/(Ys-Yd);
b(i) = Xs-(Xs-Xd(i))/(Ys-Yd)*Ys;

%(2)交点
x = 0;
y = 0;

% 交点のx,y座標
Yc1(i) = 0:YMAX;
[Xc1(i),Yc1(i)] = [a(i) * Yc1 + b(i),Yc1];

Xc2 = 0:XMAX;
Yc2(i) = 1/a(i) * Xc2 - b(i)/a(i);


% j = 1;
% 
% flg_in_up = 0;
% flg_in_left = 0;
% flg_in_down = 0;
% flg_in_right = 0;
% 
% flg_out_up = 0;
% flg_out_left = 0;
% flg_out_down = 0;
% flg_out_right = 0;
% 
% for i = 1:DMAX
%     for x = 0:XMAX-1
%         for y = 0:YMAX-1
%             
%             a(i) = (Ys-Yd)/(Xs-Xd(i));
%             b(i) = Ys-(Ys-Yd)/(Xs-Xd(i))*Xs;
%             
%             %交点の座標
%             xc_up = (y+1-b(i))/a(i);
%             yc_left = a(i)*x+b(i);
%             xc_down = (y-b(i))/a(i);
%             yc_right = a(i)*(x+1)+b(i);
%             
%             %出口
%             if a(i) <= 0 % 傾きが負のとき、下辺/右辺から入って上辺/左辺から出る
%                 if (xc_up >= x) && (xc_up <= x+1)
%                     flg_out_up = 1;
%                     flg_out_left = 0;
%                 elseif (yc_left >= y) && (yc_left <= y+1)
%                     flg_out_up = 0;
%                     flg_out_left = 1;
%                 else
%                     flg_out_up = 0;
%                     flg_out_left = 0;
%                 end
%             else %傾きが正のとき左辺/下辺から入って上辺/右辺から出る
%                 if (xc_up >= x) && (xc_up <= x+1)
%                     flg_out_up = 1;
%                     flg_out_right = 0;
%                 elseif (yc_right >= y) && (yc_right <= y+1)
%                     flg_out_up = 0;
%                     flg_out_right = 1;
%                 else
%                     flg_out_up = 0;
%                     flg_out_right = 0;
%                 end
%             end
%             
%             %入口
%             if a(i) <= 0 %傾きが負のとき、下辺/右辺から入って上辺/左辺から出る
%                 if (xc_down >= x) && (xc_down <= x+1)
%                     flg_in_down = 1;
%                     flg_in_right = 0;
%                 elseif (yc_right >= y) && (yc_right <= y+1)
%                     flg_in_down = 0;
%                     flg_in_right = 1;
%                 else
%                     flg_in_down = 0;
%                     flg_in_right = 0;
%                 end
%             else %傾きが正のとき左辺/下辺から入って上辺/右辺から出る
%                 if (xc_down >= x) && (xc_down <= x+1)
%                     flg_in_down = 1;
%                     flg_in_left = 0;
%                 elseif (yc_left >= y) && (yc_left <= y+1)
%                     flg_in_down = 0;
%                     flg_in_left = 1;
%                 else
%                     flg_in_down = 0;
%                     flg_in_left = 0;
%                 end
%             end
%             
%             %長さの計算
%             %入り口と出口によって場合分け
%             
%             %下辺->上辺
%             if (flg_in_down == 1) && (flg_out_up == 1)
%                 l(j) = sqrt((xc_down-xc_up)^2 + 1);
%                 j = j + 1;
%             end
%             
%             %下辺->左辺
%             if (flg_in_down == 1) && (flg_out_left == 1)
%                 l(j) = sqrt((yc_left-y)^2 + (xc_down-x)^2);
%                 j = j + 1;
%             end
%             
%             %下辺->右辺
%             if (flg_in_down == 1) && (flg_out_right == 1)
%                 l(j) = sqrt((yc_right-y)^2 + (xc_down-x)^2);
%                 j = j + 1;
%             end
%             
%             %右辺->上辺
%             if (flg_in_right == 1) && (flg_out_up == 1)
%                 l(j) = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
%                 j = j + 1;
%             end
%             
%             %右辺->左辺
%             if (flg_in_right == 1) && (flg_out_left == 1)
%                 l(j) = sqrt((yc_right-yc_left)^2 + 1);
%                 j = j + 1;
%             end
%             
%             %左辺->上辺
%             if (flg_in_left == 1) && (flg_out_up == 1)
%                 l(j) = sqrt((xc_up-(x+1))^2+(yc_left-(y+1))^2);
%                 j = j + 1;
%             end
%             
%             %左辺->右辺
%             if (flg_in_left == 1) && (flg_out_right == 1)
%                 l(j) = sqrt((yc_left-yc_right)^2 + 1);
%                 j = j + 1;
%             end
%             
%             %交点なし
%             if (flg_in_left == 0) && (flg_in_down == 0) && (flg_in_right == 0)
%                 l(j) = 0;
%                 j = j + 1;
%             end
%         end
%     end
% end
% 
% L = reshape(l,[1024,32]);