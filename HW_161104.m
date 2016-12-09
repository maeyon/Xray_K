% 16/11/3 L�v�Z
% 32x32�}�X�A(0,0)-(32,32)
% Xs(t):t���ڂ�X������X���W
% Ys(t):t���ڂ�X������Y���W
% Xd(t,i):t����i�Ԗڂ̌��o���X���W
% Yd(t,i):t����i�Ԗڂ̌��o���Y���W
% X������i�Ԗڂ̌��o������Ԓ����͉��L�̘A����������������
% Ys = aXs + b
% Ydi = aXdi + b
% ����āAa = (Ys-Ydi)/(Xs-Xdi), b = Ys-(Ys-Ydi)/(Xs0-Xdi)*Xs
%
%%%%%%
%%%%%%
% ���O�w��
i = 1;
XMAX = 32;
YMAX = 32;
DMAX = 32;

Xs0 = XMAX/2;   %��]�p�[����X�������W
Ys0 = -5;       %��]�p�[����X�������W
TMAX = 30;     %���e����
XCENTER = XMAX/2;      %��]���S
YCENTER = YMAX/2;      %��]���S

%��]�p�[���̌��o����W
Xd0 = XMAX/2;
Yd0 = YMAX-Ys0;

%��]�p���l������t���ڂ̓��e��X�����̍��W(���W�A������)
t = 1:TMAX;
Xs(t) = (Xs0-XCENTER)*cos(t*2*pi/TMAX)-(Ys0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER;
Xs = transpose(repmat(Xs,DMAX,1));
Ys(t) = (Xs0-XCENTER)*sin(t*2*pi/TMAX)+(Ys0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER;
Ys = transpose(repmat(Ys,DMAX,1));


%��]�p���l������t����,i�Ԗڂ̓��e�̌��o��̍��W(���W�A������)
t = 1:TMAX;
    for i = 1:DMAX
        Xd(t,i) = (Xd0-XCENTER)*cos(t*2*pi/TMAX)-(Yd0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER-(DMAX/2-i+0.5)*cos(t*2*pi/TMAX);
        Yd(t,i) = (Xd0-XCENTER)*sin(t*2*pi/TMAX)+(Yd0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER+(DMAX/2-i+0.5)*sin(t*2*pi/TMAX);
    end

%%%%%%
%%%%%%
%1x1�̘g�̒�������ʂ邩�H
%�e�ӂ��Ƃ�if��������
%�X���ŕ���

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
                
                %��_�̍��W
                xc_up = (y+1-b(t,i))/a(t,i);
                yc_left = a(t,i)*x+b(t,i);
                xc_down = (y-b(t,i))/a(t,i);
                yc_right = a(t,i)*(x+1)+b(t,i);

                %����
                if Ys(t,i) <= YCENTER %X���������o���艺�ɂ���Ƃ�
                    if a(t,i) <= 0 %�X�������̂Ƃ��A����/�E�ӂ�������ď��/���ӂ���o��
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
                    else %�X�������̂Ƃ�����/���ӂ�������ď��/�E�ӂ���o��
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
                else %X���������o�����ɂ���Ƃ�
                    if a(t,i) <= 0 %�X�������̂Ƃ��A���/���ӂ�������ĉ���/�E�ӂ���o��
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
                    else %�X�������̂Ƃ����/�E�ӂ�������ĉ���/���ӂ���o��
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
                
                %�o��
                if Ys(t,i) <= YCENTER %X���������o���艺�ɂ���Ƃ�
                    if a(t,i) <= 0 % �X�������̂Ƃ��A����/�E�ӂ�������ď��/���ӂ���o��
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
                    else %�X�������̂Ƃ�����/���ӂ�������ď��/�E�ӂ���o��
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
                else %X���������o�����ɂ���Ƃ�
                    if a(t,i) <= 0 %�X�������̂Ƃ��A���/���ӂ�������ĉ���/�E�ӂ���o��
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
                    else %�X�������̂Ƃ����/�E�ӂ�������ĉ���/���ӂ���o��
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

                %�����̌v�Z
                %������Əo���ɂ���ďꍇ����
                
                %����->���
                if (flg_in_down == 1) && (flg_out_up == 1)
                    l(j) = sqrt((xc_down-xc_up)^2 + 1);
                    j = j + 1;
                %����->����
                elseif (flg_in_down == 1) && (flg_out_left == 1)
                    l(j) = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %����->�E��
                elseif (flg_in_down == 1) && (flg_out_right == 1)
                    l(j) = sqrt((yc_right-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %�E��->���
                elseif (flg_in_right == 1) && (flg_out_up == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                    j = j + 1;
                %�E��->����
                elseif (flg_in_right == 1) && (flg_out_left == 1)
                    l(j) = sqrt((yc_right-yc_left)^2 + 1);
                    j = j + 1;
                %�E��->����
                elseif (flg_in_right == 1) && (flg_out_down == 1)
                    l(j) = sqrt((yc_right-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %����->���
                elseif (flg_in_left == 1) && (flg_out_up == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_left-(y+1))^2);
                    j = j + 1;
                %����->�E��
                elseif (flg_in_left == 1) && (flg_out_right == 1)
                    l(j) = sqrt((yc_left-yc_right)^2 + 1);
                    j = j + 1;
                %����->����
                elseif (flg_in_left == 1) && (flg_out_down == 1)
                    l(j) = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                    j = j + 1;
                %���->����
                elseif (flg_in_up == 1) && (flg_out_down == 1)
                    l(j) = sqrt((xc_down-xc_up)^2 + 1);
                    j = j + 1;
                %���->����
                elseif (flg_in_up == 1) && (flg_out_left == 1)
                    l(j) = sqrt((xc_up-(x+1))^2+(yc_left-(y+1))^2);
                    j = j + 1;
                %���->�E��
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