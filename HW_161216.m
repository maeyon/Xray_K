% 16/12/14 L�v�Z
% [16/12/9�̏h������ƂɕύX�_]
% (��)(1)l(j)->l�ɕύX
% (2)�X��������̑΍�
% (3)���e���̕\��
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
XMAX = 3;  %X��f��
YMAX = 3;  %Y��f��
DMAX = 2; %���o��̐�
DWIDTH = 1;

Xs0 = XMAX/2;   %��]�p�[����X������X���W
Ys0 = -5;       %��]�p�[����X������Y���W
TMAX = 3;     %���e����
XCENTER = XMAX/2;      %��]���S
YCENTER = YMAX/2;      %��]���S

%Xd0 = XMAX/2;   %��]�p�[���̌��o���X���W
Xd0 = [0:XMAX/(DMAX-1):XMAX];   %��]�p�[���̌��o���X���W
Yd0 = (YMAX-Ys0)*ones(1,DMAX);  %��]�p�[���̌��o���Y���W
Genten_Hosei_D = [Xd0 - XCENTER; Yd0 - YCENTER];    %��]�p�[���̌��o��̍��W(���_�Ɉړ�)

%��]�p���l������t���ڂ̓��e��X�����̍��W(���W�A������)
clear Xs Ys Xd Yd
for t = 1:TMAX
    % ��]�s��R�@�́@�쐬
    Rt = [cos((t-1)*2*pi/TMAX) -sin((t-1)*2*pi/TMAX);sin((t-1)*2*pi/TMAX) cos((t-1)*2*pi/TMAX)];
    %�@���_�␳��̌��o����W�̉�]
    temp = Rt*Genten_Hosei_D;
    %�@���Ƃ̍��W�ʒu�ɕ���
    Xd(t,:) = temp(1,:) + XCENTER;
    Yd(t,:) = temp(2,:) + YCENTER;
    
    %�@���_�␳���X�������W�̉�]
    temp = Rt*[Xs0-XCENTER;Ys0-YCENTER];
    % ���Ƃ̍��W�ʒu�ɕ���
    Xs(t) = temp(1) + XCENTER; 
    Ys(t) = temp(2) + YCENTER;
end

Xs = repmat(Xs(:),1,DMAX);
Ys = repmat(Ys(:),1,DMAX);

% tic;
% % ��]�s��R�@�́@�쐬
% Rall = [cos([1:TMAX]'*2*pi/TMAX) -sin([1:TMAX]'*2*pi/TMAX);sin([1:TMAX]'*2*pi/TMAX) cos([1:TMAX]'*2*pi/TMAX)];
% %�@���_�␳��̍��W�̉�]
% temp = Rall*Genten_Hosei_D;
% %�@���_�␳��̍��W����]�������̂����Ƃ̍��W�ʒu�ɕ���
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
% %��]�p���l������t����,i�Ԗڂ̓��e�̌��o��̍��W(���W�A������)
% %�����C��
% t = 1:TMAX;
%     for i = 1:DMAX
%         Xd(t,i) = (Xd0-XCENTER)*cos(t*2*pi/TMAX)-(Yd0-YCENTER)*sin(t*2*pi/TMAX)+XCENTER-(DMAX/2-i+0.5)*cos(t*2*pi/TMAX);
%         Yd(t,i) = (Xd0-XCENTER)*sin(t*2*pi/TMAX)+(Yd0-YCENTER)*cos(t*2*pi/TMAX)+YCENTER+(DMAX/2-i+0.5)*sin(t*2*pi/TMAX);
%     end


    
%%%%%%
%%%%%%
%1x1�̘g�̒�������ʂ邩�H
%�e�ӂ��Ƃ�if��������
%�X���ŕ���

x = 0;
y = 0;
j = 1;

a = (Ys-Yd)./(Xs-Xd);   %�X��
b = Ys-(Ys-Yd)./(Xs-Xd).*Xs;    %�ؕ�

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

                %��_�̍��W
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

                %����
                if Ys(t,i) <= YCENTER %X���������o���艺�ɂ���Ƃ�
                    if a(t,i) <= 0 %�X�������̂Ƃ��A���ӂ������͉E�ӂ������
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
                    else %�X�������̂Ƃ����ӂ������͉��ӂ������
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
                else %X���������o�����ɂ���Ƃ�
                    if a(t,i) <= 0 %�X�������̂Ƃ��A���/���ӂ������
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
                    else %�X�������̂Ƃ����/�E�ӂ������
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
                
                %�o��
                if Ys(t,i) <= YCENTER %X���������o���艺�ɂ���Ƃ�
                    if a(t,i) <= 0 % �X�������̂Ƃ����/���ӂ���o��
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
                    else %�X�������̂Ƃ����/�E�ӂ���o��
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
                else %X���������o�����ɂ���Ƃ�
                    if a(t,i) <= 0 %�X�������̂Ƃ�����/�E�ӂ���o��
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
                    else %�X�������̂Ƃ�����/���ӂ���o��
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

                %�����̌v�Z
                %������Əo���ɂ���ďꍇ����
                
                %����->���
                if (flg_in_down == 1) && (flg_out_up == 1)
                    l = sqrt((xc_down-xc_up)^2 + 1);
                %����->����
                elseif (flg_in_down == 1) && (flg_out_left == 1)
                    l = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                %����->�E��
                elseif (flg_in_down == 1) && (flg_out_right == 1)
                    l = sqrt((yc_right-y)^2 + (xc_down-(x+1))^2);
                %�E��->���
                elseif (flg_in_right == 1) && (flg_out_up == 1)
                    l = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                %�E��->����
                elseif (flg_in_right == 1) && (flg_out_left == 1)
                    l = sqrt((yc_right-yc_left)^2 + 1);
                %�E��->����
                elseif (flg_in_right == 1) && (flg_out_down == 1)
                    l = sqrt((yc_right-y)^2 + (xc_down-(x+1))^2);
                %����->���
                elseif (flg_in_left == 1) && (flg_out_up == 1)
                    l = sqrt((xc_up-x)^2+(yc_left-(y+1))^2);
                %����->�E��
                elseif (flg_in_left == 1) && (flg_out_right == 1)
                    l = sqrt((yc_left-yc_right)^2 + 1);
                %����->����
                elseif (flg_in_left == 1) && (flg_out_down == 1)
                    l = sqrt((yc_left-y)^2 + (xc_down-x)^2);
                %���->����
                elseif (flg_in_up == 1) && (flg_out_down == 1)
                    l = sqrt((xc_down-xc_up)^2 + 1);
                %���->����
                elseif (flg_in_up == 1) && (flg_out_left == 1)
                    l = sqrt((xc_up-x)^2+(yc_left-(y+1))^2);
                %���->�E��
                elseif (flg_in_up == 1) && (flg_out_right == 1)
                    l = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
                else
                    l = 0;
                end
                
                j=sub2ind([YMAX,XMAX], y+1, x+1);
                td=sub2ind([DMAX,TMAX], i, t);
                if l~=0 %l(j-1)~=0 % ��[���v�f������Ȃ�X�p�[�X��ɒǉ�
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
GRAY_D = double(gray); %�{���x�ϊ�

imshow(gray)    %���摜�̕\��


%L = reshape(l,[XMAX*YMAX*DMAX,TMAX]);