% ���{��������Ă݂�e�X�g 16/11/11 L�v�Z(�C����)
% 32x32�}�X�A(0,0)-(32,32)
% Xs:X������X���W
% Ys:X������Y���W
% Xdi:i�Ԗڂ̌��o���X���W
% Ydi:i�Ԗڂ̌��o���Y���W
% X������i�Ԗڂ̌��o������Ԓ����͉��L�̘A����������������
%(y=ax+b�Ƃ����a���傫���Ȃ肷���ĕs����ɂȂ邽�߁Ax=ay+b�̂������ɒu��������)
% Xs = aYs + b
% Xdi = aYdi + b
% ����āAa = (Xs-Xdi)/(Ys-Ydi), b = Xs-(Xs-Xdi)/(Ys-Ydi)*Ys
% (y=ax+b�Ƃ���ꍇ�͏�L��X��Y��u��������΂悢)
%%%%%%
%%%%%%
% ���O�w��
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
%1x1�̘g�̒�������ʂ邩�H
%�e�ӂ��Ƃ�if��������
%�X���ŕ���

%��_�̍��W�����߂�
%(1)�X��
a(i) = (Xs-Xd(i))/(Ys-Yd);
b(i) = Xs-(Xs-Xd(i))/(Ys-Yd)*Ys;

%(2)��_
x = 0;
y = 0;

% ��_��x,y���W
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
%             %��_�̍��W
%             xc_up = (y+1-b(i))/a(i);
%             yc_left = a(i)*x+b(i);
%             xc_down = (y-b(i))/a(i);
%             yc_right = a(i)*(x+1)+b(i);
%             
%             %�o��
%             if a(i) <= 0 % �X�������̂Ƃ��A����/�E�ӂ�������ď��/���ӂ���o��
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
%             else %�X�������̂Ƃ�����/���ӂ�������ď��/�E�ӂ���o��
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
%             %����
%             if a(i) <= 0 %�X�������̂Ƃ��A����/�E�ӂ�������ď��/���ӂ���o��
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
%             else %�X�������̂Ƃ�����/���ӂ�������ď��/�E�ӂ���o��
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
%             %�����̌v�Z
%             %������Əo���ɂ���ďꍇ����
%             
%             %����->���
%             if (flg_in_down == 1) && (flg_out_up == 1)
%                 l(j) = sqrt((xc_down-xc_up)^2 + 1);
%                 j = j + 1;
%             end
%             
%             %����->����
%             if (flg_in_down == 1) && (flg_out_left == 1)
%                 l(j) = sqrt((yc_left-y)^2 + (xc_down-x)^2);
%                 j = j + 1;
%             end
%             
%             %����->�E��
%             if (flg_in_down == 1) && (flg_out_right == 1)
%                 l(j) = sqrt((yc_right-y)^2 + (xc_down-x)^2);
%                 j = j + 1;
%             end
%             
%             %�E��->���
%             if (flg_in_right == 1) && (flg_out_up == 1)
%                 l(j) = sqrt((xc_up-(x+1))^2+(yc_right-(y+1))^2);
%                 j = j + 1;
%             end
%             
%             %�E��->����
%             if (flg_in_right == 1) && (flg_out_left == 1)
%                 l(j) = sqrt((yc_right-yc_left)^2 + 1);
%                 j = j + 1;
%             end
%             
%             %����->���
%             if (flg_in_left == 1) && (flg_out_up == 1)
%                 l(j) = sqrt((xc_up-(x+1))^2+(yc_left-(y+1))^2);
%                 j = j + 1;
%             end
%             
%             %����->�E��
%             if (flg_in_left == 1) && (flg_out_right == 1)
%                 l(j) = sqrt((yc_left-yc_right)^2 + 1);
%                 j = j + 1;
%             end
%             
%             %��_�Ȃ�
%             if (flg_in_left == 0) && (flg_in_down == 0) && (flg_in_right == 0)
%                 l(j) = 0;
%                 j = j + 1;
%             end
%         end
%     end
% end
% 
% L = reshape(l,[1024,32]);