function im_result= PoissonEditing(im1,im2,msk,bb)


[Lh Lv] = imgrad(im1);
[Gh Gv] = imgrad(im2);

X = im1;
Fh = Lh;
Fv = Lv;

msk=repmat(msk,[1 1 3]);

msk2=msk(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:);

X(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:) =  (~msk2).*X(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:) + im2.*msk2;
Fh(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:) = (~msk2).*Fh(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:) + Gh.*msk2;
Fv(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:) = (~msk2).*Fv(bb(2):bb(2)-1+bb(4),bb(1):bb(1)-1+bb(3),:) + Gv.*msk2;

%msk = zeros(size(X));
%msk(LY:LY+h,LX:LX+w,:) = 1;

%imwrite(uint8(X),'X.png');

%tic;
%Y = PoissonJacobi( X, Fh, Fv, msk );
%toc
%imwrite(uint8(Y),'Yjc.png');
tic;
im_result = PoissonGaussSeidel_( X, Fh, Fv, msk );
toc
%imwrite(uint8(Y),'Ygs.png');
