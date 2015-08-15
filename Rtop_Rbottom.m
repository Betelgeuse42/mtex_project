function [Rtop,Rbottom] =Rtop_Rbottom(h,ind,dx,dy)
        %Top coordinate matrix
        Rtop=zeros(3,h+1);
        %Bottom coordinate matrix
        Rbottom=zeros(3,h+1);
for c=1:h
                    g=ind(c,1); 
                    %coordinate matrix
                    Rtop(1,(c+0))=dx(g);
                    Rtop(2,(c+0))=dy(g);
                    Rtop(3,(c+0))=1;
                    Rbottom(1,(c+0))=dx(g);
                    Rbottom(2,(c+0))=dy(g);
                    Rbottom(3,(c+0))=-1;
end;
                    Rtop(1,c+1)=0;
                    Rtop(2,c+1)=0;
                    Rtop(3,c+1)=1;
                    Rbottom(1,c+1)=0;
                    Rbottom(2,c+1)=0;
                    Rbottom(3,c+1)=-1;
                   