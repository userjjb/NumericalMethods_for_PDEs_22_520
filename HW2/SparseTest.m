tic
Neta=5;
Nxi=5;
NNodes=Nxi*Neta;

NodeNumber=zeros(Nxi,Neta);
counter = 0;
for j = 1:Neta
	for i = 1:Nxi
		counter = counter + 1;
		NodeNumber(i,j) = counter; 
    end
end

%A = spalloc(NNodes,NNodes, NNodes*10);
%A=zeros(NNodes,NNodes);
F=ones(NNodes);

% for j = 2:(Neta - 1)
%     for i = 2:(Nxi-1)
%             NodeN = NodeNumber(i,j);
%             A(NodeN, NodeNumber(i+0,j+0)) =  19947; 
% 
%             A(NodeN, NodeNumber(i+1,j+0)) =  -4987;
%             A(NodeN, NodeNumber(i-1,j+0)) =  -4987;
% 
%             A(NodeN, NodeNumber(i+0,j+1)) =  -4987;
%             A(NodeN, NodeNumber(i+0,j-1)) =  -4987;
% 
%             A(NodeN, NodeNumber(i+1,j+1)) =  1054;
%             A(NodeN, NodeNumber(i+1,j-1)) =  -1054;
%             A(NodeN, NodeNumber(i-1,j+1)) =  -1054;
%             A(NodeN, NodeNumber(i-1,j-1)) =  1054;
%     end
% end

DiagValue = ones(NNodes,9);
% DiagValue([1:0, NNodes-2*Nxi-1:NNodes],1)=0;
% DiagValue([1:1, NNodes-2*Nxi+0:NNodes],2)=0;
% DiagValue([1:2, NNodes-2*Nxi+1:NNodes],3)=0;
% 
% DiagValue([1:0+Nxi, NNodes-Nxi-1:NNodes],4)=0;
% DiagValue([1:1+Nxi, NNodes-Nxi+0:NNodes],5)=0;
% DiagValue([1:2+Nxi, NNodes-Nxi+1:NNodes],6)=0;
% 
% DiagValue([1:0+2*Nxi, NNodes-1:NNodes],7)=0;
% DiagValue([1:1+2*Nxi, NNodes-0:NNodes],8)=0;
% DiagValue([1:2+2*Nxi, NNodes+1:NNodes],9)=0;


for cousin=0:2
    for brother=0:2
        for iter=1:Neta-3
            R=Nxi*iter-1:Nxi*iter;
            DiagValue(R+brother+cousin*Nxi,cousin*3+brother+1)=0;
        end
        DiagValue([1:brother+cousin*Nxi, NNodes-(2-cousin)*Nxi+brother-1:NNodes],cousin*3+brother+1)=0;
    end
end
    

DiagValue = DiagValue.*99;
Q=spdiags(DiagValue, [-Nxi-1 -Nxi -Nxi+1 -1 0 1 Nxi-1 Nxi Nxi+1], Nxi*Neta, Nxi*Neta);
Q=full(Q)

Solution = A\F;
toc
beep
