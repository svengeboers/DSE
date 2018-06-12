function a  = elNodeReaction(qy,qz,L)
% 	Description: generates array of contribution of a discretized 
% 				 distributed force acting in y for an 8 DoF beam
% 				 element
% 	
% 	Inputs:
% 		q       =  Distributed load magnitude [N/m]
% 		L       =  Young's Modulus [Pa]
% 	
% 	Outputs:
% 		array   =  q force array
		
	
	a = [qy*L/2,qy*L*L/12,qz*L/2,-qz*L*L/12,...
        qy*L/2,-qy*L*L/12,qz*L/2,qz*L*L/12]';
    
end
