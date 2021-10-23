clear;
syms t1 t2 t3
syms l0 l1

psi1 = [0;0;0;0;0;1];
psi2 = [0;-l0;  0 ; -1 ; 0 ;0];
psi3 = [-l0;0;0;0;1;0]; 
s1 = [ cos(t1) -sin(t1) 0 0 ; sin(t1) cos(t1) 0 0 ; 0 0 1 0 ; 0 0 0 1] ;
s2 = [1 0 0 0 ; 0 cos(t2) sin(t2) -l0*sin(t2); 0 -sin(t2) cos(t2) l0*(1-cos(t2)) ; 0 0 0 1 ] ;
s3 = [ cos(t3) 0 sin(t3) -l0*sin(t3) ; 0 1 0 0 ; -sin(t3) 0 cos(t3) l0*(1-cos(t3)); 0 0 0 1] ;
gst = [ 1 0 0 0 ; 0 1 0 l1 ; 0 0 1 l0 ; 0 0 0 1] ; 


% product of transformation matrices 
final = s1*s2*s3*gst ;
inter1 = s1*s2  ; 
inter2 = inter1*s3 ; 
final = simplify(final);
inter1 = simplify(inter1) ;
inter2 = simplify(inter2);

% Forward Kinematics map
final



% Getting the adjoint transformation matrices
adj_s1 = createAdj(s1);
adj_s1_s2 = createAdj(inter1);


% Modified twists after adjoint transformation for calculating spatial Jacobian 
npsi1 = psi1 ;
npsi2 = simplify(adj_s1*psi2) ;
npsi3 = simplify(adj_s1_s2*psi3);

spatial_jacobian = [npsi1 npsi2 npsi3]

% Inverse transformations for the purpose of calculating inverse adjoint
% transformation

inv_1 = s3*gst;
inv_2 = s2*inv_1;
inv_3 = s1*inv_2;


% Inverse Adjoint Transformation.
inv_adj_s3 = createInvAdj(inv_1);
inv_adj_s2_s3 = createInvAdj(inv_2);
inv_adj_s1_s2_s3 = createInvAdj(inv_3);


%Modified twists after inverse adjoint transformation for calculating the
%Body Jacobian.
i_psi_1 = simplify(inv_adj_s1_s2_s3*psi1);
i_psi_2 = simplify(inv_adj_s2_s3*psi2);
i_psi_3 = simplify(inv_adj_s3*psi3);

body_jacobian = [i_psi_1 i_psi_2 i_psi_2]


% Function for calculating the Inverse adjoint transformation matrix corresponding
% to a given homogeneous transformation given as input.
function i_adj_ = createInvAdj(T)
   rotation_matrix = sym([]);
    i_adj_ = sym([]);
    for i  = 1:3
        for j = 1:3
            i_adj_(j,i) = simplify(T(i,j));
            i_adj_(j+3,i+3) = simplify(T(i,j));
            rotation_matrix(j,i) = simplify(T(i,j));
        end
    end
    p_hat = [ 0 -simplify(T(3,4)) simplify(T(2,4)) ; simplify(T(3,4)) 0 -simplify(T(1,4)) ; -simplify(T(2,4)) simplify(T(1,4)) 0];
    adj_p =  -rotation_matrix*p_hat;
    for i = 1 :3
        for j = 4:6
            i_adj_(i,j) = simplify(adj_p(i,j-3));
        end
    end

end




% Function for calculating the adjoint transformation matrix corresponding
% to a given homogeneous transformation given as input.
function adj_  = createAdj(T)
    rotation_matrix = sym([]);
    adj_ = sym([]);
    for i  = 1:3
        for j = 1:3
            adj_(i,j) = simplify(T(i,j));
            adj_(i+3,j+3) = simplify(T(i,j));
            rotation_matrix(i,j) = simplify(T(i,j));
        end
    end
    p_hat = [ 0 -simplify(T(3,4)) simplify(T(2,4)) ; simplify(T(3,4)) 0 -simplify(T(1,4)) ; -simplify(T(2,4)) simplify(T(1,4)) 0];
    adj_p = p_hat*rotation_matrix;
    for i = 1 :3
        for j = 4:6
            adj_(i,j) = simplify(adj_p(i,j-3));
        end
    end
end
