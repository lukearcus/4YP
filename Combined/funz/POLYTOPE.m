classdef POLYTOPE <handle
    %POLYTOPE defined as 
    % P = {x: Ax <= b, Aeq x= beq, x>=lb, x<= ub }
    %(unless some cond are set to [],  eg if Aeq= beq = [] 
    %then def is P = {x: Ax <= b, Aeq x= beq, x>=lb, x<= ub }
    properties
    A =[];
    b =[];
    Aeq =[];
    beq =[];
    lb =[];
    ub =[];
    options=[];
    end
    
    methods
    
    
    function obj = POLYTOPE
            %obj.options = optimset('linprog');
            obj.options = optimset('lsqlin');
    %obj.options.MaxIter =9000000;
    obj.options.Display ='none';
     obj.options.Diagnostics ='off';
     obj.options.LargeScale ='off';
    end
       
          function answer = iselement(obj,x)
                answer =  (obj.project_onto(x) == x);
          end
    
      
        function xprojected = project_onto(obj,x)
         
             C = eye(length(x));
             
             if iscolumn(x) 
                  xprojected = lsqlin(C,x,obj.A,obj.b,obj.Aeq,obj.beq,obj.lb,obj.ub,[],obj.options);                 
                 
             else
                 xprojected = lsqlin(C,x',obj.A',obj.b',obj.Aeq',obj.beq',obj.lb,obj.ub,[],obj.options);
                 xprojected = xprojected';
             
             end
            
              
          end
   
    
        end
    
end


