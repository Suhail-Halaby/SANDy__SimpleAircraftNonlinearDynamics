classdef VehicleModel < matlab.System
    % VehicleModel  Ground-truth vehicle dynamics wrapper around runVehicle

    properties (Nontunable)
        Plane        % the full 'plane' struct, StripModel objects and all
        ModelName    % string, e.g. "simplePlate"
        Rho = 1.22   % air density
    end

    methods (Access = protected)

        function setupImpl(obj)
            % Called once at sim start — nothing required here since
            % Plane is passed in as a property, not built here.
        end

        function [F, M, Mass, I, CG] = stepImpl(obj, CTRL, BodyRates, BodyVel, OrQuat, Height)
            [F, M, Mass, CG, I] = runVehicle(obj.Plane, CTRL, BodyRates, ...
                BodyVel, OrQuat, Height, ...
                obj.Rho, obj.ModelName);
        end

        function num = getNumInputsImpl(~)
            num = 5;   % CTRL, BodyRates, BodyVel, OrQuat, Height
        end

        function num = getNumOutputsImpl(~)
            num = 5;   % F, M, Mass, I, CG
        end

        function [n1,n2,n3,n4,n5] = getOutputNamesImpl(~)
            n1 = 'F'; n2 = 'M'; n3 = 'Mass'; n4 = 'I'; n5 = 'CG';
        end

        function flag = isInputSizeMutableImpl(~,~)
            flag = false;
        end

        function [s1,s2,s3,s4,s5] = getOutputSizeImpl(obj)
            s1 = [3 1];   % F    (3x1 force vector)
            s2 = [3 1];   % M    (3x1 moment vector)
            s3 = [1 1];   % Mass (scalar)
            s4 = [3 3];   % I    (3x3 inertia tensor)
            s5 = [3 1];   % CG   (3x1 position vector)
        end

        function [t1,t2,t3,t4,t5] = getOutputDataTypeImpl(obj)
            t1 = 'double';
            t2 = 'double';
            t3 = 'double';
            t4 = 'double';
            t5 = 'double';
        end

        function [c1,c2,c3,c4,c5] = isOutputComplexImpl(obj)
            c1 = false;
            c2 = false;
            c3 = false;
            c4 = false;
            c5 = false;
        end

        function [f1,f2,f3,f4,f5] = isOutputFixedSizeImpl(obj)
            f1 = true;
            f2 = true;
            f3 = true;
            f4 = true;
            f5 = true;
        end

    end
end