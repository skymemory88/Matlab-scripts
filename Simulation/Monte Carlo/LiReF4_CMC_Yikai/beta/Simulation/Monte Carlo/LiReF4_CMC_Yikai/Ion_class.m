classdef Ion_class
    properties
        name % magnetic ion species
        prop % ion proportion
        cfRot % crystal field rotation angle relative to crystal lattice
        hyp % hyerfine isotope concentration
        J % total magnetic moment
        L % orbital moment
        S % spin moment
        I % nuclear spin moment
        nLande % nuclear Lande factor
        gLande % electronic lande factor
        A % hyperfine coupling strength
        abc % lattice constant
        B % CEF constants
        h4 % h4 anisotropy constant
        renorm % total magnetic moment renormalization factor
        ex % exchange interaction strength
    end
    
    methods
        function obj = Ion(name, prop, cfRot, hyp, J, L, S, I, nLande, gLande, A, abc, B, h4, renorm, ex)
            % Constructor to initialize the Ion object
            obj.name = name;
            obj.prop = prop;
            obj.cfRot = cfRot;
            obj.hyp = hyp;
            obj.J = J;
            obj.L = L;
            obj.S = S;
            obj.I = I;
            obj.nLande = nLande;
            obj.gLande = gLande;
            obj.A = A;
            obj.abc = abc;
            obj.B = B;
            obj.h4 = h4;
            obj.renorm = renorm;
            obj.ex = ex;
            % Initialize magnetic matrices
            obj = obj.initializeMagneticMatrices();
        end

        function obj = setProp(obj, propValue)
            % Set property with validation
            if propValue >= 0 && propValue <= 1
                obj.prop = propValue;
            else
                error('Property value must be between 0 and 1.');
            end
        end
    end
end
