%computational assignment che221
%Equation of state - Soave-Redlich-Kwong (SRK)
%system -Cumene
%varun khatri
%190946

classdef srk
    properties
        Tc_exp;
        Pc_exp;
        Vc_exp;
        omega;
        Tc;
        Pc;
        Vc;
        a;
        b;
        alpha;
        m;
        Tr;
        R = 0.0083;  
    end
    methods
        %constructor
        function [self] = srk(T,P,V,omega)
            self.Tc_exp = T;
            self.Pc_exp = P;
            self.Vc_exp = V;
            self.omega = omega;
        end
        function [self] = calculate_constant_for_EOS(self)
            %for SRK
            %https://www.e-education.psu.edu/png520/m10_p5.html
            %a = 0.427480*(R^2*Tc^2)/Pc
            %b = 0.086640*(R*Tc/Pc)
            
            self.a = 0.427480 * (self.R * self.Tc_exp) ^ 2 / self.Pc_exp;
            self.b = 0.086640 * (self.R * self.Tc_exp / self.Pc_exp);
            self.m = 0.48508 + 1.55171 * self.omega - 0.1561 * self.omega ^ 2;
            
        end
        function [self] = calculate_alpha(self,T)
            self.Tr = T/self.Tc_exp;
            self.alpha = (1 + self.m * (1 - sqrt(self.Tr))) ^ 2;
        end
        
        function [z] = calculate_z(self,T,v)
            rho = 1/v;
            z = 1 / (1 - rho* self.b) - (self.alpha * self.a * rho / (self.R * T *(1 + rho * self.b)));
        end
        function [d] = departure_function(self,T,v)
            rho = 1/v;
            z = self.calculate_z(T,v);
            d = -1 * log(1 - rho * self.b) -1 * (self.alpha * self.a / (self.R * self.b * T)) * log(1 + rho * self.b) + (z - 1) - log(z);
        end
        
        function [v] = find_volume(self,P,T)
            
            A = P;
            B = -1 * self.R * T;
            C = -1 * (self.R * T * self.b + P * (self.b)^2 - self.alpha * self.a);
            D = -1 * self.alpha * self.a * self.b;
            
            coeff  = [A B C D];
            
            v = roots(coeff);
            %only real roots required
            v = v(imag(v) ==0);
        end
        
        function[p] = find_pressure(self,T,v) 
            p = self.R * T /(v-self.b) - (self.alpha * self.a) / (v * (v+ self.b));
        end

        function [value] = chem_potential_diff(self,T,v_l,v_g) 
            value = self.departure_function(T,v_l) - self.departure_function(T,v_g);
        end

        function [value] = find_value(self,T,P)
            v = self.find_volume(P,T);
            if( size(v) == 1) 
                value = self.departure_function(T,v);
            else 
                v = sort(v);
                value = self.chem_potential_diff(T,v(1),v(3));
            end
        end
        function plot(self,T,l_limit,u_limit)
            y = @(x) ((self.R*T)/(x-self.b))-(self.alpha* self.a /((x)*(x + self.b)));
            fplot(@(x) y(x),[l_limit, u_limit],'r','LineWidth',1)
            ylim([0,6]);
        end
        function [limit] = find_lower_limit(self,T)
            P = 0.000001;
            while 1
                v = self.find_volume(P,T);
                s = size(v);
                if( s(1)== 3)
                    break
                end
                if( P > self.Pc_exp)
                    break
                end
                P = P + 0.1;
            end
            limit = P;
        end
        function [limit] = find_upper_limit(self,T,iterator,limit)
            P = self.Pc_exp;
            while 1 
                v = self.find_volume(P,T);
                s = size(v);
                if( s(1)== 3)
                    break
                end
                if( P < limit)
                    break
                end
                P = P - iterator;
            end
            limit = P;
        end
                
        function [Pc, Tc, Vc , upper_dome_region] = find_critical_temp_and_pressure(self)
            % let guess be 400 kelvin
            Tc = 400;
            self = self.calculate_alpha(Tc);
            Pc = self.find_psat(Tc);
            v = sort(self.find_volume(Pc , Tc));
            error =  abs(v(1) - v(3));
            while 1
                self = self.calculate_alpha(Tc);
                Pc = self.find_psat(Tc);
                v = sort(self.find_volume(Pc , Tc));
                error =  abs(v(1) - v(3));
                if error > 1
                    Tc = Tc + 10;
                elseif (error > 0.65)
                    Tc = Tc +5;
                else
                %in this region srk and other eos do not behave properly 
                %as it is very near to criticlal point
                    break;
                end
                
            end
            %now , Tc is very close to critical temperature , so we can now adopt 
            %a very intesive computational method to calculate Tc and Pc
            upper_dome_region = [];
            while error > 0.35
                Tc = Tc + 1;
                self = self.calculate_alpha(Tc);
                Pc = 0.5*(self.find_upper_limit(Tc,0.1,Pc) + self.find_lower_limit(Tc));
                v = sort(self.find_volume(Pc,Tc));
                error = abs(v(1) - v(3));
                upper_dome_region = [upper_dome_region ; [Tc Pc v(1) v(3)] ];
            end
            while error > 0.1
                Tc = Tc + 0.1;
                self = self.calculate_alpha(Tc);
                Pc = self.find_upper_limit(Tc,0.0001,Pc);
                v = sort(self.find_volume(Pc,Tc));
                error = abs(v(1) - v(3));
                upper_dome_region = [upper_dome_region ; [Tc Pc v(1) v(3)] ];
            end 
             while error > 0.001
                Tc = Tc + 0.001;
                self = self.calculate_alpha(Tc);
                Pc = self.find_upper_limit(Tc,0.0001,Pc);
                v = sort(self.find_volume(Pc,Tc));
                if size(v) == 1
                    break;
                end
                error = abs(v(1) - v(3));
                upper_dome_region = [upper_dome_region ; [Tc Pc v(1) v(3)] ];
             end 
            upper_dome_region = [upper_dome_region ; [Tc Pc v(1) v(1)] ];
            Vc = v(1);  
            self.Pc = Pc;
            self.Tc = Tc;
            self.Vc = Vc;
            

        end
                    
        function [mid_value] = find_psat(self,T)
            % upper limit cannot be more than critical temperature of the system
            upper_limit = self.Pc_exp; 
            %lower limit can be equal to minimum pressure at which we get three real roots i.e start of two phase region
            lower_limit = self.find_lower_limit(T); 
            
            max_error = 0.0001; 
            error = 1000000; % any arbitrary value of error for the start
            mid_value = upper_limit;
            
            % bisection method 
            while(error > max_error) 
                prev_mid_value = mid_value;
                mid_value = ( upper_limit + lower_limit ) / 2;
                value_upper_limit = self.find_value(T,upper_limit);
                value_lower_limit = self.find_value(T,lower_limit);
                value_mid =         self.find_value(T,mid_value);

                if value_upper_limit * value_lower_limit > 0 
                     disp("error : wrong inital boundary condition are choosen");
                    return;
                end

                if abs(value_mid) < 0.0000000000000000001 %nearly equal to zero 
                    break;
                elseif value_upper_limit * value_mid < 0
                    lower_limit = mid_value;
                else
                    upper_limit = mid_value;
                end
                
                error = (abs(mid_value - prev_mid_value)/abs(mid_value)) * 100;
            end

        end
                
    end
end

    
 