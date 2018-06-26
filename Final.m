classdef Final
   properties (Constant)
       T=1;     %:= expiry date
       r=0.04;  %:= continiously compounded interest rate
       
       sig=0.2; %:= volatility of the asset
       S0=100;  %:= asset price at time 0 in €
       K=100;   %:= strike price in €
       
       N=10^4;  %:= number of simulations
       n=50;    %:= discretization steps
       mu=1.1   %:= parameter for excercise 4
      end
   methods (Static)
        %function to compute positive parts
        function PoP = positivePart(ST)
            for j = 1:length(ST)
                if ST(j)>Final.K
                    PoP(j) = ST(j)-Final.K;
                else
                    PoP(j) = 0;
                end
            end
        end
        %brownian motion approximation by gaussian random walk
        function B = gaussRW(n)
           Z = randn(1,n); %normally distributet random variables
           S = zeros(1,n); %summation
           B = zeros(1,n); %approximation of brownian motion
           for j = 1:n
               if j==1
                   S(j)=Z(j);
               else
                   S(j) = S(j-1)+Z(j);
               end
               B(j) = S(j)/sqrt(n); 
           end
        end
        %creates a huge array of brownian motins
        function SIMU = simulations(ammount)
           SIMU = zeros(ammount,50);
           for j = 1:ammount
               SIMU(j,:)=Final.gaussRW(50);
           end
       end
        %compute geometric brownian motion by euler sceme
        function STvec = geometric(amount, steps, B)
           y = ones(1,steps);
           SE = zeros(amount,steps);
           t = linspace(1/steps,1,steps);
           B_delta = zeros(amount,steps);
           for b = 1:amount
               B_delta(b,1) = B(b,1);
               SE(b,1) = 1+ Final.sig*B_delta(b,1);
               for k = 2:steps
                   B_delta(b,k) = B(b,k) - B(b,k-1);
                   SE(b,k) = SE(b,k-1) + Final.sig*SE(b,k-1)*B_delta(b,k);
               end
           end
        %We need just the ST so forget about the rest:
        STvec = SE(:,50);
        end
        function [Teta , left , right] = Exc1_i()
            %generate 10^4 N(0,T) distributet random variables WT
            WT = randn(1,Final.N)*sqrt(Final.T);
            ST = Final.S0*exp((Final.r-Final.sig^2/2)*Final.T + Final.sig*WT); 
            PoP = Final.positivePart(ST);
            MC = 0;
            for j = 1:Final.N
                MC = MC + PoP(j);
            end
            MC = MC / Final.N;
            Teta = exp(-Final.r * Final.T)*MC;
            %compute confidence Intervalls for alpha = 0.95
            X = exp(-Final.r * Final.T) .*PoP;
            [left,right] = Final.ConvInt(X,Teta);
        end
        function [left,right] = ConvInt(X,MC)
            X = X - MC;
            Var = sum(X.^2);
            Var = Var/(length(X)-1); 
            StD = sqrt(Var);
            %always take alpha = 0.95. Then we have C_alpha = 1.96
            left = MC - 1.96* (StD/sqrt(length(X)));
            right = MC + 1.96* (StD/sqrt(length(X)));
        end
        function [Teta , left , right] = Exc1_ii()
            %Number of Generations
            ammount = Final.N;
            %generate brownian motion with 50 discretisation steps
            SIMU = Final.simulations(ammount);
            %compute geometric brownian motion with euler sceme
            STvec = Final.geometric(ammount, Final.n, SIMU);
            %scale to our assed model
            STvec = Final.S0 * exp(Final.r * Final.T)*STvec;
            PoP = Final.positivePart(STvec);
            MC = 0;
            for j = 1:ammount
                MC = MC + PoP(j);
            end
            MC = MC / ammount;
            Teta = exp(-Final.r * Final.T)*MC;
            
            %compute confidence Intervalls for alpha = 0.95
            X = exp(-Final.r * Final.T) .*PoP;
            [left,right] = Final.ConvInt(X,Teta);
            
        end
       function C = Exc2
           d1 = 1/(Final.sig*sqrt(Final.T))*(log(Final.S0/Final.K) + (Final.r+Final.sig^2/2)*Final.T);
           d2 = d1 - Final.sig*sqrt(Final.T);
           %formula:
           C = Final.S0*normcdf(d1)-Final.K*exp(-Final.r*Final.T)*normcdf(d2);
       end
       function [Theta , left , right] = Exc3
          %Z = ST ; create 10^4 realisations
          ammount = Final.N;
          %generate brownian motion with 50 discretisation steps
          SIMU = Final.simulations(ammount);
          %compute geometric brownian motion with euler sceme
          Z = Final.geometric(ammount, Final.n, SIMU);
          %scale to our assed model Z = ST
          Z = Final.S0 * exp(Final.r * Final.T)*Z;
          
          %Y = exp(-rT) * PositivePart(ST-K) ; create 10^4 realisations
          Y = exp(-Final.r*Final.T)* Final.positivePart(Z);
          
          %%%%%%%%%%%STEP1%%%%%%%%%%%%%
          %approximation of c with smal integer ammount of realisations
          p = 400;
          EZ = 0;
          EY = 0;
          Var = 0;
          Cov = 0;
          for j = 1:p
              EZ = EZ + Z(j);
              EY = EY + Y(j);
          end
          EZ = EZ / p;
          EY = EY / p;
          for k = 1:p
              Var = Var + ((Z(k)-EZ)*(Z(k)-EZ));
              Cov = Cov + ((Y(k)-EY)*(Z(k)-EZ));
          end
          Var = Var / (p-1);
          Cov = Cov / (p-1);
          c = - (Cov/Var);
          
          %%%%%%%%%%STEP2%%%%%%%%%%%%%%%
          Theta = 0;
          for l = 1:Final.N
              Theta = Theta + Y(l)-c*(Z(l)-EZ);
          end
          Theta = Theta / Final.N;
          
          %compute confidence Intervalls
          X = zeros(1,Final.N);
          for l = 1:Final.N
              X(l) = Y(l) - c*(Z(l)-EZ);
          end
          [left,right] = Final.ConvInt(X,Theta);
       end
       function [Teta , left , right] = Exc4_i(my)
          %set
          WT = randn(1,Final.N);
          WT = WT * sqrt(Final.T);
          %compute the function f
          f = Final.S0*exp((Final.r- (Final.sig^2/2))*Final.T + Final.sig*(WT+my*Final.T));
          f = Final.positivePart(f);
          f = exp(-Final.r * Final.T).*f;
          %Monte Carlo
          X = zeros(1,Final.N);
          for j = 1: Final.N
              X(j) = exp(-my*WT(j)-(my^2/2)*Final.T)*f(j);
          end
          Teta = mean(X);
          %Confidence Intervalls
          [left,right] = Final.ConvInt(X,Teta);
       end
       function [best,record] = Exc4_ii()
           %let us do tests for mu going from 0.5 to 2 in steps of 0.1 with 20 sampels eace
           mu = 0.5;
           best = 0;
           record = 1000;
           while mu<2
               length = 0;
               for k = 1:20
                   [Theta,left,right] = Final.Exc4_i(mu);
                   length = length + right - left;
               end
               length = length / 20;
               if length < record
                   best = mu;
                   record = length;
               end
               mu = mu + 0.1;
           end
       end
       function [Teta , left , right] = Exc5_i(mst)
          [Teta , left , right] = Final.strata(mst);
       end
       function StD = Exc5_ii()
          Theta = zeros(1,7);
          left = zeros(1,7);
          right = zeros(1,7);
          StD = zeros(1,7);
          [Theta(1) , left(1) , right(1)]  = Final.strata(20);
          [Theta(2) , left(2) , right(2)]  = Final.strata(50);
          [Theta(3) , left(3) , right(3)]  = Final.strata(100);
          [Theta(4) , left(4) , right(4)]  = Final.strata(200);
          [Theta(5) , left(5) , right(5)]  = Final.strata(500);
          [Theta(6) , left(6) , right(6)]  = Final.strata(1000);
          [Theta(7) , left(7) , right(7)]  = Final.strata(2000);
          StD = right - Theta;
          StD = StD / 1.96;
          StD = StD * sqrt(Final.N);
          
       end
       function [Teta , left , right] = strata(m_st)
          %Using code of the presentation slides of Variance Reduction
          %Methods
           m_st = 500;
           p = 1/m_st;
           Nst = Final.N/m_st;
           M = zeros(1,m_st);
           Teta = 0;
           Var = 0;
           for j = 1 : m_st
               %make uniform random variables in the intervalls
               U = (j-1)/m_st + rand(1,Nst)/m_st;
               %generate normally distributet random variables for sampling
               %WT
               WT = norminv(U,0,Final.T);
               %generate E[(ST-K)+|WT in I]
               ST_I = Final.S0 * exp((Final.r - Final.sig^2/2)*Final.T + Final.sig.*WT);
               PoP = exp(-Final.r*Final.T)*Final.positivePart(ST_I);
               %Monte Carlo
               M(j) = mean(PoP);
               %Confidence
               Sum = sum(PoP);
               Sum_sq = sum(PoP.^2);
               Sig = (Sum_sq-(Sum^2)/Nst)/(Nst-1);
               Var = Var +Sig*p^2/Nst;
           end
           %Monte Carlo
           Teta = mean(M);
           left = Teta -1.96*sqrt(Var);
           right = Teta + 1.96*sqrt(Var);
       end
   end
end
