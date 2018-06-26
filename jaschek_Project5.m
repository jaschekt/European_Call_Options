function main

    FirstExc;
    SecondExc;
    ThirdExc;
    FourthExc;
    FifthExc;
    
    disp('********************THE**END*****************')
    friendlyPrint('I hope you had a great experience and will never ever have problems with the estimation of premiums for european call options.');
    close all
    
    function FirstExc
        disp('***Welcome to the FIRST exercise***');
        friendlyPrint('press a button to see result of exercise 1.1');
        disp('Plain Monte Carlo with geometric Brownian motion:');
        [Theta,left,right] = Final.Exc1_i;
        Theta
        ConfI = [left,right]
        length = right - left
        
        friendlyPrint('press a button to see result of exercise 1.2');
        disp('Plain Monte Carlo with Euler scheme:');
        [Theta,left,right] = Final.Exc1_i;
        Theta
        ConfI = [left,right]
        length = right - left
    end
    
    function SecondExc
        friendlyPrint('***Welcome to the SECOND exercise***\npress a button to start');
        disp('Analytic Solution: ')
        Theta = Final.Exc2
    end
    
    function ThirdExc
        friendlyPrint('***Welcome to the THIRD exercise***\npress a button to start');
        disp('Control variates method: ')
        [Theta,left,right] = Final.Exc3;
        Theta
        ConfI = [left,right]
        length = right - left
    end
    
    function FourthExc
        disp('***Welcome to the FOURTH exercise***');
        friendlyPrint('press a button to see result of exercise 4.1');
        disp('Importance Sampling with mu = 1.1:')
        [Theta,left,right] = Final.Exc4_i(1.1);
        Theta
        ConfI = [left,right]
        length = right - left
        friendlyPrint('press a button to see result of exercise 4.2');
        disp('please wait a moment...')
        [optimal_mu,record]  = Final.Exc4_ii()
    end
    
    function FifthExc
        disp('***Welcome to the FIFTH exercise***');
        friendlyPrint('press a button to see result of exercise 5.1');
        disp('Stratified Sampling with m_st = 500')
        [Theta,left,right] = Final.Exc5_i(500);
        Theta
        ConfI = [left,right]
        length = right - left
        friendlyPrint('press a button to see result of exercise 5.2');
        disp('Length of ConfI Stratified Sampling different mst:')
        mst = ['    mst_20' , '   mst_50', '   mst_100', '   mst_200', '   mst_500','   mst_1000', '   mst_2000']
        StD = Final.Exc5_ii()
    end

    function friendlyPrint(message)
        fprintf(strcat('\n',message,'\n\n'));
        pause();
     end
end