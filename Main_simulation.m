% Simulation. Need to call encode.m and decode.m as function.
% Obtain the size of matrix U after the decoding,
% and the decoding time. Results are averaged.
% 3 algorithms are tested.'standard','Kim','random'.
% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------
K_list = [1000:1000:7000, 8192];
% K_list = 200:200:2000;
frames = 50; 
Simu_U = zeros(3, length(K_list));
Simu_time = zeros(3, length(K_list));
mark = {'-ko', '-k*', '-k^','-kv'};

% run 
for jj = 1 : length(K_list)
    K = K_list(jj);
    [G_LT_all, Symbol_all] = encode (K);
    for method_index = 1:3
        tmp_u = [];
        tmp_time = [];
        for iii = 1: frames
            [suc_flag, time_pass, u ] = decode (K, G_LT_all, Symbol_all, method_index);
            if suc_flag
                tmp_u = [tmp_u,u];
                tmp_time = [tmp_time,time_pass];
            end
        end
        Simu_U(method_index,jj) = mean(tmp_u);
        Simu_time(method_index,jj) = mean(tmp_time);
    end
end

% plot
figure(1);
for ii = 1: 3
    plot(K_list,Simu_U(ii,:),mark{ii}); hold on;
end
xlabel('K');ylabel('u');
title('u ' );
legend('standard','Kim','random','location','northwest');
grid on;
figure(2);
for ii = 1: 3
    plot(K_list,Simu_time(ii,:),mark{ii}); hold on;
end
xlabel('K');ylabel('Decoding Time(sec)');
title('time ' );
legend('standard','Kim','random','location','northwest');
grid on;




    
    