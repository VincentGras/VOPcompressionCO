function VOP=VOP_compression_Lee_SOR(filename_in,eps_G,Sglobal,name_VOP,use_Sglobal,show_waitbar)
% This code was written by Stephan Orzada, German Cancer Research Center (DKFZ).
% In 2020; email: stephan.orzada@dkfz.de
% VOP-compression algorithm according to Lee et al. (DOI 10.1002/mrm.23140)
% This alogrithm uses a pre-determined Sglobal matrix for enhanced
% performance at low SAR values. Needs R2016b or newer.
% eps_G: Overestimation control according to Lee.
% filename_in: Points to a .mat file containting a variable "matrices" containing SAR matrices (Q-matrices); (N_channel, N_channel, N_matrices).
% Sglobal: Provide Sglobal for algorithm (N_channel,N_channel). If you do not use it, provide
% anything and set "use_Sglobal" to 1 (for Sglobal according to Lee) or 2 for diagonal matrix.
% name_VOP: provide a string for a name for save-file.
% use_Sglobal: set to 0 to use your own Sglobal matrix, set to 1 to
% calculate it according to Lee.
% show_waitbar: set to 1 to see a waitbar, set to 0 to disable waitbar.
% VOP contains the calculated VOPs (N_channel, N_channel, N_VOPs).
% The results are also saved in a file named "VOP_SOR_'name_VOP'_'eps_G'"
%
% This script needs R2016b or later.
% This script is heavily optimized for speed. It uses parallel computing and a
% lot of vecorization. SAR matrices are explicitly made symmetric within the code
% to enhance the eigenvalue decomposition. The coefficient update is implemented in a non-intuitive
% way to speed up the calculations. Please refer to the corresponding
% comments in the code.
% The code uses Cholesky-Factorization as proposed by Andre Kuehne (Proc.
% Intl. Soc. ISMRM 2017, #0478.)
% Most of this code was used in Orzada et al. 2020: https://doi.org/10.1007/s10334-020-00890-0
% This version does not use the iterative approach as proposed by Lee, but
% iterates to make the minimum Eigenvalue of P (see below) >=0. This
% converges much better.
%
%  IMPORTANT NOTE:
%  The author takes no responsibility for the validity of the SAR models
%  generated with this tool. Validating the results is the responsibility
%  of the user!

tic %Start timer for total Elapsed time.
%constants for internal purpose:
Num_iterations=500; %Iterations for upperbound determination loop
Num_inner_iterations=100; %Iterations for updating Cwv.
update_waitbar_iter=1000; %Update waitbar after checking x matrices.
convergence_scheme=0; %how to update c_wv. 0: multiply with random value (slow, better compression, recommended), 1:insert random value (fast). 


disp('Loading matrices...')
load(filename_in,'matrices');
disp('Loading complete.')
try
[x,~,N]=size(matrices); %Determine size of Omega
catch
    error('File does not contain a variable named "matrices"!')
end


if use_Sglobal==1 %If the user wants to use the actual Sglobal, calculate it
    Sglobal=sum(matrices,3);
elseif use_Sglobal==2 %If user wants a diagonal matrix, calculate it.
    Sglobal=diag(ones(x,1));
end

eigenvalues_max=zeros(N,1);

disp('Calculating all Eigenvalues...')
parfor a=1:N %Using parfor is faster than using for in this case. If the user cannot use the parallelized form, "for" can be used instead. This is only a minor reduction in speed for large problems, since this part is only calculated ones.
    matrices(:,:,a)=(matrices(:,:,a)+matrices(:,:,a)')/2; %Make Sv symmetric. Due to numerical errors, Sv can be slightly non-symmetric, which causes eig to be slower and give non-real eigenvalues.
    eig_val_temp=eig(matrices(:,:,a));%Calculate maximum Eigenvalues
    eigenvalues_max(a)=max(real(eig_val_temp)); %Calculate maximum of real part.
end

disp('Sorting Eigenvalues...')
[B,I_eigen]=sort(eigenvalues_max,'descend'); %Sort maximum Eigenvalues. I_eigen contains the indexes.
disp('Reordering matrices...')
matrices=matrices(:,:,I_eigen); %Resorting all matrices enhances speed when fetching matrices for cholesky factorization. Probably due to prefetch. You can omit this line and the following one if you run into memory issues.
I_eigen=1:N; %The new order is now 1,2,3,...,N

Sglobal=(Sglobal+Sglobal')/2; %make sure Sglobal is symmetric. (faster eig and only real eigenvalues).

Sglobal=Sglobal/min(real(eig(Sglobal)))*B(1); %Normalize Sglobal. We use the minimum eigenvalue here, because it determines the minimum value the overestimation term can assume, and this is crucial for the number of VOPs.

Sglobal=Sglobal*eps_G; %multiply with eps_G for overestimation control.

overestimation_factor=max(real(eig(Sglobal)))/real(B(1))*100; %Calculate the maximum overestimation for information purposes.

disp(['Starting: ' num2str(eps_G) ' Maximum Overestimation: ' num2str(overestimation_factor) '% of Worst Case SAR'])

clear B %not needed anymore. We only need I, since it contains the order in which we look at the matrices Sv.



V_sub=matrices(:,:,I_eigen(1))+Sglobal; %Put first Voxel in subvolume and add Sglobal for Overestimation

if show_waitbar
    f = waitbar(0,'Starting');
end

ran_out_counter=0; %Counter for matrices which are included as VOPs without definitive descision.
num_vops=1; %initialize num vops
N_start=N; %Save number of matrices.
a=1;
while a<N %Evaluate all matrices (voxels). 1st has already been put in.
    a=a+1;
    num_vops_new=size(V_sub,3); %Number of vops this round.
    if num_vops_new>num_vops || a==2 %Check whether number of VOPs has changed. Allocate Memory for SAR computation
        num_vops=num_vops_new;
        V_sub_pre=V_sub/num_vops;
        P_pre=sum(V_sub_pre,3); %pre-calculate because it is the same for the first iteration for EVERY matrix.
    end
    can_be_upperbounded=0;

    %1. Initialize Coefficients
    c_wv=ones(1,1,num_vops)/num_vops; %use all ones as starting vector.

    Sv_current=matrices(:,:,I_eigen(a)); %Get current SAR matrix
    
    P=P_pre-Sv_current; %Calculate P for the first iteration step. For later steps, the calculation is done under point 4.
    ax_count=0;
    V_sub_temp1=V_sub_pre;
    for b=1:Num_iterations
        
        %2. Can it be upperbounded?
        
        [V,D]=eig(P,'vector');   %Calculate eigenvalues and -vectors
        V_eigenvalues=real(D); %take real part for calculation of inequalities and minimum.
        eigen_vector_analysis=sum(V_eigenvalues<0); %Find out if there is a negative Eigenvalue.

        if eigen_vector_analysis == 0 %if eigenvalues are positive, it can be upperbounded. Break.
            can_be_upperbounded=1; %it can be upperbounded, so it does not have to be put in Subset.
            kappa=a;
            N_save=N;
            try
                P_pre2=sum_V_sub_temp2/sum(C_WV_temp); %May not exist at first VOP.
            catch
                P_pre2=P_pre;
            end
            delete_pos=zeros(N,1); %This is an array to define which matrices are dominated and can be taken out of the calculation.
            while kappa<N
                kappa=kappa+1;
                Sv_current=matrices(:,:,I_eigen(kappa));
                P=P_pre2-Sv_current;  %If this is Positive Semi Definite, the matrix Sv_current is dominated
                [~,flag] = chol(P); %We only need to know whether a Cholesky decomposition is possible, we do not need the actual result.
                if flag==0
                    delete_pos(kappa)=1;
                end
            end
            dropped_mats=sum(delete_pos);
            I_eigen(logical(delete_pos))=[]; %Delete the entries from this vector. Only matrices which this vector still points to are included in further calculations.
            N=numel(I_eigen); %Determine new number of matrices we are still working with.
            disp(['Dropped ' num2str(dropped_mats) ' of ' num2str(N_save) ' matrices using coefficients from last step.'])
            break
        end
        
        %3. Can't it be upperbounded?
        [~,I]=min(V_eigenvalues); %find minimum eigenvalue
        V_eigen_min=V(:,I); %Select eigenvector corresponding to minimal eigenvalue
                
        %Calculate SAR in subset. This is done in a vectorized way to enhance speed.        
        A_2=V_sub.*conj(V_eigen_min); %Multiply slices with scalars (1st step of vector-matrix multiplication) (>=R2016b)
        A_3=sum(A_2,1); %Add up the slices (2nd step of vector-matrix multiplication)
        A_4=A_3.*V_eigen_min.'; %Multiply vectors with scalars (Works from R2016b onward)
        A_5=sum(A_4,2); %add up the vectors.
        result_SAR=max(abs(A_5)); %find maximum value over all VOPs.
        
        if abs(V_eigen_min'*Sv_current*V_eigen_min)>result_SAR
            %It can't be upperbounded, so break. 'can_be_upperbounded' remains 0
            break
        end
        
        %4. Update Coefficients c_wv (It is not clear in Lee's Paper how this
        %is done. This is my heuristic implementation which works quite well.
        notnegative=0; 
        test_value_temp=-1e100; %intialize with very low value
        P_temp=P;
        success_counter=0;
        for ax=1:Num_inner_iterations
            C_WV_temp=c_wv; %prepare temporary c_wv
            V_sub_temp2=V_sub_temp1; %save V_sub_temp1
            pos_change=ceil(num_vops*rand());
            if convergence_scheme==0
                C_WV_temp(pos_change)=C_WV_temp(pos_change)*abs(randn); %change one value. Fun fact: ceil(num_vops*rand()) is about 10 times faster than randi(num_vops). (R2018b).
            else
                C_WV_temp(pos_change)=abs(randn)*sum(c_wv);
            end
            %after execution of the for loop, normalization of c_wv and V_sub_temp is necessary.
            V_sub_temp2(:,:,pos_change)=V_sub(:,:,pos_change).*C_WV_temp(pos_change); %Only apply the CHANGED coefficient. Significantly reduces the number of multiplications. The necessary normalization is done later.
            sum_V_sub_temp2=sum(V_sub_temp2,3);
            P=sum_V_sub_temp2/sum(C_WV_temp)-Sv_current; %Do the normalization here.
            test_value=min(real(eig(P))); %Instead of Lee's approach, try to make the minimum Eigenvalue larger. Converges more quickly.
            if test_value>=0    %if test value is non-negative, break and tell outer loop
                success_counter=success_counter+1;
                if success_counter==1000 %Extra iterations help drop more matrices during the cholesky-step.
                    c_wv=C_WV_temp; %c_wv can be updated, since a suitable result has been found.
                    V_sub_temp1=V_sub_temp2; %We need to save V_sub, since the coefficients have been applied already.
                    notnegative=1; %This flag can be used to display that a suitable result has been found.
                    P_temp=P;
                    break
                end
            end
            if test_value>=test_value_temp %if test_value is larger than test_value of previous iteration, keep value, if not, it is discarded.
                c_wv=C_WV_temp; %update c_wv temporarily. c_wv can be updated, since the new coefficients are better than the old ones.
                V_sub_temp1=V_sub_temp2; %We need to save V_sub, since the coefficients have been applied already.
                test_value_temp=test_value;  %keep test_value in mind.
                P_temp=P;
            end
        end
        P=P_temp;
        ax_count=ax_count+ax; %cumulatively save the number of inner loop counts, just to inform the user.
        if notnegative==0
            %disp('still negative'); %The used might like to know that no suitable result has been reached. This is not a problem actually, but the user might like to be informed.
        end
        
        %do the normalization to prevent overflow:
        sum_cwv=sum(c_wv);
        if sum_cwv>1e3 %As long as the sum is not too big, we do not need to do the normalization, since no overflow will occur.
            V_sub_temp1=V_sub_temp1/sum_cwv;
            c_wv=c_wv/sum_cwv;
        end
        %Repeat until iterations are all done. Should only take a few
        %iterations.
    end
    
    if can_be_upperbounded==0 %If broken out of loop and not upperbounded, add it to subset.
        if b==Num_iterations %if number of iterations was too low, display that
            disp('Ran out of iterations');
            ran_out_counter=ran_out_counter+1;
        end
        if show_waitbar
            waitbar((N_start-N+a)/N_start,f,[name_VOP ', ' num2str((N_start-N+a)/N_start*100,'%4.2f') ' %. #VOPs: ' num2str(num_vops)]); %Update Waitbar.
        end
        V_sub(:,:,num_vops+1)=Sv_current+Sglobal; %Add current matrix to subset.
        disp([num2str(eps_G) ', Finished: ' num2str((N_start-N+a)/N_start*100) '% VOP: ' num2str(num_vops+1) ' #Iterations: ' num2str(b) ' avg. inner loops: ' num2str(floor(ax_count/b))])
    end
    %disp(['Number of iterations: ' num2str(b) ' Inner loops: ' num2str(ax_count) ' Avg. inner loops: ' num2str(ax_count/b)])
    if mod(a,update_waitbar_iter)==0 && show_waitbar
        waitbar((N_start-N)/N_start,f,[name_VOP ', ' num2str((N_start-N)/N_start*100,'%4.2f') ' %. #VOPs: ' num2str(num_vops)]); %Update Waitbar.
    end
    if mod(a,1000)==0 %give a progress every 1000 checked matrices
        disp(['Progress eps_G = ' num2str(eps_G) '. Finished: ' num2str((N_start-N+a)/N_start*100) '% #VOPs: ' num2str(num_vops)])
    end
end
elapsed_time=toc;
disp([num2str(ran_out_counter/num_vops*100) '% VOPs with insufficient iterations for descision.'])
disp(['Finished with ' num2str(num_vops) ' VOPs. Elapsed time is: ' num2str(uint32(elapsed_time)) ' s.'])

VOP=V_sub; %place all matrices in VOP variable

for a=1:num_vops
    VOP(:,:,a)=V_sub(:,:,a); %call SAR-matrices VOPs for saving them
end

if show_waitbar
    close(f)
end

save(['VOP_Lee_SOR_' name_VOP '_' num2str(eps_G) '.mat'],'VOP','Sglobal','elapsed_time','eps_G','ran_out_counter') %save VOPs in file with a name that tells what was compressed and how.