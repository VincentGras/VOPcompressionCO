function [fileout]=VOP_compression_iterative_VG(filename_in,eps_G,R,Sglobal,name_VOP,use_Sglobal,max_num_VOPs, min_epsG)
% VG : Minor modification for the purpose of this comparative study 
% VG : (see VOP_compression_iterative for comparison with the original)
% This code was written by Stephan Orzada, German Cancer Research Center (DKFZ).
% October 2020 
% email: stephan.orzada@dkfz.de
% Enhanced VOP-Compression algorithm (Paper under Review).
% This compression algorithm uses an iterative approach to reduce the
% number of VOPs and the calculation time for a given overestimation. 
% The VOP-compression algorithm according to Lee et al. (DOI
% 10.1002/mrm.23140) is run in an iterative way, each iteration starting
% with the sub set (VOPs) derived from the previous iteration and a
% decreased overestimation compared to the previous iteration.
% Needs R2016b or newer.
%
% filename_in: Points to a .mat file containting a variable "matrices" containing SAR matrices (Q-matrices); (N_channel, N_channel, N_matrices).
% eps_G: Starting value for Overestimation, relative to worst case local SAR. (Note: The MINIMUM value of the overestimaton term is eps_G*SAR_wc)
% R: Factor to reduce the overestimation after each iteration step. (The overestimation factor eps_G is multiplied by R after each step)
% Sglobal: Provide Sglobal for algorithm (N_channel,N_channel). If you do not use it, provide anything and set "use_Sglobal" to 1 (for Sglobal according to Lee) or 2 for diagonal matrix.
% name_VOP: provide a string for a name for save-file.
% use_Sglobal: set to 0 to use your own Sglobal matrix, set to 1 to calculate it according to Lee.
% max_num_VOPs: maximum number of VOPs. At this number of VOPs, the algorithm is stopped.
% 
% VOP contains the calculated VOPs (N_channel, N_channel, N_VOPs) including overestimation.
% The results are saved after each iteration in a file named "VOP_SOR_'name_VOP'_'eps_G'.mat". Here eps_G is the overestimation factor of respective iteration.
%
% This script needs R2016b or later.
% This script is heavily optimized for speed in Matlab. It uses parallel computing and a
% lot of vecorization. SAR matrices are explicitly made symmetric within the code
% to enhance the eigenvalue decomposition. The coefficient update is implemented in a non-intuitive
% way to speed up the calculations. Please refer to the corresponding
% comments in the code.
% The code uses Cholesky-decomposition as proposed by Andre Kuehne (Proc.
% Intl. Soc. ISMRM 2017, #0478.)
% Sections of this code are derived from Orzada et al. 2020: https://doi.org/10.1007/s10334-020-00890-0
% To find the coefficient c_wv this version does not use the iterative approach as proposed by Lee, but
% iterates to make the minimum Eigenvalue of P (see below) >=0. This
% converges much better and yields better results for the
% Cholesky-decomposition part proposed by Kuehne.
%
%  IMPORTANT NOTE:
%  The author takes no responsibility for the validity of the SAR models
%  generated with this tool. Validating the results is the responsibility
%  of the user!



tic %Start timer for total Elapsed time.

fileout = cell(0,1);

%constants for internal purpose:
Num_iterations=500; %Iterations for upperbound determination loop
Num_inner_iterations=100; %Iterations for updating Cwv.
convergence_scheme=0; %how to update c_wv. 0: multiply with random value (slow, better compression, recommended), 1:insert random value (fast).
eps_decrease=R;

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
parfor a=1:N %Using parfor is faster than using for in this case. If the user cannot use the parallelized form, "for" can be used instead. This is only a minor reduction in speed for large problems, since this part is only calculated once.
    matrices(:,:,a)=(matrices(:,:,a)+matrices(:,:,a)')/2; %Make Sv symmetric. Due to numerical errors, Sv can be slightly non-symmetric, which causes eig to be slower and give non-real eigenvalues.
    eig_val_temp=eig(matrices(:,:,a)); %Calculate maximum Eigenvalues
    eigenvalues_max(a)=max(real(eig_val_temp)); %Calculate maximum of real part.
end

disp('Sorting Eigenvalues...')
[B,I_eigen]=sort(eigenvalues_max,'descend'); %Sort maximum Eigenvalues. I_eigen contains the indexes.
disp('Reordering matrices...')
matrices=matrices(:,:,I_eigen); %Resorting all matrices enhances speed when fetching matrices for cholesky factorization. Probably due to prefetch. You can omit this line and the following one if you run into memory issues.
I_eigen=1:N; %The new order is now 1,2,3,...,N
I_eigen_save=I_eigen; %save I_eigen for later runs.
N_start=N; %Save number of matrices.

Sglobal=(Sglobal+Sglobal')/2; %make sure Sglobal is symmetric. (faster eig() and only real eigenvalues).

Sglobal=Sglobal/min(real(eig(Sglobal)))*B(1); %Normalize Sglobal. We use the minimum eigenvalue here, because it determines the minimum value the overestimation term can assume, and this is crucial for the number of VOPs.

Sglobal=Sglobal*eps_G; %multiply with eps_G for overestimation control.

overestimation_factor=max(real(eig(Sglobal)))/real(B(1))*100; %Calculate the maximum overestimation for information purposes.

disp(['Starting: ' num2str(eps_G) ' Maximum Overestimation: ' num2str(overestimation_factor) '% of Worst Case SAR'])

clear B %not needed anymore. We only need I, since it contains the order in which we look at the matrices Sv.



V_sub=matrices(:,:,I_eigen(1))+Sglobal; %Put first Voxel in subvolume and add Sglobal for Overestimation

ran_out_counter=0; %Counter for matrices which are included as VOPs without definitive descision.
num_vops=1; %initialize num vops
run_number=1;
too_many_vops=false; %to check whether a run has finished correctly or was broken due to too many vops.
while num_vops<max_num_VOPs && eps_G >= min_epsG
    I_eigen=I_eigen_save; %Since this is a new iteration step, we need all matrices again.
    
    V_sub=V_sub-Sglobal; %Due to memory and calculation time constraints, the sub set does not contain the original matrices as described in Lee's paper, but the matrices + overestimation. Here we recover the original matrices.
    if run_number>1 %If this is not the first iteration step, reduce the overestimation
        Sglobal=Sglobal*eps_decrease; %Change the overestimation term
        eps_G=eps_G*eps_decrease; %Change Overestimation factor (this is only used for display reasons hereafter).
    end
    disp(['Starting run #' num2str(run_number) '. eps_G = ' num2str(eps_G)])
    run_number=run_number+1;
    V_sub=V_sub+Sglobal; %Add the overestimation term to the matrices in the sub set. This saves calculation time, because it is only done once for each matrix in each iteration step.
    a=1; %Counter for the matrix number
    N=N_start; %Total number of matrices in the full set.
    
    while a<N %Evaluate all matrices (voxels). 1st has already been put in.
        a=a+1; %select next matrix.
        num_vops_new=size(V_sub,3); %Number of vops this round.
        if num_vops_new>num_vops || a==2 %Check whether number of VOPs has changed. Allocate Memory for SAR computation
            num_vops=num_vops_new;
            V_sub_pre=V_sub/num_vops;
            P_pre=sum(V_sub_pre,3); %pre-calculate because it is the same for the first iteration for EVERY matrix.
            steps_since_last_VOP=0;
        end
        steps_since_last_VOP=steps_since_last_VOP+1;
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
                kappa=a; %Initialize counter for checking other matrices with Cholesky.
                N_save=N;
                try
                    P_pre2=sum_V_sub_temp2/sum(C_WV_temp); %May not exist at first. VOP.
                catch
                    P_pre2=P_pre; %is the sum of all VOPs multiplied by the coefficients c_wv.
                end
                delete_pos=zeros(N,1); %This is an array to define which matrices are dominated and can be taken out of the calculation.
                while kappa<N
                    kappa=kappa+1;
                    Sv_current=matrices(:,:,I_eigen(kappa));
                    P=P_pre2-Sv_current; %If this is Positive Semi Definite, the matrix Sv_current is dominated
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
            
            %do the normalization to prevent overflow:
            sum_cwv=sum(c_wv);
            if sum_cwv>1e3 %As long as the sum is not too big, we do not need to do the normalization, since no overflow will occur.
                V_sub_temp1=V_sub_temp1/sum_cwv;
                c_wv=c_wv/sum_cwv;
            end
            %Repeat until iterations are all done.
        end
        
        if can_be_upperbounded==0 %If broken out of loop and not upperbounded, add it to subset.
            if b==Num_iterations %if number of iterations was too low, display that
                disp('Ran out of iterations');
                ran_out_counter=ran_out_counter+1; %Count the number of times the algorithm couldn't decide whether a matrix is a VOP or not. Purely for information purpose.
            end
            V_sub(:,:,num_vops+1)=Sv_current+Sglobal; %Add current matrix to subset.
            disp([num2str(eps_G) ', Finished: ' num2str((N_start-N+a)/N_start*100) '% VOP: ' num2str(num_vops+1) ' #Iterations: ' num2str(b) ' avg. inner loops: ' num2str(floor(ax_count/b))])
            
            if num_vops==max_num_VOPs
                too_many_vops = true;
                break
            end
            
            
            
        end
    end
    elapsed_time=toc;
        
    VOP=V_sub; %place all matrices in VOP variable

    disp([num2str(ran_out_counter/num_vops*100) '% VOPs with insufficient iterations for descision.'])
    fileout = cat(1, fileout, fullfile(fileparts(filename_in), ['VOP_SOR_' name_VOP '_' num2str(eps_G) '.mat']));
    disp(['VOP_SOR_' name_VOP '_' num2str(eps_G) '.mat'])
    pause(1);
    save(fileout{end},'VOP','Sglobal','elapsed_time','eps_G','ran_out_counter') %save VOPs in file with a name that tells what was compressed and how.
    disp('Too many VOPs. Finishing.')
    disp(['Finished with ' num2str(size(VOP,3)) ' VOPs. Elapsed time is: ' num2str(uint32(elapsed_time)) ' s.'])
end
