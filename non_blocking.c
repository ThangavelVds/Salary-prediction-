#include <mpi.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#define n 50 
  
float a[] = {9,9,13,15,8,6,13,5,15,6,14,3,5,9,14,5,14,14,13,13,8,7,4,13,15,13,12,5,4,3,8,5,4,9,3,5,10,3,9,12,12,6,14,12,7,6,15,12,4,4};
float b[] = {32646,34043,46039,46790,30346,32498,49915,33129,49876,34799,48182,34579,32096,33596,47256,33759,47134,46207,46765,49019,30394,32728,31897,49098,49190,49002,46755,30422,34779,33567,31720,32810,30277,34274,32275,32170,34959,33768,30330,48682,46075,32353,47932,49652,32924,30380,48870,47663,34849,34849};
float a2[1000],b2[1000];
int main(int argc, char* argv[]) 
{ 
    int pid, np,elements_per_process, n_elements_recieved; 
    
    MPI_Status status;
    MPI_Request request;
    int request_complete = 0;
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
    MPI_Comm_size(MPI_COMM_WORLD, &np); 
    if (pid == 0) { 
        int index, i; 
        elements_per_process = n / np; 
        if (np > 1) { 
            for (i = 1; i < np - 1; i++) { 
                index = i * elements_per_process; 
  
                MPI_Isend(&elements_per_process, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD,&request);  
                MPI_Isend(&a[index], elements_per_process, MPI_FLOAT, i, 0, MPI_COMM_WORLD,&request);
	        MPI_Isend(&b[index],elements_per_process,MPI_FLOAT, i, 0,MPI_COMM_WORLD,&request);  
            } 

            index = i * elements_per_process; 
            int elements_left = n - index; 

            MPI_Isend(&elements_left, 1, MPI_FLOAT, i, 0,  MPI_COMM_WORLD,&request);         
            MPI_Isend(&a[index], elements_left,  MPI_FLOAT, i, 0,  MPI_COMM_WORLD,&request);
	    MPI_Isend(&b[index],elements_per_process,MPI_FLOAT, i, 0,MPI_COMM_WORLD,&request);
	   if (!request_complete)
	      MPI_Test(&request, &request_complete, &status);
	   if (!request_complete)
              MPI_Wait(&request, &status); 
        } 
        float sum1 = 0; 
        float sum2 = 0;
	for (i = 0; i < elements_per_process; i++) 
	{
		sum1 += a[i]; 
 		sum2 += b[i];
	}	
	
        float tmp1; 
	float tmp2;
        for (i = 1; i < np; i++) 
	{ 
            MPI_Irecv(&tmp1, 1, MPI_FLOAT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD, &request);           
	    MPI_Irecv(&tmp2, 1, MPI_FLOAT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&request); 
            MPI_Wait(&request, &status);
	    int sender = status.MPI_SOURCE; 
  
            sum1 += tmp1;
	    sum2 += tmp2;
        } 
  	float avg_1 = sum1/50;
	float avg_2 = sum2/50;
        printf("\n****************************************************************\n");
        printf("\nAverage of array 1 is : %f\n", avg_1);
        printf("\nAverage of array 2 is : %f\n", avg_2);
	printf("\n****************************************************************\n");
        float top,down,m,c,rmse=0,y_pred;
        float ss_tot=0,ss_res=0,r2;
        for(i=0;i<n;i++)
        {
          top += (a[i] - avg_1) * (b[i] - avg_2);
          down += pow((a[i] - avg_1), 2);
          m = top / down;
          c = avg_2 - (m * avg_1);
        }
	
        printf("\nThe coefficients are:\n%f\t%f\n",m,c);
        for(i=0;i<n;i++)
         {
           y_pred = c + m * a[i];
           rmse += pow((b[i] - y_pred) , 2);
         }
        printf("\nThe RMSE value is :%f\n",sqrt(rmse/n));
        for(i=0;i<n;i++)
        {
          y_pred = c + m * a[i];
          ss_tot += pow((b[i] - avg_2), 2);
          ss_res += pow((b[i] - y_pred) ,2);
        }
    r2 = 1 - (ss_res/ss_tot);
    printf("\nR-Square Value :%f\n",r2);
    printf("\n****************************************************************\n");
	}
        
 else { 
        MPI_Irecv(&n_elements_recieved,1, MPI_INT, 0, 0,MPI_COMM_WORLD,&request); 
        MPI_Irecv(&a2, n_elements_recieved,MPI_INT, 0, 0,MPI_COMM_WORLD,&request); 
  	MPI_Irecv(&b2, n_elements_recieved,MPI_INT, 0, 0,MPI_COMM_WORLD,&request);
        MPI_Wait(&request, &status);	
        // calculates its partial sum 
        float partial_sum_1 = 0;
	float partial_sum_2 = 0;
        for (int i = 0; i < n_elements_recieved; i++) 
	{
            partial_sum_1 += a2[i];
	    partial_sum_2 += b2[i]; 
	}
        MPI_Isend(&partial_sum_1, 1, MPI_FLOAT,0, 0, MPI_COMM_WORLD,&request); 
	MPI_Isend(&partial_sum_2, 1, MPI_FLOAT,0, 0, MPI_COMM_WORLD,&request);
    }    
    MPI_Finalize(); 
  
    exit (0); 
}
