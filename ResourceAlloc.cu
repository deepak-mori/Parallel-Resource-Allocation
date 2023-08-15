#include <iostream>
#include <stdio.h>
#include <cuda.h>

#define max_N 100000
#define max_P 30
#define BLOCKSIZE 1024
#define SLOTS 24

using namespace std;

//*******************************************
// Write down the kernels here
__global__ void find_success(int N, int *d_req_id, int *d_req_fac, int *d_req_cen, int *d_req_start, int *d_req_slots, int *d_capacity, int *d_center_start, int *d_facility, int *d_tot_req, int *d_succ_req){

  int id = blockIdx.x * blockDim.x + threadIdx.x; // Finding id

  if(id<N){
    int cnt_capacity[max_P*SLOTS]; // array that checks maximum capacity for the facility
    int k = 0;
    // Initializing the cnt_capacity array with zero
    while(k < d_facility[id]*SLOTS){
      cnt_capacity[k] = 0;
      k++;
    }

    // Finding request is successfull or not
    int i = d_tot_req[id];
    while(i<d_tot_req[id+1]){
      int s_index = d_req_start[i] - 1 + d_req_fac[i]*SLOTS;
      int e_index = d_req_slots[i] + s_index;
      int l = s_index;
      int cnt = 0;

      while(l < e_index){
        if(cnt_capacity[l] < d_capacity[d_center_start[id]+d_req_fac[i]] && e_index <= d_req_fac[i]*SLOTS + SLOTS){
          cnt++;
        }
        l++;
      }

      if(cnt == e_index - s_index){
        l = s_index;
        while(l < e_index){
          cnt_capacity[l]++;
          l++;
        }
        d_succ_req[id]++;
      }
      i++;
    }
  }
}

//***********************************************

// Function to sort on center ids
int partition_array(int *req_id, int *req_cen, int *req_fac, int *req_start, int *req_slots, int S, int R){

  int i = (S-1);
  for (int j=S; j<=R-1; j++)  
  {  
    if ((req_cen[j] < req_cen[R]) || ((req_cen[j] == req_cen[R]) && req_id[j] < req_id[R])){  
      i++; 
      swap(req_cen[j], req_cen[i]);
      swap(req_fac[j], req_fac[i]);
      swap(req_id[j], req_id[i]);
      swap(req_slots[j], req_slots[i]);
      swap(req_start[j], req_start[i]); 
    }  
  }  
  swap(req_cen[R], req_cen[i+1]);
  swap(req_fac[R], req_fac[i+1]);
  swap(req_id[R], req_id[i+1]);
  swap(req_slots[R], req_slots[i+1]);
  swap(req_start[R], req_start[i+1]); 
  return (i+1); 
} 

// Sorting request array wrt req_cen
void sort_center(int *req_id, int *req_cen, int *req_fac, int *req_start, int *req_slots, int S, int R){
  if (S<R) {  
    int p = partition_array(req_id, req_cen, req_fac, req_start, req_slots, S, R);
    sort_center(req_id, req_cen, req_fac, req_start, req_slots, S, p - 1);  
    sort_center(req_id, req_cen, req_fac, req_start, req_slots, p + 1, R);  
  }  
}

int main(int argc,char **argv)
{
	// variable declarations...
    int N,*centre,*facility,*capacity,*fac_ids, *succ_reqs, *tot_reqs;

    FILE *inputfilepointer;
    
    //File Opening for read
    char *inputfilename = argv[1];
    inputfilepointer = fopen( inputfilename , "r");

    if ( inputfilepointer == NULL )  {
        printf( "input.txt file failed to open." );
        return 0; 
    }

    fscanf( inputfilepointer, "%d", &N ); // N is number of centres
	
    // Allocate memory on cpu
    centre=(int*)malloc(N * sizeof (int));  // Computer  centre numbers
    facility=(int*)malloc(N * sizeof (int));  // Number of facilities in each computer centre
    fac_ids=(int*)malloc(max_P * N  * sizeof (int));  // Facility room numbers of each computer centre
    capacity=(int*)malloc(max_P * N * sizeof (int));  // stores capacities of each facility for every computer centre 

    int success = 0;  // total successful requests
    int fail = 0;   // total failed requests
    tot_reqs = (int *)malloc(N*sizeof(int));   // total requests for each centre
    succ_reqs = (int *)malloc(N*sizeof(int)); // total successful requests for each centre

    // Input the computer centres data
    int k1=0 , k2 = 0;
    for(int i=0;i<N;i++)
    {
      fscanf( inputfilepointer, "%d", &centre[i] );
      fscanf( inputfilepointer, "%d", &facility[i] );
      
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &fac_ids[k1] );
        k1++;
      }
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &capacity[k2]);
        k2++;     
      }
    }

    // variable declarations
    int *req_id, *req_cen, *req_fac, *req_start, *req_slots;   // Number of slots requested for every request
    
    // Allocate memory on CPU 
	  int R;
	  fscanf( inputfilepointer, "%d", &R); // Total requests
    req_id = (int *) malloc ( (R) * sizeof (int) );  // Request ids
    req_cen = (int *) malloc ( (R) * sizeof (int) );  // Requested computer centre
    req_fac = (int *) malloc ( (R) * sizeof (int) );  // Requested facility
    req_start = (int *) malloc ( (R) * sizeof (int) );  // Start slot of every request
    req_slots = (int *) malloc ( (R) * sizeof (int) );   // Number of slots requested for every request
    
    // Input the user request data
    for(int j = 0; j < R; j++)
    {
       fscanf( inputfilepointer, "%d", &req_id[j]);
       fscanf( inputfilepointer, "%d", &req_cen[j]);
       fscanf( inputfilepointer, "%d", &req_fac[j]);
       fscanf( inputfilepointer, "%d", &req_start[j]);
       fscanf( inputfilepointer, "%d", &req_slots[j]);
       tot_reqs[req_cen[j]]+=1;  
    }
		
    // call for sorting the array wrt center id
    sort_center(req_id, req_cen, req_fac, req_start, req_slots, 0, R-1);
    
    // variable declarations on CPU
    int *center_start;
    int *total_req;
    // variable declarations on GPU
    int *d_req_id, *d_req_fac, *d_req_cen, *d_req_start, *d_req_slots, *d_capacity, *d_center_start, *d_facility, *d_tot_req, *d_succ_req;

    // Allocate memory on CPU
    center_start = (int *)malloc(N*sizeof(int));
    total_req = (int *)malloc((N+1)*sizeof(int));

    //Finding starting index of each center and starting index of request of each center
    int sum = 0;
    total_req[0] = 0;
    for(int i=0; i<N; i++){
      total_req[i+1] = total_req[i] + tot_reqs[i];
      center_start[i] = sum;
      sum += facility[i];
    }

    // Allocate memory on GPU
    cudaMalloc(&d_req_id, R*sizeof(int));
    cudaMalloc(&d_req_fac, R*sizeof(int));
    cudaMalloc(&d_req_cen, R*sizeof(int));
    cudaMalloc(&d_req_start, R*sizeof(int));
    cudaMalloc(&d_req_slots, R*sizeof(int));
    cudaMalloc(&d_capacity, max_P*N*sizeof(int));
    cudaMalloc(&d_center_start, N*sizeof(int));
    cudaMalloc(&d_facility, N*sizeof(int));
    cudaMalloc(&d_tot_req, (N+1)*sizeof(int));
    cudaMalloc(&d_succ_req, N*sizeof(int));

    // copying memory from CPU to GPU
    cudaMemcpy(d_req_id, req_id, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_fac, req_fac, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_cen, req_cen, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_start, req_start, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_slots, req_slots, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_capacity, capacity, max_P*N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_center_start, center_start, N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_facility, facility, N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tot_req, total_req, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_succ_req, succ_reqs, N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemset(d_succ_req, 0, N*sizeof(int));


    //*********************************
    // Call the kernels here
    int blocksize = ceil((float)N/1024);
    find_success<<<blocksize, BLOCKSIZE>>>(N, d_req_id, d_req_fac, d_req_cen, d_req_start, d_req_slots, d_capacity, d_center_start, d_facility, d_tot_req, d_succ_req);
    cudaMemcpy(succ_reqs, d_succ_req, N*sizeof(int), cudaMemcpyDeviceToHost);
    //********************************

    // Calculating total number of success and failure
    for(int i=0; i<N; i++){
      success += succ_reqs[i];
    }
    fail = total_req[N] - success;

    // Output
    char *outputfilename = argv[2]; 
    FILE *outputfilepointer;
    outputfilepointer = fopen(outputfilename,"w");

    fprintf( outputfilepointer, "%d %d\n", success, fail);
    for(int j = 0; j < N; j++)
    {
        fprintf( outputfilepointer, "%d %d\n", succ_reqs[j], tot_reqs[j]-succ_reqs[j]);
    }
    fclose( inputfilepointer );
    fclose( outputfilepointer );
    cudaDeviceSynchronize();
	return 0;
}