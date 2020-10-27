#include <Rcpp.h>
##include <math.h>
using namespace Rcpp;
//' Describe the function
//'
//' @param describe the params
//' 
//' @return describe the return 
// [[Rcpp::export]]
double AfsariCorcpp(NumericVector x, NumericMatrix mydist){
	int n = x.size(),i,j,k, Max, Sec ; // Size of vector
	
	double maptoSec[2][2][2] =
	{{ {2, 3}, {1, 1}, },
	{ {1, 1}, {3, 2} }};
	
	//double maptoMax[2][2][2] = { 
	// 1 -> (i,j) is the furthest -> median k -> 3
	// 2 -> (i,k) is the furthest -> median j -> 2
	// 3 -> (j,k) is the furthest -> median 
	//       { {(i,j), (i,j)}, {(i,j), (j,k)}, },
	//       { {(i,k), (i,j)}, {(i,k), (j,k)} }};
	
	double maptoMedian[2][2][2] = { 
		{ {3, 3}, {3,1}, },
		{ {2, 3}, {2,1} }
	};
	double sum = 0; // Sum value

	// For loop, note cpp index shift to 0
	for(i = 0; i < n; i++)
	for(j = i + 1 ; j < n; j++)
	for(k = j + 1 ; k < n; k++)
	{
		Sec = maptoSec[x[i]<x[j]][x[i]<x[k]][x[j]<x[k]];
		Max = maptoMedian[mydist(i,j)<mydist(i,k)][mydist(i,j)<mydist(j,k)][mydist(i,k)<mydist(j,k)];
		sum += ( Max == Sec  );
		//Rcout << i << j << k << Max << Sec << std::endl;
	}
	return(6.0*sum/n/(n-1)/(n-2)); // Obtain and return the Mean
}


//' Describe the function
//'
//' @param describe the params
//' 
//' @return describe the return 
// [[Rcpp::export]]
double AfsariCorVarcpp(NumericMatrix mydist){
	int n = mydist.nrow(), Max; // Size of vector
	int MedianCount = 0; 
	//Count how many times i the median of both (i,j,k) and (i,l,o)  
	
	
	int ijk[3]; 
	
	//double maptoSec[2][2][2] =
	//        {{ {2, 3}, {1, 1}, },
	//         { {1, 1}, {3, 2} }};
	
	//double maptoMax[2][2][2] = { 
	// 1 -> (i,j) is the furthest -> median k -> 3
	// 2 -> (i,k) is the furthest -> median j -> 2
	// 3 -> (j,k) is the furthest -> median 
	//       { {(i,j), (i,j)}, {(i,j), (j,k)}, },
	//       { {(i,k), (i,j)}, {(i,k), (j,k)} }};
	
	int maptoMedian[2][2][2] = { 
		{ {3, 3}, {3,1}, },
		{ {2, 3}, {2,1} }
	};
	//double sum = 0; // Sum value
	
	
	
	// For loop, note cpp index shift to 0
	for(int i = 0; i < n; i++)
		for(int j = i + 1 ; j < n; j++)
			for(int k = j + 1 ; k < n; k++) {
				ijk[0] = i;
				ijk[1] = j;
				ijk[2] = k;
				
				Max = ijk[maptoMedian[mydist(i,j)<mydist(i,k)]
				[mydist(i,j)<mydist(j,k)]
				[mydist(i,k)<mydist(j,k)]-1];
	
	
				for( int l = 0 ; l < n ; l ++ ){
					if( i == l || j == l || k == l ) continue;
					for( int o = l + 1 ; o < n; o++){
						if( i == o || j == o || k == o ) continue;
	
						if(maptoMedian[mydist(Max,l)<mydist(Max,o)]
							[mydist(Max,l)<mydist(l,o)]
							[mydist(Max,o)<mydist(l,o)] == 1) MedianCount++;
					}
		
				}
			}
	//return(6.0*sum/n/(n-1)/(n-2)); // Obtain and return the Mean
	return((double)MedianCount);
}




//' Describe the function
//'
//' @param describe the params
//' 
//' @return describe the return 
// [[Rcpp::export]]
NumericMatrix AfsariCorVarcpp2(NumericMatrix mydist){
	int n = mydist.nrow(), Max1, Max2; // Size of vector
	
	//int Medianalpha = 0 ,
	//    Medianbeta = 0,
	//    Mediangamma = 0,
	//    Mediantheta = 0;
	
	int MedianCount[3][3] = {{0,0,0},
	{0,0,0},
	{0,0,0}};
	NumericMatrix Out(3,3);
	
	//Counts how many times i the median of both (i,j,k) and (i,j,l)
	//Counts how many times k is the median of (i,j,k) 
	//and l is the median of(i,j,l)
	
	
	
	//double maptoMax[2][2][2] = { 
	// 1 -> (i,j) is the furthest -> median k -> 3
	// 2 -> (i,k) is the furthest -> median j -> 2
	// 3 -> (j,k) is the furthest -> median 
	//       { {(i,j), (i,j)}, {(i,j), (j,k)}, },
	//       { {(i,k), (i,j)}, {(i,k), (j,k)} }};
	
	int maptoMedian[2][2][2] = { 
	{ {2, 2}, {2,0}, }, //{{ {3, 3}, {3,1}, },
	{ {1, 2}, {1,0} }};//{ {2, 3}, {2,1} }};
	
	
	
	// For loop, note cpp index shift to 0
	for(int i = 0; i < n; i++)
	for(int j = 0 ; j < n ; j++) {
		if( i == j)	continue;
	
		for(int k = 0; k < n ; k++) {
			if( (k == i) || (k == j) ) continue;
			//Rcout << i <<\" \" << j << \" \" << 
			//           k << std::endl;
			for(int l = 0 ; l < n ; l++) {
				if( (l == i) || (l == j) || (l == k) ) continue;
	
				Max1 = maptoMedian[mydist(i,j)<mydist(i,k)]
					[mydist(i,j)<mydist(j,k)]
					[mydist(i,k)<mydist(j,k)];
	
				Max2 = maptoMedian[mydist(i,j)<mydist(i,l)]
					[mydist(i,j)<mydist(j,l)]
					[mydist(i,l)<mydist(j,l)];
	
				MedianCount[Max1][Max2]++;
	
				// //Rcout << i <<\" \" << j << \" \" << 
				// //         k <<\" \" << l << \":\" <<
				// //         Max1 <<\" \" << Max2 << std::endl;
				// 
				// if( ( Max1 == Max2 )  && ( Max1 < 2 ))
				//   Medianalpha ++;
				// else if(( Max1 == 2 )  && ( Max2 == 2 ))
				//     Mediantheta ++;
				// else if(( Max1 == 2 ) || ( Max2 == 2 ))
				//     Mediangamma ++;
				// else Medianbeta++;
			}
		}
	}
	for( int i = 0 ; i < 3; i++ )
		for( int j =0; j < 3; j++ )
			Out(i,j) = MedianCount[i][j];
	
	return(Out);
}

