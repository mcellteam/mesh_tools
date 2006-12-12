
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

main(){

	int j;
	int k;
	for (j=0;j<3;j++) {
		k = (j+1)%3; 
		cout << j << " " << k << endl;
	}

}

