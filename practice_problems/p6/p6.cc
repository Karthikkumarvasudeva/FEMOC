#include<iostream>
#include<vector>
#include "triangle.h"



int main(){

Triangle tri1(1.,5), tri2(3,2);

std::cout<<tri1.area()<<std::endl;

tri1.base = 3.;

std::cout<<tri1.area()<<std::endl;


 return 0;
}
