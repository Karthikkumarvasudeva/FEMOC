#include<iostream>
#include<vector>
#include "triangle.h"



int main(){

Triangle tri1(1.,8), tri2(3,2);

std::cout<<tri1.area()<<std::endl;

tri1.base = 4.;

std::cout<<tri1.area()<<std::endl;


 return 0;
}
