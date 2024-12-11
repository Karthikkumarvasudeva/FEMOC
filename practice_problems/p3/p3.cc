#include<iostream>
#include<vector>

int main(){

int a, *b; //pointer b
a=3;
b = &a;
std::cout <<*b <<std::endl;

a = 4;
std::cout <<*b <<std::endl;



std::vector<double> vec1(2,3.1), *vec2;

vec2 = &vec1;

std::cout<<vec1.size()<<" "<<vec2->size()<<std::endl;

vec1.resize(3,2.1);

std::cout<<vec1.size()<<" "<<vec2->size()<<std::endl;


std::vector<double> vec(2,7.77);

for(std::vector<double>::iterator  it = vec.begin(); it != vec.end(); it++){ //iterator it = vector first element; it not equal to vector last element (doesnt run out) ; increment to next element)

 std::cout<< *it <<std::endl;

}





 return 0;
}
