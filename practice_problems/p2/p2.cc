#include<iostream>
#include<vector>


// Overload the << operator to print a std::vector
std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i < vec.size() - 1) os << ", ";
    }
    os << "]";
    return os;
}

int main(){

//std::cout<<"Hello World!"<<std::endl;
 int a=3, b=3;

 // if and else if  conditions

 if (a == 2){
  std::cout<<"the first if loop is executed with a=2"<<std::endl;
 }
 else if(a == 3){

  std::cout<<"if else loop is executed with a=3"<<std::endl;
 }
 else {
  std::cout<<"no conditions are met"<<std::endl;
 }


// if and else if with && and || conditions

 if (a == 2 && b == 3){
  std::cout<<"if with a=2 and b=3"<<std::endl;
 }
 else if(a == 4 || b == 3){

  std::cout<<"else if with a=4 or b=3 "<<std::endl;
 }
 else {
  std::cout<<"no conditions are met"<<std::endl;
 }


 // for condition

 for (int i=0; i<=10; i++){
  std::cout<<i<<std::endl;
 }


std::vector<double> vec(10,2); //vector of size 10 with elements 2

// for loop that adds i value to the vector defined above
for(int i=0 ; i<10; i++){
 vec[i]+=i;
 std::cout<<vec[i]<<std::endl; //vec elements
 std::cout<<vec[i]<<std::endl; //vec with 2's incremented witin for loop
}

int c=4;
std::vector<double> vec1(5,c);
std::cout<<"vector (5,c) with c=4 : "<<vec1<<std::endl;

 return 0;
}
