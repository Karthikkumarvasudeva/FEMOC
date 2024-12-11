#include<iostream>
#include<vector>

double triangle_area1(int base, int height){
 //double area;
 //area = 0.5*base*height;
 // return area; // or
 return 0.5*base*height;
}

void triangle_area2(int base, int height, double &area){
 area = 0.5*base*height;
}

void square(int length, double &Area1, double &perimeter){
 Area1 = length*length;
 perimeter = 4*length;
}

void print_vec(std::vector<double> vec1){
    for (int i =0 ; i<vec1.size();i++){
     std::cout<<vec1[i]<<std::endl;
    }

}


int main(){
    double area;

    double Area;
    Area = triangle_area1(5,2);
    std::cout<<Area<<std::endl;


    triangle_area2(7,3, area);
    std::cout<<area<<std::endl;

    double Area1, perimeter;
    square(4,Area1, perimeter);
    std::cout<<Area1<<std::endl;
    std::cout<<perimeter<<std::endl;


    std::vector<double> vectest(6,3.);
    print_vec(vectest);






 return 0;
}
