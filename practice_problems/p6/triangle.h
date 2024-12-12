class Triangle{

public:
    Triangle(double b, double h); //constructor
    //~Triangle(){}; //destructor
    double area();
    double base, height;

};

//defining a function to the class

Triangle::Triangle(double b, double h){
 base = b;
 height = h;
}

double Triangle::area(){
 return 0.5*base*height;

}
