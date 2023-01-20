/**
*  
* Solution to course project # 3
* Introduction to programming course
* Faculty of Mathematics and Informatics of Sofia University
* Winter semester 2022/2023
*
* @author Ivaylo Andreev
* @idnumber 4MI0600186* @compiler GCC
*
*
*
*/


#include <iostream>
#include <cmath>
#include <cstring>
#pragma warning( disable : 4996)
using namespace std;

const char beg = '(';
const char fin = ')';
const char mid = ',';
const char sp = ' ';
const char root_beg = '(';
const char root_end = ')';
const char over = '/';
const char mult = '*';
const char pl = '+';
const char mi = '-';
const char sq[5] = "sqrt";
const int epsilon = 0.001;

const unsigned int number_of_options = 12;
const unsigned int max_size = 1024;

void Absolute_Value (int& num) {
    if (num < 0) {
        num *= -1;
    }
}

void Absolute_Value (double& num) {
    if (num < 0) {
        num *= -1;
    }
}

void Swap (int& num1, int& num2) {
    int temp = num1;
    num1 = num2;
    num2 = temp;
}

bool Are_Equal (double num, int m) {
    double diff = m - num;
    Absolute_Value (diff);
    if (diff < 0) {
        diff *= -1;
    }
    if (epsilon >= diff) {
        return true;
    }
    return false;
}

bool Are_Equal (double num1, double num2) {
    double diff = num1 - num2;
    if (diff < 0) {
        diff *= -1;
    }
    if (diff <= epsilon) {
        return true;
    }
    return false;
}

bool Is_An_Int (double num) {
    int test = num/1;
    if ( Are_Equal (num, test)) {
        return true;
    }
    return false;
}

unsigned int Greatest_Common_Deviser (int num1, int num2) {
    if (num1 == 0 || num2 == 0) {
        return 1;
    }
    Absolute_Value (num1);
    Absolute_Value (num2);
    if (num2 > num1) {
        Swap (num1, num2);
    }
    while (num2 > 0) {
        int temp = num1 % num2;
        num1 = num2;
        num2 = temp;
    }
    return num1;
}

int Get_Smallest (int a, int b) {
    if (a <= b) {
        return a;
    }
    return b;
}

int Get_Smallest (int a, int b, int c) {
    if (a <= Get_Smallest (b, c)) {
        return a;
    }
    else if (b <= Get_Smallest (a, c)) {
        return b;
    }
    return c;
}

void Put_In_Front_Of_Root (double& num_in_root, double& num_in_front_of_root) {
    if (!Is_An_Int (num_in_root)) {
        return;
    }
    int copy_of_num_in_root = num_in_root;
    int copy_of_num_in_front_of_root = num_in_front_of_root;
    if (copy_of_num_in_front_of_root != -1) {
        copy_of_num_in_front_of_root = 1;
    }
    if (copy_of_num_in_root < 0) {
        return;
    }
    double value_of_root = sqrt (copy_of_num_in_root);
    for ( unsigned int i = 2 ; i <= value_of_root ; ) {
        if (copy_of_num_in_root % (i*i) ==  0) {
            copy_of_num_in_front_of_root *= i;
            copy_of_num_in_root /= (i*i);
        }
        else {
            ++i;
        }
    }
    num_in_front_of_root = copy_of_num_in_front_of_root;
    num_in_root = copy_of_num_in_root;
}

void Reduce_Coefficients (double& a, double& b, double& c) {
    if (!Is_An_Int (a) && !Is_An_Int (b) && !Is_An_Int (c)) {
        return;
    }
    int gcd_of_a_and_b = Greatest_Common_Deviser (a, b);
    int gcd_of_a_and_c = Greatest_Common_Deviser (a, c);
    int gcd_of_b_and_c = Greatest_Common_Deviser (b, c);
    
    int smallest_gcd = Get_Smallest (gcd_of_a_and_b, gcd_of_a_and_c, gcd_of_b_and_c);
    if (smallest_gcd != 1) {
        a /= smallest_gcd;
        b /= smallest_gcd;
        c /= smallest_gcd;
    }
}

void Change_Signs (double& a, double& b, double& c) {
    if (a > 0) {
        return;
    }
    else if (a < 0 || (Are_Equal (a, 0) && b < 0)) {
        a *= -1;
        b *= -1;
        c *= -1;
    }
}

void Fix_Coefficients (double& a, double& b, double& c) {
    Reduce_Coefficients ( a, b, c );
    Change_Signs ( a, b, c );
}

void Reduce_Fraction (double& num1, double& num2) {
    if (!Is_An_Int (num1) && !Is_An_Int (num2)) {
        return;
    }
    int gcd = Greatest_Common_Deviser (num1, num2);
    num1 /= gcd;
    num2 /= gcd;
}

void Change_Signs (double& numerator, double& denominator) {
    if (denominator >= 0) {
        return;
    }
    else if (denominator < 0) {
        numerator *= -1;
        denominator *= -1;
    }
}

void Fix_Fraction (double& numerator, double& denominator) {
    Change_Signs (numerator, denominator);
    Reduce_Fraction (numerator, denominator);
}

void Print_Simple_Fraction (double numerator, double denominator) {
    Fix_Fraction (numerator, denominator);
    cout << numerator;
    if (!Are_Equal (numerator, 0) && !Are_Equal (denominator, 1)) {
        cout << over << denominator;
    }
}

void Print_Simple_Fractions (double numerator1, double denominator1, double numerator2, double denominator2) {
    cout << beg;
    Print_Simple_Fraction (numerator1, denominator1);
    cout << mid << sp;
    Print_Simple_Fraction (numerator2, denominator2);
    cout << fin;
}

void Print_Equation_Of_Line (double a, double b, double c) {
    Fix_Coefficients (a, b, c);
    if (Are_Equal (a, 0) && Are_Equal (b, 0)) {
        cout << "Invalid coefficients! ";
    }
    else if (Are_Equal (a, 0)) {
        cout << "y = ";
        if (c > 0) {
            cout << mi;
        }
        else if (c < 0) {
            Absolute_Value (c);
        }
        Print_Simple_Fraction (c, b);
    }
    else if (Are_Equal (b, 0)) {
        cout << "x = ";
        if (c > 0) {
            cout << mi;
        }
        else if (c < 0) {
            Absolute_Value ( c );
        }
        Print_Simple_Fraction ( c, a );
    }
    else {
        if (!Are_Equal (a, 1)) {
            cout << a;
        }
        cout << "x" << sp;
        if (b > 0) {
            cout << pl << sp;
            if ( b != 1 ) {
                cout << b;
            }
            cout << "y" << sp;
        }
        else if (b < 0) {
            cout << mi << sp;
            Absolute_Value ( b );
            if ( !Are_Equal (b, 1) ) {
                cout << b;
            }
            cout << "y" << sp;
        }
        if (c > 0) {
            cout << pl << sp << c << sp;
        }
        else if (c < 0) {
            cout << mi << sp;
            Absolute_Value ( c );
            cout << c << sp;
        }
        cout << "= 0";
    }
}

void Print_Equation_Of_Parabola (double a, double b, double c) {
    cout << "y = ";
    if (!Are_Equal (a, 1)) {
        cout << a;
    }
    cout << "x^2" << sp;
    if (b > 0) {
        cout << pl << sp;
        if (!Are_Equal (b, 1)) {
            cout << b;
        }
        cout << "x" << sp;
    }
    else if (b < 0) {
        cout << mi << sp;
        Absolute_Value (b);
        if (!Are_Equal (b, 1)) {
            cout << b;
        }
        cout << "x" << sp;
    }
    if (c > 0) {
        cout << pl << sp << c << sp;
    }
    else if (c < 0) {
        cout << mi << sp;
        Absolute_Value (c);
        cout << c;
    }
}

void Fill_Coefficients_Of_Lines_Crossing_TWo_Points (double x1, double y1, double x2, double y2, double& a, double& b, double& c) {
    a = 0;
    b = 0;
    c = 0;
    a = y1 - y2;
    b = x2 - x1;
    c = x1*y2 - x2*y1;
}

void Print_Line_Crossing_Two_Points (double x1, double y1, double x2, double y2) {
    double a, b, c;
    Fill_Coefficients_Of_Lines_Crossing_TWo_Points (x1, y1, x2, y2, a, b, c);
    Print_Equation_Of_Line (a, b, c);
}

bool Are_Overlapping_Lines (double a1, double b1, double c1, double a2, double b2, double c2) {
    Fix_Coefficients (a1, b1, c1);
    Fix_Coefficients (a2, b2, c2);
    
    if (Are_Equal (a1, a2) && Are_Equal (b1, b2) && Are_Equal (c1, c2)) {
        return true;
    }
    return false;
}

bool Are_Parallel_Lines (double a1, double b1, double c1, double a2, double b2, double c2) {
    Fix_Coefficients (a1, b1, c1);
    Fix_Coefficients (a2, b2, c2);
    if (Are_Equal (a1, 0) && Are_Equal (a2, 0) || Are_Equal (b1, 0) && Are_Equal (b2, 0)) {
        return true;
    }
    double coef1 = a1 / a2;
    double coef2 = b1 / b2;
    if ( Are_Equal (coef1, coef2))
    {
        a2 *= coef2;
        b2 *= coef2;
        c2 *= coef2;
    }
    if (Are_Equal (a1, a2) && Are_Equal (b1, b2) && !Are_Equal (c1, c2)) {
        return true;
    }
    return false;
}

bool Are_Penperdicular (double a1, double b1, double c1, double a2, double b2, double c2) {
    Fix_Coefficients (a1, b1, c1);
    Fix_Coefficients (a2, b2, c2);
    if (Are_Overlapping_Lines (a1, b1, c1, a2, b2, c2) || Are_Parallel_Lines (a1, b1, c1, a2, b2, c2)) {
        return false;
    }
    if ( Are_Equal (a1*a2, -b1*b2)) {
        return true;
    }
    return false;
}

bool Are_On_The_Same_Line ( double x1, double y1, double x2, double y2, double x3, double y3 ) {
    Fix_Fraction (x1, y1);
    Fix_Fraction (x2, y2);
    Fix_Fraction (x3, y3);

    if (Are_Equal (x1, x2) && Are_Equal (y1, y2) || Are_Equal (x1, x3) && Are_Equal (y1, y3) || Are_Equal (x2, x3) && Are_Equal (y2, y3)) {
        return true;
    }
    return false;
}

bool A_Point_Is_On_A_Line ( double x0, double y0, double a, double b, double c ) {
    double check = a*x0 + b*y0 + c;
    if ( Are_Equal ( check, 0 ) ) {
        return true;
    }
    return false;
}

void Fill_Point_Of_Interception_Of_Two_Lines (double a1, double b1, double c1, double a2, double b2, double c2, double& i1, double& i2, double& j1, double& j2) {
    i1 = b1*c2 - b2*c1;
    i2 = a1*b2 - a2*b1;
    Fix_Fraction ( i1, i2 );

    j1 = a2*c1 - a1*c2;
    j2 = a1*b2 - a2*b1;
    Fix_Fraction ( j1, j2 );
}

bool Are_Intercepting_Lines (double a1, double b1, double c1, double a2, double b2, double c2) {
    if (Are_Overlapping_Lines (a1, b1, c1, a2, b2, c2)) {
        return false;
    }
    else if (Are_Parallel_Lines (a1, b1, c1, a2, b2, c2)) {
        return false;
    }
    return true;
}

void Print_Point_Of_Interception_Of_Two_Lines (double a1, double b1, double c1, double a2, double b2, double c2) {
    if (Are_Overlapping_Lines (a1, b1, c1, a2, b2, c2)) {
        cout << "The two lines are overlapping! " << endl;
        return;
    }
    else if (Are_Parallel_Lines (a1, b1, c1, a2, b2, c2)) {
        cout << "The two lines are parallel! " << endl;
        return;
    }
    double i1 = b1*c2 - b2*c1;
    double i2 = a1*b2 - a2*b1;
    Fix_Fraction (i1, i2);

    double j1 = a2*c1 - a1*c2;
    double j2 = a1*b2 - a2*b1;
    Fix_Fraction (j1, j2);

    Print_Simple_Fractions (i1, i2, j1, j2);
}

void Parallel_Line_Passsing_Throuh_A_Point (double x0, double y0, double a, double b, double c) {
    double a0, b0, c0;
    if (A_Point_Is_On_A_Line (x0, y0, a, b, c)) {
        a0 = a;
        b0 = b;
        c0 = c;
    }
    else {
        a0 = a;
        b0 = b;
        c0 = - (a*x0 + b*y0);
    }
    Print_Equation_Of_Line (a0, b0, c0);
}

void Perpendicular_Line_Passsing_Throuh_A_Point (double x0, double y0, double a, double b, double c) {
    double a0 = b;
    double b0 = -a;
    double c0 = a*y0 - b*x0;
    Print_Equation_Of_Line (a0, b0, c0);
}

void Height (double x1, double y1, double x2, double y2, double x3, double y3, int decider) {
    double a0, b0, c0;
    if (decider == 1) {
        a0 = x2 - x3;
        b0 = y2 - y3;
        c0 = x1*(x3 - x2) + y1*(y3 - y2);
    }
    else if (decider == 2) {
        a0 = x1 - x3;
        b0 = y1 - y3;
        c0 = x2*(x3 - x1) + y2*(y3 - y1);
    }
    else if (decider == 3) {
        a0 = x1 - x2;
        b0 = y1 - y2;
        c0 = x3*(x2 - x1) + y3*(y2 - y1);
    }
    Print_Equation_Of_Line (a0, b0, c0);
}

void Median (double x1, double y1, double x2, double y2, double x3, double y3, int decider) {
    double a0, b0, c0;
    if (decider == 1) {
        a0 = y2 + y3 - 2*y1;
        b0 = 2*x1 - x2 - x3;
        c0 = y1*(x2 + x3) - x1*(y2 + y3);
    }
    else if (decider == 2) {
        a0 = y1 + y3 - 2*y2;
        b0 = 2*x2 - x1 - x3;
        c0 = y2*(x1 + x3) - x2*(y1 + y3);
    }
    else if (decider == 3) {
        a0 = y1 + y2 - 2*y3;
        b0 = 2*x3 - x1 - x2;
        c0 = y3*(x1 + x2) - x3*(y1 + y2);
    }
    Print_Equation_Of_Line (a0, b0, c0);
}

void Bisection (double x1, double y1, double x2, double y2, double x3, double y3, int decider) {
    double a0, b0, c0;
    if (decider == 1) {
        a0 = 2*x2 - 2*x3;
        b0 = 2*y2 - 2*y3;
        c0 = x3*x3 - x2*x2 + y3*y3 - y2*y2;
    }
    else if (decider == 2) {
        a0 = 2*x1 - 2*x3;
        b0 = 2*y1 - 2*y3;
        c0 = x3*x3 - x1*x1 + y3*y3 - y1*y1;
    }
    else if (decider == 3) {
        a0 = 2*x1 - 2*x2;
        b0 = 2*y1 - 2*y2;
        c0 = x2*x2 - x1*x1 + y2*y2 - y1*y1;
    }
    Print_Equation_Of_Line (a0, b0, c0);
}

void Fill_Point_Of_Interception_Of_A_Line_And_A_Parabola (double a1, double b1, double c1, double a2, double b2, double c2, double& i1, double& i2, double& i3, double& i4, double& i5, double& j1, double& j2, double& j3, double& j4, double& j5, int point) {
    Fix_Coefficients (a1, b1, c1);
    Fix_Coefficients (a2, b2, c2);
    i1 = 0; j1 = 0;
    i2 = 0; j2 = 0;
    i3 = 1; j3 = 1;
    i4 = 0; j4 = 0;
    i5 = 0; j5 = 0;
    
    double discriminant = (a2 + b1*b2)*(a2 + b1*b2) - 4*a1*b2*(b2*c1 + c2);
    if (discriminant < 0) {
        return;
    }
    else if (Are_Equal (discriminant, 0)) {
        i1 = - (a2 + b1*b2);
        i2 = 2*a1*b2;
        Fix_Fraction (i1, i2);
    
        j1 = a2*a2 + a2*b1*b2 - 2*a1*b2*c2;
        j2 = 2*a1*b2*b2;
        Fix_Fraction (j1, j2);
        
        i4 = -1;
        j4 = -1;
    }
    else if (discriminant > 0 && point == 1) {
        i1 = - (a2 + b1*b2);
        i2 = 2*a1*b2;
        i4 = discriminant;
        Put_In_Front_Of_Root (i4, i3);
        i5 = i2;
        if (Are_Equal (i4, 1)) {
            i1 += i3;
            i4 = -1;
            i3 = 0;
        }
        Fix_Fraction (i1, i2);
        Fix_Fraction (i3, i5);
        
        j1 = a2*a2 + a2*b1*b2 - 2*a1*b2*c2;
        j2 = 2*a1*b2*b2;
        j4 = discriminant;
        Put_In_Front_Of_Root (j4, j3);
        j3 *= -a2;
        j5 = j2;
        if (Are_Equal (j4, 1)) {
            j1 += j3;
            j4 = -1;
            j3 = 0;
        }
        Fix_Fraction (j1, j2);
        Fix_Fraction (j3, j5);
    }
    else if (discriminant > 0 && point == 2) {
        i1 = - (a2 + b1*b2);
        i2 = 2*a1*b2;
        i4 = discriminant;
        Put_In_Front_Of_Root (i4, i3);
        i3 *= -1;
        i5 = i2;
        if (i4 == 1) {
            i1 += i3;
            i4 = -1;
            i3 = 0;
        }
        Fix_Fraction (i1, i2);
        Fix_Fraction (i3, i5);
        
        j1 = a2*a2 + a2*b1*b2 - 2*a1*b2*c2;
        j2 = 2*a1*b2*b2;
        j4 = discriminant;
        Put_In_Front_Of_Root (j4, j3);
        j3 *= a2;
        j5 = j2;
        if (j4 == 1) {
            j1 += j3;
            j4 = -1;
            j3 = 0;
        }
        Fix_Fraction (j1, j2);
        Fix_Fraction (j3, j5);
    }
}

void Print_Fraction_With_Root (double i1, double i2, double i3, double i4, double i5, double j1, double j2, double j3, double j4, double j5) {
    if (Are_Equal (i4, -1) && Are_Equal (j4, -1)) {
        Print_Simple_Fractions (i1, i2, j1, j2);
        return;
    }
    cout << beg;
    if (!Are_Equal (i1, 0)) {
        cout << i1;
        if (!Are_Equal (i2, 1)) {
            cout << over << i2 << sp;
        }
        else if (Are_Equal (i2, 1)) {
            cout << sp;
        }
    }
    if (i3 < 0) {
        cout << mi;
        Absolute_Value (i3);
        if (!Are_Equal (i1, 0)) {
            cout << sp;
        }
    }
    else if (i3 > 0 && !Are_Equal (i1, 0)) {
        cout << pl << sp;
    }
    if (!Are_Equal (i3, 1)) {
        cout << i3 << mult;
    }
    cout << sq << root_beg << i4 << root_end;
    if (!Are_Equal (i5, 1)) {
        cout << over << i5;
    }
    cout << mid << sp;
    if (!Are_Equal (j1, 0)) {
        cout << j1;
        if (!Are_Equal (j2, 1)) {
            cout << over << j2 << sp;
        }
        else if (Are_Equal (j2, 1)) {
            cout << sp;
        }
    }
    if (j3 < 0) {
        cout << mi;
        Absolute_Value (j3);
        if (!Are_Equal (j1, 0)) {
            cout << sp;
        }
    }
    else if ( j3 > 0 && !Are_Equal (j1, 0)) {
        cout << pl << sp;
    }
    if (!Are_Equal (j3, 1)) {
        cout << j3 << mult;
    }
    cout << sq << root_beg << j4 << root_end;
    if (!Are_Equal (j5, 1)) {
        cout << over << j5;
    }
    cout << fin;
}

void Print_Point_Of_Interception_Of_A_Line_And_A_Parabola (double a1, double b1, double c1, double a2, double b2, double c2) {
    double i1, i2, i3, i4, i5, j1, j2, j3, j4, j5;
    int point = 0;

    if ( Are_Equal ( a1, 0 ) && Are_Equal ( b1, 0 ) ) {
        cout << "Invalid coefficients for the parabola!";
        return;
    }
    if ( Are_Equal ( a2, 0 ) && Are_Equal ( b2, 0 ) ) {
        cout << "Invalid coefficients for the line!";
        return;
    }
    if ( Are_Equal ( b2, 0 ) ) {
        double numerator1 = -c2;
        double denominator1 = a2;
        double numerator2 = a1*c2*c2 - a2*b1*c2 + a2*a2*c1;
        double denominator2 = a2*a2; 
        Print_Simple_Fractions ( numerator1, denominator1, numerator2, denominator2 );
        return;
    }
    
    double discriminant = ( a2 + b1*b2 )*( a2 + b1*b2 ) - 4*a1*b2*( b2*c1 + c2 );
    if ( discriminant < 0 ) {
        cout << "No real solutions! " << endl;
        return;
    }
    else if ( Are_Equal ( discriminant, 0 ) ) {
        Fill_Point_Of_Interception_Of_A_Line_And_A_Parabola ( a1, b1, c1, a2, b2, c2, i1, i2, i3, i4, i5, j1, j2, j3, j4, j5, point );
        Print_Fraction_With_Root ( i1, i2, i3, i4, i5, j1, j2, j3, j4, j5 );
        return;
    }
    point = 1;
    cout << endl;
    Fill_Point_Of_Interception_Of_A_Line_And_A_Parabola ( a1, b1, c1, a2, b2, c2, i1, i2, i3, i4, i5, j1, j2, j3, j4, j5, point );
    Print_Fraction_With_Root ( i1, i2, i3, i4, i5, j1, j2, j3, j4, j5 );
    point = 2;
    Fill_Point_Of_Interception_Of_A_Line_And_A_Parabola ( a1, b1, c1, a2, b2, c2, i1, i2, i3, i4, i5, j1, j2, j3, j4, j5, point );
    cout << endl;
    Print_Fraction_With_Root ( i1, i2, i3, i4, i5, j1, j2, j3, j4, j5 );
}

void Fill_Coefficients_Of_Lines_Crossing_A_Point_Tanget_To_A_Parabola ( double x0, double y0, double a, double b, double c, double& t1, double& t2, double& t3, double& t4, double& t5, double& t6, double& t7, double& t8, int decider ) {
    t1 = 0;
    t2 = 1;
    t3 = 0;
    t4 = -1;
    t5 = 0;
    t6 = 0;
    t7 = -1;
    t8 = 0;
    
    double reduced_discriminant = a*a*x0*x0 - a*y0 + a*b*x0 + a*c;
    t3 = reduced_discriminant;
    t5 = reduced_discriminant;
    t8 = reduced_discriminant;
    t1 = a*a*x0*x0 + a*b*x0 + a*c + reduced_discriminant - a*y0;
    t6 = a*x0*y0 - a*a*x0*x0*x0 - a*b*x0*x0 - a*c*x0 - x0*reduced_discriminant;
    
    if ( t1 < 0 || ( Are_Equal ( t1, 0 ) && decider == 2 ) ) {
        t1 *= -1;
        t2 *= -1;
        t4 *= -1;
        t6 *= -1;
        t7 *= -1;
    }
    Put_In_Front_Of_Root ( t3, t2 );
    Put_In_Front_Of_Root ( t5, t4 );
    Put_In_Front_Of_Root ( t8, t7 );
    t2 *= 2*a*x0 + b;
    t7 *= (2*a*x0*x0 + b*x0 - y0) ;
    
    if ( decider == 2 ) {
        t2 *= -1;
        t4 *= -1;
        t7 *= -1;
    }
}

void Print_Equation_Of_Line_With_Roots (double t1, double t2, double t3, double t4, double t5, double t6, double t7, double t8) {
    if (Are_Equal (t3, 0) && Are_Equal (t5, 0) && Are_Equal (t8, 0)) {
        Print_Equation_Of_Line (t1, t4, t5);
        return;
    }
    if (Are_Equal (t3, 1) && Are_Equal (t5, 1) && Are_Equal (t8, 1)) {
        double new_a = t1 + t2;
        double new_b = t4;
        double new_c = t6 + t7;
        Print_Equation_Of_Line (new_a, new_b, new_c);
        return;
    }
    cout << beg;
    if (!Are_Equal (t1, 0)) {
        cout << t1 << sp;
    }
    if (t2 > 0 && !Are_Equal (t1, 0)) {
        cout << pl << sp;
    }
    else if (t2 < 0 && !Are_Equal (t1, 0)) {
        cout << mi << sp;
    }
    Absolute_Value ( t2 );
    if (!Are_Equal (t2, 1)) {
        cout << t2 << mult;
    }
    cout << sq << root_beg << t3 << root_end << fin << "x" << sp;
    if (t4 > 0) {
        cout << pl;
    }
    else if (t4 < 0) {
        cout << mi;
    }
    Absolute_Value ( t4 );
    if (!Are_Equal (t4, 1)) {
        cout << t4 << mult;
    }
    cout << sq << root_beg << t5 << root_end << "y" << sp;
    if ( t6 > 0 ) {
        cout << pl << sp << t6 << sp;
    }
    else if (t6 < 0) {
        cout << mi << sp;
        Absolute_Value ( t6 );
        cout << t6 << sp;
    }
    if (t7 > 0) {
        cout << pl << sp;
    }
    else if (t7 < 0) {
        cout << mi << sp;
    }
    Absolute_Value (t7);
    if (!Are_Equal (t7, 1)) {
        cout << t7 << mult;
    }
    cout << sq << root_beg << t8 << root_end << sp << "= 0";
}

void Print_Lines_Crossing_A_Point_Tangent_To_A_Parabola (double x0, double y0, double a, double b, double c) {
    double t1, t2, t3, t4, t5, t6, t7, t8;
    Fix_Coefficients (a, b, c);

    double reduced_discriminant = a*a*x0*x0 - a*y0 + a*b*x0 + a*c;
    if (Are_Equal (a, 0)) {
        cout << "Invalid equation for a parabola!" << endl;
        return;
    }
    if (Are_Equal (a*x0*x0 + b*x0 + c, y0)) {
        cout << endl;
        Print_Equation_Of_Line (2*a*x0 + b, -1, c - a*x0*x0);
        return;
    }

    int decider = 1;
    Fill_Coefficients_Of_Lines_Crossing_A_Point_Tanget_To_A_Parabola (x0, y0, a, b, c, t1, t2, t3, t4, t5, t6, t7, t8, decider);
    cout << endl;
    Print_Equation_Of_Line_With_Roots (t1, t2, t3, t4, t5, t6, t7, t8);
    if (Are_Equal (reduced_discriminant, 0)) {
        return;
    }
    decider = 2;
    Fill_Coefficients_Of_Lines_Crossing_A_Point_Tanget_To_A_Parabola (x0, y0, a, b, c, t1, t2, t3, t4, t5, t6, t7, t8, decider);
    cout << endl;
    Print_Equation_Of_Line_With_Roots (t1, t2, t3, t4, t5, t6, t7, t8);
}

bool The_Two_Points_Lie_On_The_Same_Side_Of_A_Line (double a, double b, double c, double x0, double y0, double x1, double y1) {
    if (A_Point_Is_On_A_Line (x0, y0, a, b, c) || A_Point_Is_On_A_Line (x1, y1, a, b, c)) {
        return false;
    }

    double height_of_line_relative_to_x0 = ( -c - a*x0 )*b;
    double height_of_line_relative_to_x1 = ( -c - a*x1 )*b;
    
    if (y0 > height_of_line_relative_to_x0 && y1 > height_of_line_relative_to_x1) {
        return true;
    }
    if (y0 < height_of_line_relative_to_x0 && y1 < height_of_line_relative_to_x1) {
        return true;
    }
    return false;
}

bool A_Point_Is_In_A_Triangle (double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
    double a1, b1, c1, a2, b2, c2, a3, b3, c3;
    Fill_Coefficients_Of_Lines_Crossing_TWo_Points (x1, y1, x2, y2, a1, b1, c1);
    Fill_Coefficients_Of_Lines_Crossing_TWo_Points (x1, y1, x3, y3, a2, b2, c2);
    Fill_Coefficients_Of_Lines_Crossing_TWo_Points (x2, y2, x3, y3, a3, b3, c3);

    if (!The_Two_Points_Lie_On_The_Same_Side_Of_A_Line (a1, b1, c1, x3, y3, x0, y0)) {
        return false;
    }
    if (!The_Two_Points_Lie_On_The_Same_Side_Of_A_Line ( a2, b2, c2, x2, y2, x0, y0)) {
        return false;
    }
    if (!The_Two_Points_Lie_On_The_Same_Side_Of_A_Line (a3, b3, c3, x1, y1, x0, y0)) {
        return false;
    }
    return true;
}

bool Is_A_Concave_Quadrilateral ( double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4 ) {
    if ( A_Point_Is_In_A_Triangle ( x4, y4, x1, y1, x2, y2, x3, y3 ) ) {
        return true;
    }
    if ( A_Point_Is_In_A_Triangle ( x3, y3, x1, y1, x2, y2, x4, y4 ) ) {
        return true;
    }
    if ( A_Point_Is_In_A_Triangle ( x2, y2, x1, y1, x3, y3, x4, y4 ) ) {
        return true;
    }
    if ( A_Point_Is_In_A_Triangle ( x1, y1, x2, y2, x3, y3, x4, y4 ) ) {
        return true;
    }
    return false;
}

double Get_Length (double x1, double y1, double x2, double y2) {
    double distance = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
    return distance;
}

bool Is_A_Parallelogram (double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3, double a4, double b4, double c4) {
    if (Are_Parallel_Lines (a1, b1, c1, a3, b3, c3) && Are_Parallel_Lines (a2, b2, c2, a4, b4, c4)) {
        return true;
    }
    if (Are_Parallel_Lines (a1, b1, c1, a2, b2, c2) && Are_Parallel_Lines (a3, b3, c3, a4, b4, c4)) {
        return true;
    }
    if (Are_Parallel_Lines (a1, b1, c1, a4, b4, c4) && Are_Parallel_Lines (a3, b3, c3, a2, b2, c2)) {
        return true;
    }
    return false;
}

bool Is_A_Rectangle (double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3, double a4, double b4, double c4) {
    if (!Is_A_Parallelogram (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        return false;
    }
    if (Are_Penperdicular (a1, b1, c1, a2, b2, c2) || Are_Penperdicular (a1, b1, c1, a3, b3, c3) || Are_Penperdicular (a1, b1, c1, a4, b4, c4)) {
        return true;
    }
    return false;
}

bool Has_Equal_Sides (double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3, double a4, double b4, double c4) {
    double x1, y1, x2, y2, x3, y3;
    if (Are_Parallel_Lines (a1, b1, c1, a2, b2, c2)) {
        x1 = (b1*c3 - b3*c1)/(a1*b3 - a3*b1);
        y1 = (a3*c1 - a1*c3)/(a1*b3 - a3*b1);

        x2 = (b1*c4 - b4*c1)/(a1*b4 - a4*b1);
        y2 = (a4*c1 - a1*c4)/(a1*b4 - a4*b1);

        x3 = (b2*c4 - b4*c2)/(a2*b4 - a4*b2);
        y3 = (a4*c2 - a2*c4)/(a2*b4 - a4*b2);
    }
    else if (Are_Parallel_Lines (a1, b1, c1, a3, b3, c3)) {
        x1 = (b1*c2 - b2*c1)/(a1*b2 - a2*b1);
        y1 = (a2*c1 - a1*c2)/(a1*b2 - a2*b1);

        x2 = (b1*c4 - b4*c1)/(a1*b4 - a4*b1);
        y2 = (a4*c1 - a1*c4)/(a1*b4 - a4*b1);

        x3 = (b2*c3 - b3*c2)/(a2*b3 - a3*b2);
        y3 = (a3*c2 - a2*c3)/(a2*b3 - a3*b2);
    }
    else if (Are_Parallel_Lines (a1, b1, c1, a4, b4, c4)) {
        x1 = (b1*c2 - b2*c1)/(a1*b2 - a2*b1);
        y1 = (a2*c1 - a1*c2)/(a1*b2 - a2*b1);

        x2 = (b1*c3 - b3*c1)/(a1*b3 - a3*b1);
        y2 = (a3*c1 - a1*c3)/(a1*b3 - a3*b1);

        x3 = (b2*c4 - b4*c2)/(a2*b4 - a4*b2);
        y3 = (a4*c2 - a2*c4)/(a2*b4 - a4*b2);
    }
    double len12 = Get_Length (x1, y1, x2, y2);
    double len13 = Get_Length (x1, y1, x3, y3);
    double len23 = Get_Length (x2, y2, x3, y3);
    if (Are_Equal (len12, len23) || Are_Equal (len12, len13)) {
        return true;
    }
    return false;
}

bool Is_A_Square (double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3, double a4, double b4, double c4) {
    if (!Is_A_Rectangle (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        return false;
    }
    if (Has_Equal_Sides (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        return true;
    }
    return false;
}

bool Is_A_Rhombus (double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3, double a4, double b4, double c4) {
    if (!Is_A_Parallelogram (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        return false;
    }
    if (Has_Equal_Sides (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        return true;
    }
    return false;
}

bool Is_A_Trapezoid (double a1, double b1, double c1, double a2, double b2, double c2, double a3, double b3, double c3, double a4, double b4, double c4) {
    if (Are_Parallel_Lines (a1, b1, c1, a3, b3, c3) && Are_Intercepting_Lines (a2, b2, c2, a4, b4, c4)) {
        return true;
    }
    if (Are_Parallel_Lines (a1, b1, c1, a2, b2, c2) && Are_Intercepting_Lines (a3, b3, c3, a4, b4, c4)) {
        return true;
    }
    if (Are_Parallel_Lines (a1, b1, c1, a4, b4, c4) && Are_Intercepting_Lines (a3, b3, c3, a2, b2, c2)) {
        return true;
    }
    return false;
}

bool Is_A_Valid_Symbol (char symbol) {
    if (symbol >= '0' && symbol <= '9' || symbol >= 'a' && symbol <= 'z' || symbol >= 'A' && symbol <= 'Z' || symbol == '_') {
        return true;
    }
    return false;
}

void Give_Name (char* name) {
    cin.ignore ();
    cin.getline(name, max_size);
    bool is_wrong = false;
    for (unsigned int i = 0 ; name[i] != '\0' ; ++i) {
        if (!Is_A_Valid_Symbol (name[i])) {
            is_wrong = true;
        }
    }
    if (is_wrong == true) {
        cout << "Invalid name! Enter a new name!" << endl;
        cin.ignore();
        Give_Name (name);
    }
}

void Enter_A_Number (double& n) {
    while (true) {
        cin >> n;
        if (cin.good()) {
            break; // n is a valid number
        }
        else {
            cout << "Invalid number! Enter a new number!" << endl;
            cin.clear();
            cin.ignore(100000, '\n');
        }
    }
}

void Enter_A_Parabola (char* parabola, double& a, double& b, double& c) {
    cout << "Enter a name for the parabola: ";
    Give_Name (parabola);
    cout << "Enter the coordinates of the parabola {a b c}: ";
    Enter_A_Number (a);
    Enter_A_Number (b);
    Enter_A_Number (c);
}

void Enter_A_Point (char* point, double& x, double& y) {
    cout << "Enter a name for the point: ";
    Give_Name (point);
    cout << "Enter the coordinates of the point {x y}: ";
    Enter_A_Number (x);
    Enter_A_Number (y);
}

void Enter_Two_Points (char* point1, char* point2, double& x1, double& y1, double& x2, double& y2) {
    cout << "Enter a name for the first point: ";
    Give_Name (point1);
    cout << "Enter the coordinates of the point {x1 y1}: ";
    Enter_A_Number (x1);
    Enter_A_Number (y1);

    cout << "Enter a name for the second point: ";
    Give_Name (point2);
    cout << "Enter the coordinates of the point {x2 y2}: ";
    Enter_A_Number (x2);
    Enter_A_Number (y2);
}

void Enter_Three_Points (char* point1, char* point2, char* point3, double& x1, double& y1, double& x2, double& y2, double& x3, double& y3) {
    cout << "Enter a name for the first point: ";
    Give_Name (point1);
    cout << "Enter the coordinates of the point {x1 y1}: ";
    Enter_A_Number (x1);
    Enter_A_Number (y1);

    cout << "Enter a name for the second point: ";
    Give_Name (point2);
    cout << "Enter the coordinates of the point {x2 y2}: ";
    Enter_A_Number (x2);
    Enter_A_Number (y2);

    cout << "Enter a name for the third point: ";
    Give_Name (point3);
    cout << "Enter the coordinates of the point {x3 y3}: ";
    Enter_A_Number (x3);
    Enter_A_Number (y3);
}

void Enter_A_Line (char* line, double& a, double& b, double& c) {
    cout << "Enter a name for the line: ";
    Give_Name (line);
    cout << "Enter the coordinates of the line {a b c}: ";
    Enter_A_Number (a);
    Enter_A_Number (b);
    Enter_A_Number (c);
}

void Enter_Two_Lines (char* line1, char* line2, double& a1, double& b1, double& c1, double& a2, double& b2, double& c2) {
    cout << "Enter a name for the first line: ";
    Give_Name (line1);
    cout << "Enter the coordinates of the line {a1 b1 c1}: ";
    Enter_A_Number (a1);
    Enter_A_Number (b1);
    Enter_A_Number (c1);

    cout << "Enter a name for the second line: ";
    Give_Name (line2);
    cout << "Enter the coordinates of the line {a2 b2 c2}: ";
    Enter_A_Number (a2);
    Enter_A_Number (b2);
    Enter_A_Number (c2);
}

void Enter_Four_Lines (char* line1, char* line2, char* line3, char* line4, double& a1, double& b1, double& c1, double& a2, double& b2, double& c2, double& a3, double& b3, double& c3, double& a4, double& b4, double& c4) {
    cout << "Enter a name for the first line: ";
    Give_Name (line1);
    cout << "Enter the coordinates of the line {a1 b1 c1}: ";
    Enter_A_Number (a1);
    Enter_A_Number (b1);
    Enter_A_Number (c1);

    cout << "Enter a name for the second line: ";
    Give_Name (line2);
    cout << "Enter the coordinates of the line {a2 b2 c2}: ";
    Enter_A_Number (a2);
    Enter_A_Number (b2);
    Enter_A_Number (c2);

    cout << "Enter a name for the third line: ";
    Give_Name (line3);
    cout << "Enter the coordinates of the line {a3 b3 c3}: ";
    Enter_A_Number (a3);
    Enter_A_Number (b3);
    Enter_A_Number (c3);

    cout << "Enter a name for the fourth line: ";
    Give_Name (line4);
    cout << "Enter the coordinates of the line {a4 b4 c4}: ";
    Enter_A_Number (a4);
    Enter_A_Number (b4);
    Enter_A_Number (c4);
}

void Determine_Vertex (const char* point1, const char* point2, const char* point3, char* name_of_vertex, int& decider, double x1, double y1, double x2, double y2, double x3, double y3, double& x0, double& y0) {
    Give_Name (name_of_vertex);
    if (strcmp (point1, name_of_vertex) == 0) {
        decider = 1;
        x0 = x1;
        y0 = y1;
        return;
    }
    else if (strcmp (point2, name_of_vertex) == 0) {
        decider = 2;
        x0 = x2;
        y0 = y2;
        return;
    }
    else if (strcmp (point3, name_of_vertex) == 0) {
        decider = 3;
        x0 = x3;
        y0 = y3;
        return;
    }
    else {
        cout << "Not a name of a vertex!" << endl;
        Determine_Vertex (point1, point2, point3, name_of_vertex, decider, x1, y1, x2, y2, x3, y3, x0, y0);
    }
}

void Initialisation () {
    cout << "1. Find the equation of a line from 2 given points. " << endl;
    cout << "2. Check if a given point lies on a given line. " << endl;
    cout << "3. Find the equation of a line, that intersects a given point and is parallel to another given line " << endl;
    cout << "4. Find the equation of a line, that intersects a given point and is perpendicular to another given line " << endl;
    cout << "5. Find the coordinates of the intersection of two given lines " << endl;
    cout << "6. Find the equation of a height of a triangle by its vertices given as coordinates of points " << endl;
    cout << "7. Find the equation of a median of a triangle by its vertices given as coordinates of points " << endl;
    cout << "8. Find the equation of a bisection of a triangle by its vertices given as coordinates of points " << endl;
    cout << "9. Find the equations of lines that are tangent to a parabola, given by its coefficients, and intersects a given point " << endl;
    cout << "10. Find the coordinates of the points of intersection of a parabola and a line, given by their coeficinets " << endl;
    cout << "11. Determine the type of a quadrilateral, made from the intersection of 4 given lines " << endl;
    cout << "12. Quit" << endl;
}

void Initial_Choice (double& N) {
    cout << endl << "Enter the number of one of the options above: ";
    Enter_A_Number (N);
    if (N < 1 || N > number_of_options) {
        cout << "Invalid number! Enter a new number!" << endl;
        Initial_Choice (N);
    }
}

void Option_1 () {
    char point1[max_size];
    char point2[max_size];
    double x1, y1, x2, y2;
    Enter_A_Point (point1, x1, y1);
    Enter_A_Point (point2, x2, y2);

    if (Are_Equal (x1, x2) && Are_Equal (y1, y2)) {
        cout << endl << point1 << " and " << point2 << " must be diferent!" << endl;
        return;
    }
    cout << endl << "The equation of the line, crossing points " << point1 << " " << beg << x1 << ", " << y1 << fin ;
    cout << " and " << point2 << " " << beg << x2 << ", " << y2 << fin << ", is: " << endl;
    Print_Line_Crossing_Two_Points (x1, y1, x2, y2);
    cout << endl;
}

void Option_2 () {
    char point[max_size];
    char line[max_size];
    double x, y, a, b, c;
    Enter_A_Point (point, x, y);
    Enter_A_Line (line, a, b, c);

    cout << point << " " << beg << x << ", " << y << fin;
    if (A_Point_Is_On_A_Line (x, y, a, b, c)) {
        cout << " lies on line ";
    }
    else {
        cout << " does not lie on line ";
    }
    cout << line << ": ";
    Print_Equation_Of_Line (a, b, c);
    cout << endl;
}

void Option_3 () {
    char point[max_size];
    char line[max_size];
    double x, y, a, b, c;
    Enter_A_Point (point, x, y);
    Enter_A_Line (line, a, b, c);

    cout << endl << "The equation of the line, crossing point " << point << " " << beg << x << ", " << y << fin ;
    cout << " and is parallel to line " << line << ": ";
    Print_Equation_Of_Line (a, b, c);
    cout << ", is: " << endl;
    Parallel_Line_Passsing_Throuh_A_Point (x, y, a, b, c);
    cout << endl;
}

void Option_4 () {
    char point[max_size];
    char line[max_size];
    double x, y, a, b, c;
    Enter_A_Point (point, x, y);
    Enter_A_Line (line, a, b, c);

    cout << endl << "The equation of the line, crossing point " << point << " " << beg << x << ", " << y << fin ;
    cout << " and is perpendicular to line " << line << ": ";
    Print_Equation_Of_Line (a, b, c);
    cout << ", is: " << endl;
    Perpendicular_Line_Passsing_Throuh_A_Point (x, y, a, b, c);
    cout << endl;
}

void Option_5 () {
    char line1[max_size];
    char line2[max_size];
    double a1, b1, c1, a2, b2, c2;
   
    Enter_Two_Lines (line1, line2, a1, b1, c1, a2, b2, c2);

    if ( Are_Overlapping_Lines ( a1, b1, c1, a2, b2, c2 ) ) {
        cout << "The two lines are overlapping! " << endl;
        return;
    }
    if ( Are_Parallel_Lines ( a1, b1, c1, a2, b2, c2 ) ) {
        cout << "The two lines are parallel! " << endl;
        return;
    }

    cout << "The coordinates of the point of interception of lines " << line1 << ": ";
    Print_Equation_Of_Line (a1, b1, c1);
    cout << " and " << line2 << ": ";
    Print_Equation_Of_Line (a2, b2, c2);
    cout << " are " << endl;
    Print_Point_Of_Interception_Of_Two_Lines (a1, b1, c1, a2, b2, c2);
    cout << endl;
}

void Option_6 () {
    char point1[max_size];
    char point2[max_size];
    char point3[max_size];
    double x1, y1, x2, y2, x3, y3;
    Enter_Three_Points (point1, point2, point3, x1, y1, x2, y2, x3, y3);
    if (Are_On_The_Same_Line (x1, y1, x2, y2, x3, y3)) {
        cout << "The points do not form a triangle!" << endl;
        return;
    }
    cout << "Enter the name of the vertex which the heigth crosses: ";
    char str[max_size];
    int decider;
    double x0, y0;
    Determine_Vertex (point1, point2, point3, str, decider, x1, y1, x2, y2, x3, y3, x0, y0);

    cout << "The equation of the height from the vertex " << str << " " << beg << x0 << ", " <<  y0 << fin << " is:" << endl;
    Height (x1, y1, x2, y2, x3, y3, decider);
    cout << endl;
}

void Option_7 () {
    char point1[max_size];
    char point2[max_size];
    char point3[max_size];
    double x1, y1, x2, y2, x3, y3;
    Enter_Three_Points (point1, point2, point3, x1, y1, x2, y2, x3, y3);
    if (Are_On_The_Same_Line (x1, y1, x2, y2, x3, y3)) {
        cout << "The points do not form a triangle!" << endl;
        return;
    }
    cout << "Enter the name of the vertex which the median crosses: ";
    char str[max_size];
    int decider;
    double x0, y0;
    Determine_Vertex (point1, point2, point3, str, decider, x1, y1, x2, y2, x3, y3, x0, y0);

    cout << "The equation of the median from the vertex " << str << " " << beg << x0 << ", " <<  y0 << fin << " is:" << endl;
    Median (x1, y1, x2, y2, x3, y3, decider);
    cout << endl;
}

void Option_8 () {
    char point1[max_size];
    char point2[max_size];
    char point3[max_size];
    double x1, y1, x2, y2, x3, y3;
    Enter_Three_Points (point1, point2, point3, x1, y1, x2, y2, x3, y3);
    if (Are_On_The_Same_Line (x1, y1, x2, y2, x3, y3)) {
        cout << "The points do not form a triangle!" << endl;
        return;
    }
    cout << "Enter the name of the vertex which the heigth crosses: ";
    char str[max_size];
    int decider;
    double x0, y0;
    Determine_Vertex (point1, point2, point3, str, decider, x1, y1, x2, y2, x3, y3, x0, y0);

    cout << "The equation of the bisection opposite the vertex " << str << " " << beg << x0 << ", " <<  y0 << fin << " is:" << endl;
    Bisection (x1, y1, x2, y2, x3, y3, decider);
    cout << endl;
}

void Option_9 () {
    char point[max_size];
    char parabola[max_size];
    double x, y, a, b, c;
    Enter_A_Point (point, x, y);
    Enter_A_Parabola (parabola, a, b, c);
    if (Are_Equal (a, 0)) {
        cout << "The leading coeficient of the parabola cannot be zero!" << endl;
        return;
    }
    double reduced_discriminant = a*a*x*x - a*y + a*b*x + a*c;
    if ( reduced_discriminant < 0 ) {
        cout << "Such lines do not exist!" << endl;
        return;
    }
    else if (Are_Equal (reduced_discriminant, 0)) {
        cout << "The equation of the line, tangent to " << parabola << ": ";
        Print_Equation_Of_Parabola (a, b, c);
        cout << " and crossing " << point << " " << beg << x << ", " << y << fin << ", is: ";
        Print_Lines_Crossing_A_Point_Tangent_To_A_Parabola (x, y, a, b, c);
        cout << endl;
        return;
    }
    cout << "The equations of the lines, tangent to " << parabola << ": ";
    Print_Equation_Of_Parabola (a, b, c);
    cout << " and crossing " << point << " " << beg << x << ", " << y << fin << ", are: ";
    Print_Lines_Crossing_A_Point_Tangent_To_A_Parabola (x, y, a, b, c);
    cout << endl;
}

void Option_10 () {
    char parabola[max_size];
    char line[max_size];;
    double a1, b1, c1, a2, b2, c2;
    Enter_A_Parabola (parabola, a1, b1, c1);
    Enter_A_Line (line, a2, b2, c2);
    double discriminant = (a2 + b1*b2)*(a2 + b1*b2) - 4*a1*b2*(b2*c1 + c2);
    if (Are_Equal (a1, 0)) {
        cout << "The leading coeficient of the parabola cannot be zero!" << endl;
        return;
    }
    else if ( discriminant < 0 ) {
        cout << "Such points do not exist!" << endl;
        return;
    }

    cout << "The coordinates of the points of interception of " << parabola << ": ";
    Print_Equation_Of_Parabola (a1, b1, c1);
    cout << " and " << line << ": ";
    Print_Equation_Of_Line (a2, b2, c2);
    cout << " are ";
    Print_Point_Of_Interception_Of_A_Line_And_A_Parabola (a1, b1, c1, a2, b2, c2);
    cout << endl;
}

void Option_11 () {
    char line1[max_size];
    char line2[max_size];
    char line3[max_size];
    char line4[max_size];
    double a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4;
    Enter_Four_Lines (line1, line2, line3, line4, a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4);

    if (Is_A_Square (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        cout << "The quadrilateral is a square" << endl;
        return;
    }
    if (Is_A_Rhombus (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        cout << "The quadrilateral is a rhombus" << endl;
        return;
    }
    if (Is_A_Rectangle (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        cout << "The quadrilateral is a rectangle" << endl;
        return;
    }
    if (Is_A_Parallelogram (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        cout << "The quadrilateral is a parallelogram" << endl;
        return;
    }
    if (Is_A_Trapezoid (a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4)) {
        cout << "The quadrilateral is a trapezoid" << endl;
        return;
    }
    else {
        cout << "The quadrilateral does not have a special shape" << endl;
    }
}

void Option_12 (bool& want_to_quit) {
    want_to_quit = false;
}

void Choice (int N, bool& want_to_quit) {
    if (N == 1) {
        Option_1 ();
    }
    else if (N == 2) {
        Option_2 ();
    }
    else if (N == 3) {
        Option_3 ();
    }
    else if (N == 4) {
        Option_4 ();
    }
    else if (N == 5) {
        Option_5 ();
    }
    else if (N == 6) {
        Option_6 ();
    }
    else if (N == 7) {
        Option_7 ();
    }
    else if (N == 8) {
        Option_8 ();
    }
    else if (N == 9) {
        Option_9 ();
    }
    else if (N == 10) {
        Option_10 ();
    }
    else if (N == 11) {
        Option_11 ();
    }
    else if (N == 12) {
        Option_12 (want_to_quit);
    }
}

void Use_Geometry_Tool () {
    bool want_to_quit = true;
    Initialisation ();
    while (want_to_quit != false) {
        double N;
        Initial_Choice (N);
        Choice (N, want_to_quit);
    }
}

int main()
{
    Use_Geometry_Tool ();
}