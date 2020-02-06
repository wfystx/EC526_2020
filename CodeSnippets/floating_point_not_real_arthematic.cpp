/**** 

Floating point ONLY approximate real numbers.

For example the distributive law 

z * (x + y) ~ z * x + z* y

is only (a good) approximation.

Braket change the answer. Sometime  need
to figure out what is best.


  https://en.wikipedia.org/wiki/Single-precision_floating-point_format

 http://gcc.gnu.org/wiki/FloatingPointMath

***/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits.h> /* for CHAR_BIT */

using namespace std;

#define PI 3.141592653589793

/** formatted output of ieee-754 representation of float */
void show_ieee754 (float f)
{
    union {
        float f;
        uint32_t u;
    } fu = { .f = f };
    int i = sizeof f * CHAR_BIT;

    printf ("  ");
    while (i--)
        printf ("%d ", (fu.u >> i) & 0x1);

        putchar ('\n');
    printf (" |- - - - - - - - - - - - - - - - - - - - - - "
            "- - - - - - - - - -|\n");
    printf (" |s|      exp      |                  mantissa"
            "                   |\n\n");
}


int main()
{
  float x,y,z;
  float form1, form2, diff;
   
  x = 1.0/3.0;
  y = - 1.0/PI;
  z = sqrt(5);
 

 

  printf( "x =  %4.20f    y =  %4.20f  z  =  %4.20f \n\n",x,y,z); 

 form1 =  z * (x + y);
 printf(" z * (x + y)  =     %4.20f \n\n", form1);
   printf ("\nIEEE-754 Single-Precision representation of: %f\n\n", form1);
    show_ieee754 (form1);

 form2 =  (z * x)  + (z *y);
 printf("  (z * x)  + (z *y)  =   %4.20f \n\n", form2);
     printf ("\nIEEE-754 Single-Precision representation of: %f\n\n", form2);
    show_ieee754 (form2);

   diff = form1 - form2;
    printf(" ( z  +  (x + y))  -  ((z  + x) + y) =   %4.20f \n\n ", diff);
      printf ("\nIEEE-754 Single-Precision representation of: %f\n\n", diff);
    show_ieee754 (diff);
     
   
return 0;
}
