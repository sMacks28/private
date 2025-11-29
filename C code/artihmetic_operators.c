#include <stdio.h>

int main()
{

int x = 5;
float y = 2;

// + addition
// - subtraction
// * multiplication
// / division
// % modulus (gives 0 for fully divison and 1 for rest) e.g. 2 % 2 = 0 and 5 % 2 = 1
// ++ increment
// -- decrement

float z = x / y; //OR float z = x / (float) y; with keeping int y, pay attention to integer divison
printf("%f", z);
return 0;
}