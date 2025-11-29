#include <stdio.h>

int main()
{

/* variable =  Allocated space in memory to store a value. We
refer to a variables name to acces the store value. That variable now
behaves as if it was the value it contains. BUT we need to declare what type
of data we are storing.
*/
int x;// declaration
x = 123; // initialization

int y = 321; //both

int age = 24; //integer
float gpa = 2.0; // floating point number
char grade = 'A'; //stores a single character
char name[] = "Max";// array of characters (string)

//printf("Hello I am %s\n", name);
//printf("and I am %d years old", age);
printf("Hello I am %s and I am %d years old\n", name, age);
//To output a variable inside of a printf statement you have to allocate space 
// via % followed by d (decimal) and ,variable after the ""
printf("I have an %c in mathematical methods\n", grade);
printf("My average GPA ist %f", gpa);
return 0;
}