#include <stdio.h>
#include <stdbool.h> // <----- include this for boolean variables
int main()
{

    char a = 'C';       //single character %c
    char b[] = "Bro";   //array of characters %s, use "" not ''

    float c = 3.141592; //4 bytes (32 bits of precision) 6-7 digits %f
    double d = 3.141592653589793; // 8 bytes (64 bits of precision) 15-16 digits %lf

    bool e = true; //1 byte (true or false)

    char f = 120; //1 byte (-128 to 127) %d or %c depending on which output you want
    //if you put %d for the format specifier in a printf statement, it gives just the decimal number as the output
    //but if you use %c it converts the number to a character according to the ASCII table
    
    printf("this is ASCII for 120: %c\n", f);
    printf("this is just the number according to the char f: %d\n", f);

    unsigned char g =255; //non negative number extend the range from 0 to 255 %d or %c
    //Remark: by making the prespecifier const before initilizing a variable, that variable cannot be changed

    const float PI = 3.14159; //good practice case: upper case letters for constants
    return 0;
}